#!/usr/bin/env python3
# --------------------------------------------------------------------------
# plan_dada2.py
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
"""
plan_dada2.py - work out the DADA2 settings for each processing group by reading
the actual quality profiles, and write the plan DADA2 will follow.

  python scripts/plan_dada2.py raw metadata/read_inventory/runs_resolved_v2.tsv \
         --out metadata/dada2 --sample 3000

For every group (cohort x layout x region) it reports, from the reads themselves:
  - the position where quality drops below the threshold, forward and reverse
  - a suggested truncation length for each
  - whether the pair can still overlap enough to merge after truncation
  - how much of the primer to cut (0 when the submitter already trimmed it)

Outputs (in --out):
  dada2_params.tsv     one row per group: the settings run_dada2.R will use
  group_manifests/     one file list per group, for run_dada2.R to read
  plan_warnings.txt    groups that need a decision before processing
"""
import argparse, collections, gzip, os, re, sys
import pandas as pd

# length of the forward primer, cut from the start of read 1 when still present
PRIMER_LEN = {"V1-V2": 20, "V1-V3": 20, "V3-V4": 17, "V4": 19, "V4-V5": 19, "V6-V8": 17}
# approximate length of the full amplicon, used to check pairs can still overlap
AMPLICON = {"V1-V2": 310, "V1-V3": 490, "V3-V4": 465, "V4": 253, "V4-V5": 375, "V6-V8": 500}
MIN_OVERLAP = 12          # DADA2's default minimum overlap for merging
OVERLAP_MARGIN = 20       # extra bases of safety on top of the minimum


def quality_profile(path, n_reads):
    """Mean Phred quality at each cycle, and the read-length distribution."""
    totals, counts, lengths = collections.Counter(), collections.Counter(), []
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i >= n_reads * 4:
                    break
                if i % 4 == 3:
                    q = line.rstrip("\n")
                    lengths.append(len(q))
                    for pos, ch in enumerate(q):
                        totals[pos] += ord(ch) - 33
                        counts[pos] += 1
    except OSError:
        return {}, []
    return {p: totals[p] / counts[p] for p in sorted(counts)}, lengths


def trunc_at(profile, lengths, min_q, min_frac=0.8):
    """Last cycle where mean quality is still above min_q and most reads reach."""
    if not profile:
        return 0
    n = len(lengths)
    keep = 0
    for pos in sorted(profile):
        enough_reads = sum(1 for L in lengths if L > pos) >= min_frac * n
        if profile[pos] >= min_q and enough_reads:
            keep = pos + 1
        else:
            break
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_dir")
    ap.add_argument("resolved")
    ap.add_argument("--out", default="metadata/dada2")
    ap.add_argument("--sample", type=int, default=3000, help="reads to profile per file")
    ap.add_argument("--files-per-group", type=int, default=6)
    ap.add_argument("--min-q", type=float, default=25.0)
    ap.add_argument("--inventory", default=None,
                    help="read_inventory.tsv, used to see if primers are still present")
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, "group_manifests"), exist_ok=True)

    r = pd.read_csv(a.resolved, sep="\t", low_memory=False)
    r["region"] = r.region.astype(str).str.replace(r"\s*\(confirm\)", "", regex=True).str.strip()
    if "group" not in r:
        r["group"] = r.cohort_id + "__" + r.layout + "__" + r.region

    primer_frac = {}
    inv_path = a.inventory or os.path.join(os.path.dirname(a.resolved), "read_inventory.tsv")
    if os.path.exists(inv_path):
        inv = pd.read_csv(inv_path, sep="\t", low_memory=False)
        if "region_primer_frac" in inv:
            primer_frac = dict(zip(inv.run_accession, inv.region_primer_frac))

    rows, warn = [], []
    for gname, g in r.groupby("group"):
        cohort = g.cohort_id.iloc[0]
        layout = g.layout.iloc[0]
        region = g.region.iloc[0]

        # the files for this group
        man = []
        for _, row in g.iterrows():
            acc = row.run_accession
            f1 = os.path.join(a.raw_dir, f"{acc}_1.fastq.gz")
            f2 = os.path.join(a.raw_dir, f"{acc}_2.fastq.gz")
            f0 = os.path.join(a.raw_dir, f"{acc}.fastq.gz")
            if layout == "paired" and os.path.exists(f1) and os.path.exists(f2):
                man.append((acc, f1, f2))
            elif os.path.exists(f0):
                man.append((acc, f0, ""))
            elif os.path.exists(f1):
                man.append((acc, f1, ""))
        if not man:
            warn.append(f"{gname}: no files found on disk - skipped")
            continue

        mpath = os.path.join(a.out, "group_manifests", f"{gname}.tsv")
        pd.DataFrame(man, columns=["run_accession", "fwd", "rev"]).to_csv(
            mpath, sep="\t", index=False)

        # profile a handful of files
        probe = man[:a.files_per_group]
        pf, lf, pr_, lr = collections.Counter(), [], collections.Counter(), []
        pf_n = pr_n = 0
        for _, f1, f2 in probe:
            p, L = quality_profile(f1, a.sample)
            if p:
                pf_n += 1
                for k, v in p.items():
                    pf[k] += v
                lf.extend(L)
            if f2:
                p2, L2 = quality_profile(f2, a.sample)
                if p2:
                    pr_n += 1
                    for k, v in p2.items():
                        pr_[k] += v
                    lr.extend(L2)
        prof_f = {k: v / pf_n for k, v in pf.items()} if pf_n else {}
        prof_r = {k: v / pr_n for k, v in pr_.items()} if pr_n else {}

        tf = trunc_at(prof_f, lf, a.min_q)
        tr = trunc_at(prof_r, lr, a.min_q) if prof_r else 0

        # primer still present?
        fr = [primer_frac.get(x[0], 0) for x in probe if x[0] in primer_frac]
        has_primer = (sum(fr) / len(fr)) >= 0.5 if fr else False
        trim_left = PRIMER_LEN.get(region, 0) if has_primer else 0

        # can the pair still overlap after truncation?
        amp = AMPLICON.get(region, 0)
        overlap = (tf + tr - amp) if (layout == "paired" and amp) else None
        note = ""
        if layout == "paired" and amp:
            need = MIN_OVERLAP + OVERLAP_MARGIN
            if overlap < need:
                # give back length to the reverse read, which is usually the weaker one
                short = need - overlap
                tr_new = min(tr + short, int(pd.Series(lr).median()) if lr else tr)
                if tf + tr_new - amp >= need:
                    note = (f"reverse truncation raised from {tr} to {tr_new} to keep "
                            f"{MIN_OVERLAP}+ bases of overlap")
                    tr = tr_new
                    overlap = tf + tr - amp
                else:
                    note = "CANNOT MERGE: reads too short for this amplicon"
                    warn.append(f"{gname}: after quality truncation the pair cannot overlap "
                                f"({tf}+{tr} vs amplicon {amp}). Use forward reads only, "
                                f"or process with justConcatenate.")

        if layout in ("merged", "single_long"):
            note = (note + "; " if note else "") + \
                   ("454/long reads: use HOMOPOLYMER_GAP_PENALTY=-1 and BAND_SIZE=32"
                    if layout == "single_long" else
                    "already-merged amplicons: no denoising of pairs, closed-reference track only")

        rows.append(dict(
            group=gname, cohort_id=cohort, layout=layout, region=region, n_runs=len(man),
            manifest=mpath, trimLeft=trim_left, truncLen_fwd=tf, truncLen_rev=tr,
            maxEE_fwd=2, maxEE_rev=(2 if layout == "paired" else ""),
            truncQ=2, expected_amplicon=amp, overlap_after_trunc=overlap,
            primer_in_reads="yes" if has_primer else "no",
            mean_q_at_trunc_fwd=round(prof_f.get(max(tf - 1, 0), 0), 1),
            mean_q_at_trunc_rev=round(prof_r.get(max(tr - 1, 0), 0), 1) if prof_r else "",
            note=note))

        if tf < 100:
            warn.append(f"{gname}: forward truncation only {tf} bp - quality drops very early, "
                        f"inspect this group before trusting it")

    if not rows:
        sys.exit("No groups could be planned - check the raw folder and the manifest.")
    p = pd.DataFrame(rows).sort_values(["region", "layout", "n_runs"],
                                       ascending=[True, True, False])
    p.to_csv(os.path.join(a.out, "dada2_params.tsv"), sep="\t", index=False)
    open(os.path.join(a.out, "plan_warnings.txt"), "w").write("\n".join(warn) + "\n")

    print(f"{len(p)} groups planned, {p.n_runs.sum()} runs\n")
    print(p[["group", "n_runs", "layout", "region", "trimLeft", "truncLen_fwd",
             "truncLen_rev", "overlap_after_trunc"]].to_string(index=False))
    print(f"\n{len(warn)} warning(s) -> {a.out}/plan_warnings.txt")
    for w in warn[:10]:
        print("  - " + w)


if __name__ == "__main__":
    main()
