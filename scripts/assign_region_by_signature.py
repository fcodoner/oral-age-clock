#!/usr/bin/env python3
# --------------------------------------------------------------------------
# assign_region_by_signature.py
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
assign_region_by_signature.py - give a 16S region to runs whose primers were
trimmed before deposit, by matching where the read starts against runs where the
primer was found.

Idea: a trimmed read begins at the first base after the primer. That position is
region-specific and highly conserved, so the first ~30 bases of a trimmed V4 read
look like the first ~30 bases after the V4 primer in an untrimmed run. We learn
those signatures from the resolved runs and match the unresolved ones to them.

  python scripts/assign_region_by_signature.py raw \
         metadata/read_inventory/runs_resolved.tsv \
         --out metadata/read_inventory --sample 400

Writes runs_resolved_v2.tsv (same columns, regions filled in where confident)
and region_signature_report.tsv showing the match score for every unknown run.
"""
import argparse, collections, gzip, os, sys
import pandas as pd

PRIMERS = {
    "V1-V2": ["AGAGTTTGATCMTGGCTCAG", "AGAGTTTGATCCTGGCTCAG"],
    "V1-V3": ["AGAGTTTGATCMTGGCTCAG"],
    "V3-V4": ["CCTACGGGNGGCWGCAG", "CCTACGGGAGGCAGCAG", "CCTACGGGNBGCASCAG"],
    "V4":    ["GTGCCAGCMGCCGCGGTAA", "GTGCCAGCAGCCGCGGTAA"],
    "V4-V5": ["GTGYCAGCMGCCGCGGTAA"],
}
IUPAC = {"N": "ACGT", "W": "AT", "S": "GC", "M": "AC", "K": "GT", "R": "AG",
         "Y": "CT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
K = 10          # k-mer length taken from the start of the read
TOPN = 400      # k-mers kept per region signature


def strip_primer(seq):
    for region, plist in PRIMERS.items():
        for p in plist:
            if len(seq) < len(p):
                continue
            mm = 0
            ok = True
            for s, c in zip(seq, p):
                if not (s == c or s in IUPAC.get(c, "")):
                    mm += 1
                    if mm > 2:
                        ok = False
                        break
            if ok:
                return seq[len(p):], region
    return seq, None


def read_starts(path, n_reads):
    """First K bases after any primer, for the first n_reads of a file."""
    out = []
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i >= n_reads * 4:
                    break
                if i % 4 == 1:
                    s, _ = strip_primer(line.strip().upper())
                    if len(s) >= K:
                        out.append(s[:K])
    except OSError:
        pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_dir")
    ap.add_argument("resolved")
    ap.add_argument("--out", default="metadata/read_inventory")
    ap.add_argument("--sample", type=int, default=400)
    ap.add_argument("--min-score", type=float, default=0.35)
    ap.add_argument("--min-margin", type=float, default=1.5,
                    help="best region must beat the runner-up by this factor")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    r = pd.read_csv(a.resolved, sep="\t", low_memory=False)
    r["region"] = r.region.fillna("UNKNOWN").astype(str).str.replace(
        r"\s*\(confirm\)", "", regex=True).str.strip()
    known = r[(r.region != "UNKNOWN") & (r.region != "")]
    unknown = r[(r.region == "UNKNOWN") | (r.region == "")]
    if unknown.empty:
        print("Nothing to resolve.")
        return
    print(f"learning signatures from {len(known)} resolved runs, "
          f"matching {len(unknown)} unresolved runs", file=sys.stderr)

    def first_file(files):
        f = str(files).split(";")[0]
        return os.path.join(a.raw_dir, f)

    # build one signature per region, from up to 25 runs each
    sig = {}
    for region, g in known.groupby("region"):
        c = collections.Counter()
        for _, row in g.head(25).iterrows():
            p = first_file(row.files)
            if os.path.exists(p):
                c.update(read_starts(p, a.sample))
        if c:
            sig[region] = set(k for k, _ in c.most_common(TOPN))
            print(f"  {region}: {len(sig[region])} signature k-mers", file=sys.stderr)

    rows, assigned = [], 0
    for i, (_, row) in enumerate(unknown.iterrows(), 1):
        if i % 100 == 0:
            print(f"  matched {i} runs...", file=sys.stderr)
        p = first_file(row.files)
        starts = read_starts(p, a.sample) if os.path.exists(p) else []
        rec = {"run_accession": row.run_accession, "cohort_id": row.cohort_id,
               "read_len_median": row.read_len_median}
        if not starts:
            rec.update(best="", score=0.0, margin=0.0, decision="no reads")
            rows.append(rec)
            continue
        scores = {reg: sum(1 for s in starts if s in ks) / len(starts)
                  for reg, ks in sig.items()}
        ranked = sorted(scores.items(), key=lambda kv: -kv[1])
        best, bs = ranked[0]
        second = ranked[1][1] if len(ranked) > 1 else 0.0
        margin = bs / second if second > 0 else (99.0 if bs > 0 else 0.0)
        ok = bs >= a.min_score and margin >= a.min_margin
        rec.update(best=best, score=round(bs, 3), margin=round(margin, 2),
                   decision="assigned" if ok else "too weak",
                   **{f"score_{k}": round(v, 3) for k, v in scores.items()})
        if ok:
            assigned += 1
        rows.append(rec)

    rep = pd.DataFrame(rows)
    rep.to_csv(os.path.join(a.out, "region_signature_report.tsv"), sep="\t", index=False)

    # a cohort is assigned as a whole, by majority of its confidently matched runs
    out = r.copy()
    for cid, g in rep.groupby("cohort_id"):
        good = g[g.decision == "assigned"]
        if len(good) >= max(3, 0.5 * len(g)):
            region = good.best.value_counts().idxmax()
            frac = (good.best == region).mean()
            if frac >= 0.8:
                out.loc[out.cohort_id == cid, "region"] = region
                out.loc[out.cohort_id == cid, "region_evidence"] = \
                    f"read-start signature ({frac:.0%} of {len(good)} runs)"
    out.to_csv(os.path.join(a.out, "runs_resolved_v2.tsv"), sep="\t", index=False)

    print(f"\n{assigned} of {len(unknown)} runs matched a signature confidently")
    print("\nRegion after this pass:")
    print(out.region.value_counts().to_string())
    still = out[out.region.isin(["UNKNOWN", ""])]
    if len(still):
        print(f"\nStill unresolved: {len(still)} runs in {still.cohort_id.nunique()} "
              f"cohorts -> {', '.join(sorted(still.cohort_id.unique()))}")
        print("For these, read the source paper and set region_16S in studies.tsv.")


if __name__ == "__main__":
    main()
