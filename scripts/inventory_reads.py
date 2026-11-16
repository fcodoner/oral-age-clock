#!/usr/bin/env python3
# --------------------------------------------------------------------------
# inventory_reads.py
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
inventory_reads.py - look inside the downloaded FASTQ files and report, per run:
layout (single / paired / merged), read length, quality encoding, and which 16S
region was amplified. Nothing here is in the repository metadata reliably, and
all three decide how the run must be processed.

  python scripts/inventory_reads.py raw metadata/FINAL_16S_only.tsv \
         --out metadata/read_inventory --sample 2000

Outputs (in --out):
  read_inventory.tsv      one row per run: files, layout, read lengths, region
  cohort_summary.tsv      one row per cohort: layouts, regions, length ranges
  processing_groups.tsv   the groups to run DADA2 on (cohort x layout x region)
  warnings.txt            anything that needs a human decision

Standard library + pandas. Reads only the first --sample reads of each file.
"""
import argparse, collections, gzip, os, re, sys
import pandas as pd

# 16S primers, as sequences appear at the start of a read
PRIMERS = {
    "V1-V2": ["AGAGTTTGATCMTGGCTCAG", "AGAGTTTGATCCTGGCTCAG"],
    "V1-V3": ["AGAGTTTGATCMTGGCTCAG"],
    "V3-V4": ["CCTACGGGNGGCWGCAG", "CCTACGGGAGGCAGCAG", "CCTACGGGNBGCASCAG"],
    "V4":    ["GTGCCAGCMGCCGCGGTAA", "GTGCCAGCAGCCGCGGTAA"],
    "V4-V5": ["GTGYCAGCMGCCGCGGTAA"],
    "V6-V8": ["ACGCGHNRAACCTTACC"],
}
IUPAC = {"N": "ACGT", "W": "AT", "S": "GC", "M": "AC", "K": "GT", "R": "AG",
         "Y": "CT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
# typical amplicon length once merged, used to tell merged pairs from single-end reads
MERGED_LEN = {"V1-V2": 310, "V1-V3": 490, "V3-V4": 460, "V4": 253, "V4-V5": 375, "V6-V8": 500}


def matches(seq, primer, max_mismatch=2):
    if len(seq) < len(primer):
        return False
    mm = 0
    for s, p in zip(seq, primer):
        ok = (s == p) or (s in IUPAC.get(p, ""))
        if not ok:
            mm += 1
            if mm > max_mismatch:
                return False
    return True


def scan(path, n_reads):
    """Return (count, lengths, primer hits, max quality char) from the first n_reads."""
    lengths, hits, maxq = [], collections.Counter(), 0
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i >= n_reads * 4:
                    break
                k = i % 4
                if k == 1:
                    s = line.strip().upper()
                    lengths.append(len(s))
                    for region, plist in PRIMERS.items():
                        if any(matches(s, p) for p in plist):
                            hits[region] += 1
                            break
                elif k == 3:
                    q = line.strip()
                    if q:
                        maxq = max(maxq, max(ord(c) for c in q))
    except OSError as e:
        return None, [], hits, 0
    return len(lengths), lengths, hits, maxq


def med(xs):
    return int(pd.Series(xs).median()) if xs else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_dir")
    ap.add_argument("table", help="FINAL_16S_only.tsv or FINAL_sample_set_M1.tsv")
    ap.add_argument("--out", default="metadata/read_inventory")
    ap.add_argument("--sample", type=int, default=2000, help="reads to scan per file")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    t = pd.read_csv(a.table, sep="\t", low_memory=False)
    cohort = dict(zip(t.run_accession, t.cohort_id))
    fold = dict(zip(t.run_accession, t.fold_id)) if "fold_id" in t else {}

    files = collections.defaultdict(list)
    for f in sorted(os.listdir(a.raw_dir)):
        if not f.endswith((".fastq.gz", ".fq.gz")):
            continue
        m = re.match(r"([A-Z]{3}\d+)(_[12])?\.(fastq|fq)\.gz$", f)
        if m:
            files[m.group(1)].append(f)

    rows, warn = [], []
    missing = [r for r in t.run_accession if r not in files]
    if missing:
        warn.append(f"{len(missing)} runs in the table have no file on disk, "
                    f"e.g. {missing[:5]}")

    for i, (acc, fl) in enumerate(sorted(files.items()), 1):
        if i % 200 == 0:
            print(f"  scanned {i} runs...", file=sys.stderr)
        fl = sorted(fl)
        paired = any(f.endswith("_1.fastq.gz") for f in fl) and \
                 any(f.endswith("_2.fastq.gz") for f in fl)
        main_file = next((f for f in fl if f.endswith("_1.fastq.gz")), fl[0])
        n, lens, hits, maxq = scan(os.path.join(a.raw_dir, main_file), a.sample)
        if n is None:
            warn.append(f"{acc}: file could not be read (corrupt gzip?) - re-download")
            continue
        if n == 0:
            warn.append(f"{acc}: file is empty")
            continue

        region = hits.most_common(1)[0][0] if hits else ""
        frac = round(hits[region] / n, 3) if region else 0.0
        L = med(lens)
        # a single file whose reads are far longer than one sequencing read is a
        # pre-merged amplicon, not a single-end run; it must not go through DADA2
        # paired or single trimming blindly
        expected = MERGED_LEN.get(region, 0)
        if paired:
            layout = "paired"
        elif region and expected >= 350 and L >= expected * 0.85:
            # only a long amplicon (V3-V4, V1-V3, V4-V5, V6-V8) can be told apart
            # from a single read by length alone
            layout = "merged"
        elif L > 320:
            layout = "merged?"
        else:
            # a ~250 bp V4 file is ambiguous: merged V4 is also ~253 bp
            layout = "single"

        rows.append(dict(
            run_accession=acc, cohort_id=cohort.get(acc, ""), fold_id=fold.get(acc, ""),
            n_files=len(fl), files=";".join(fl), layout=layout,
            read_len_median=L, read_len_min=min(lens), read_len_max=max(lens),
            reads_scanned=n, region=region, region_primer_frac=frac,
            phred64="yes" if maxq > 74 else "no"))

    if not rows:
        sys.exit("No FASTQ files found. Check the raw folder path.")
    inv = pd.DataFrame(rows)
    inv.to_csv(os.path.join(a.out, "read_inventory.tsv"), sep="\t", index=False)

    def joinvc(s):
        return "; ".join(f"{k}:{v}" for k, v in s.value_counts().items())

    summ = inv.groupby("cohort_id").agg(
        n_runs=("run_accession", "size"), layouts=("layout", joinvc),
        regions=("region", lambda s: joinvc(s.replace("", "unknown"))),
        len_median=("read_len_median", "median"),
        len_range=("read_len_median", lambda s: f"{int(s.min())}-{int(s.max())}"),
        phred64=("phred64", joinvc))
    summ.to_csv(os.path.join(a.out, "cohort_summary.tsv"), sep="\t")

    inv["group"] = (inv.cohort_id.astype(str) + "__" + inv.layout + "__" +
                    inv.region.replace("", "UNKNOWN"))
    grp = inv.groupby("group").agg(n_runs=("run_accession", "size"),
                                   cohort_id=("cohort_id", "first"),
                                   layout=("layout", "first"), region=("region", "first"),
                                   len_median=("read_len_median", "median"))
    grp.sort_values("n_runs", ascending=False).to_csv(
        os.path.join(a.out, "processing_groups.tsv"), sep="\t")

    # things a human must decide
    for cid, g in inv.groupby("cohort_id"):
        if g.layout.nunique() > 1:
            warn.append(f"{cid}: mixed layouts ({joinvc(g.layout)}) - process separately")
        if g.region.replace("", "unknown").nunique() > 1:
            warn.append(f"{cid}: more than one 16S region detected "
                        f"({joinvc(g.region.replace('', 'unknown'))}) - check the source paper")
    amb = inv[(inv.layout == "single") & (inv.region.isin(["V4"])) &
              (inv.read_len_median >= 240)]
    if len(amb):
        warn.append(f"{len(amb)} V4 runs of about 250 bp in a single file are ambiguous: "
                    f"merged V4 amplicons are the same length as one read. Check the source "
                    f"paper for these cohorts: "
                    f"{', '.join(sorted(amb.cohort_id.unique())[:10])}")
    nohit = inv[inv.region == ""]
    if len(nohit):
        warn.append(f"{len(nohit)} runs have no detectable primer (already trimmed, or an "
                    f"unlisted primer). Cohorts: "
                    f"{', '.join(sorted(nohit.cohort_id.unique())[:10])}")
    weak = inv[(inv.region != "") & (inv.region_primer_frac < 0.5)]
    if len(weak):
        warn.append(f"{len(weak)} runs match a primer in under half their reads - "
                    f"verify the region before trusting it")
    if (inv.phred64 == "yes").any():
        warn.append("Some runs look like old Phred+64 quality encoding - convert before DADA2")
    open(os.path.join(a.out, "warnings.txt"), "w").write("\n".join(warn) + "\n")

    print(f"\n{len(inv)} runs scanned across {inv.cohort_id.nunique()} cohorts")
    print("\nLayouts:  " + joinvc(inv.layout))
    print("Regions:  " + joinvc(inv.region.replace("", "unknown")))
    print(f"\nProcessing groups: {inv.group.nunique()}")
    print(f"\n{len(warn)} warning(s) written to {a.out}/warnings.txt")
    for w in warn[:12]:
        print("  - " + w)


if __name__ == "__main__":
    main()
