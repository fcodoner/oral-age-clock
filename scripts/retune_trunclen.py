#!/usr/bin/env python3
# --------------------------------------------------------------------------
# retune_trunclen.py
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
retune_trunclen.py - choose a truncLen that does not throw away whole samples.

DADA2 discards any read SHORTER than truncLen. If a cohort's reads vary in
length, a truncLen at the median silently deletes every sample on the short
side. This measures the real length distribution of a group and suggests a
value that keeps the samples.

  python scripts/retune_trunclen.py AGP_ORAL__single__V4 LEE17_CN__single__V4 \
         --params metadata/dada2/dada2_params.tsv --apply

Without --apply it only reports. With --apply it writes the new truncLen and a
matching minLen into the parameter table (backup kept).
"""
import argparse, gzip, os, shutil, sys
import pandas as pd


def read_lengths(path, n_reads):
    out = []
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i >= n_reads * 4:
                    break
                if i % 4 == 1:
                    out.append(len(line.strip()))
    except OSError:
        pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("groups", nargs="+")
    ap.add_argument("--params", default="metadata/dada2/dada2_params.tsv")
    ap.add_argument("--sample", type=int, default=1000, help="reads per file")
    ap.add_argument("--files", type=int, default=40, help="files to measure per group")
    ap.add_argument("--keep", type=float, default=0.95,
                    help="fraction of SAMPLES that must keep most of their reads")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    p = pd.read_csv(a.params, sep="\t", low_memory=False)
    for c in ("action", "note"):
        if c in p:
            p[c] = p[c].astype(object).fillna("")
    changed = False
    for gname in a.groups:
        row = p[p.group == gname]
        if row.empty:
            print(f"{gname}: not in the parameter table")
            continue
        row = row.iloc[0]
        man = pd.read_csv(row.manifest, sep="\t")
        files = man.fwd.tolist()
        step = max(1, len(files) // a.files)
        probe = files[::step][:a.files]

        per_sample_p10, all_lens = [], []
        for f in probe:
            if not os.path.exists(f):
                continue
            L = read_lengths(f, a.sample)
            if L:
                s = pd.Series(L)
                per_sample_p10.append(s.quantile(0.10))
                all_lens.extend(L)
        if not all_lens:
            print(f"{gname}: no readable files")
            continue
        s = pd.Series(all_lens)
        ps = pd.Series(per_sample_p10)

        # a sample survives only if its reads reach truncLen, so the limit is set by
        # the SHORTEST samples, not by the overall read-length distribution
        suggested = int(ps.quantile(1 - a.keep))
        suggested = max(100, min(suggested, int(s.quantile(0.5))))
        cur = int(row.truncLen_fwd) if row.truncLen_fwd else 0

        print(f"\n{gname}  ({row.n_runs} runs, {len(probe)} files measured)")
        print(f"  read length overall: min {s.min()}, p10 {int(s.quantile(.1))}, "
              f"median {int(s.median())}, max {s.max()}")
        print(f"  per-sample 10th percentile: min {int(ps.min())}, "
              f"median {int(ps.median())}, max {int(ps.max())}")
        print(f"  current truncLen {cur}  ->  suggested {suggested}")
        would_die = (ps < cur).mean() if cur else 0
        print(f"  samples whose reads are mostly shorter than the current truncLen: "
              f"{would_die:.0%}")
        if a.apply:
            if not changed:
                shutil.copy(a.params, a.params + ".bak")
                changed = True
            # DADA2 applies trimLeft AFTER truncation, so a read ends up
            # truncLen - trimLeft bases long. minLen is checked against THAT, so a
            # minLen equal to truncLen rejects every read when trimLeft > 0.
            trim_left = int(row.trimLeft) if not pd.isna(row.trimLeft) else 0
            safe_minlen = max(50, suggested - trim_left - 5)
            p.loc[p.group == gname, "truncLen_fwd"] = suggested
            p.loc[p.group == gname, "minLen"] = safe_minlen
            if trim_left:
                print(f"  trimLeft is {trim_left}, so reads end up {suggested - trim_left} bp: "
                      f"minLen set to {safe_minlen}, not {suggested}")
            p.loc[p.group == gname, "action"] = (
                f"truncLen lowered from {cur} to {suggested} (minLen {safe_minlen}): "
                f"the original value discarded every read in {would_die:.0%} of samples")
            print(f"  applied")

    if changed:
        p.to_csv(a.params, sep="\t", index=False)
        print(f"\nUpdated {a.params} (backup at {a.params}.bak)")
        print("Delete the old output of these groups and run them again:")
        for g in a.groups:
            print(f"  rm processed/{g}_*")
    else:
        print("\n(report only - add --apply to write the new values)")


if __name__ == "__main__":
    main()
