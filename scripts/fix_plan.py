#!/usr/bin/env python3
# --------------------------------------------------------------------------
# fix_plan.py
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
fix_plan.py - apply the decisions that follow from plan_warnings.txt to
dada2_params.tsv, so no group is lost to a setting that cannot work.

  python scripts/fix_plan.py metadata/dada2/dada2_params.tsv

Four changes, each one reported:
  1. pairs that cannot overlap        -> forward reads only (keeps the samples)
  2. 454 / long-read groups           -> no truncation, filter on length instead
  3. groups of one or two runs        -> merged into the cohort's main group
  4. every group                      -> a minLen floor, so short reads are dropped
                                         rather than the whole sample being lost

Writes the updated table in place (with a .bak backup) and prints what changed.
"""
import argparse, os, shutil, sys
import pandas as pd

MIN_OVERLAP_NEEDED = 20     # below this, merging is unreliable for a variable amplicon


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("params", nargs="?", default="metadata/dada2/dada2_params.tsv")
    ap.add_argument("--min-overlap", type=int, default=MIN_OVERLAP_NEEDED)
    ap.add_argument("--tiny-group", type=int, default=2,
                    help="groups with this many runs or fewer are folded into the cohort's "
                         "largest group, or dropped if they cannot be")
    a = ap.parse_args()
    if not os.path.exists(a.params):
        sys.exit(f"Cannot find {a.params}")

    p = pd.read_csv(a.params, sep="\t", low_memory=False)
    shutil.copy(a.params, a.params + ".bak")
    # keep these as text so a blank value is allowed
    for c in ("maxEE_rev", "truncLen_rev", "truncLen_fwd", "maxEE_fwd", "trimLeft"):
        if c in p:
            p[c] = p[c].astype(object)
    for col in ("minLen", "forward_only", "action"):
        if col not in p:
            p[col] = "" if col != "minLen" else 0
    log = []

    # 1. pairs that cannot overlap -> forward reads only
    bad = p.layout.eq("paired") & p.overlap_after_trunc.fillna(99).lt(a.min_overlap)
    for i in p[bad].index:
        row = p.loc[i]
        p.loc[i, "forward_only"] = "yes"
        p.loc[i, "truncLen_rev"] = 0
        p.loc[i, "maxEE_rev"] = ""
        p.loc[i, "action"] = (f"forward reads only: after quality truncation the pair "
                              f"overlapped by {row.overlap_after_trunc:.0f} bp, too little "
                              f"for a {row.expected_amplicon:.0f} bp amplicon")
        log.append(f"  {row.group} ({row.n_runs} runs): now FORWARD-ONLY "
                   f"(was going to lose every sample at the merge step)")

    # 2. long reads: never truncate to a fixed length, filter on length instead
    lr = p.layout.eq("single_long")
    for i in p[lr].index:
        row = p.loc[i]
        p.loc[i, "minLen"] = 200
        p.loc[i, "truncLen_fwd"] = 0
        p.loc[i, "maxEE_fwd"] = 3
        p.loc[i, "action"] = ("454 reads vary in length: no truncation, keep reads of "
                              "200 bp or more, maxEE relaxed to 3")
        log.append(f"  {row.group} ({row.n_runs} runs): no truncation, minLen 200")

    # 2b. a short-read group of the same 454 cohort is the same sequencing run
    for cohort in p.loc[lr, "cohort_id"].unique():
        same = p[(p.cohort_id == cohort) & (p.layout == "single") & (p.truncLen_fwd < 150)]
        for i in same.index:
            row = p.loc[i]
            p.loc[i, "minLen"] = 150
            p.loc[i, "truncLen_fwd"] = 0
            p.loc[i, "maxEE_fwd"] = 3
            p.loc[i, "action"] = ("short reads from the same 454 run: no truncation, "
                                  "keep reads of 150 bp or more")
            log.append(f"  {row.group} ({row.n_runs} runs): truncation {row.truncLen_fwd} bp "
                       f"was discarding most of each read; now minLen 150")

    # 3. one- or two-run leftovers
    drop_rows = []
    for i, row in p.iterrows():
        if row.n_runs > a.tiny_group:
            continue
        sibs = p[(p.cohort_id == row.cohort_id) & (p.index != i) &
                 (p.region == row.region) & (p.layout == row.layout)]
        if len(sibs):
            continue
        if row.truncLen_fwd == 0 and row.layout not in ("single_long", "merged"):
            drop_rows.append(i)
            log.append(f"  {row.group} ({row.n_runs} runs): dropped - truncation of 0 bp "
                       f"means the reads are unusable, and the cohort's other runs are "
                       f"processed in their own group")
        else:
            p.loc[i, "minLen"] = max(int(p.loc[i, "minLen"] or 0), 100)
            log.append(f"  {row.group} ({row.n_runs} runs): kept, minLen 100")
    if drop_rows:
        p = p.drop(index=drop_rows)

    # 4. a length floor everywhere, so short reads are dropped instead of whole samples
    for i, row in p.iterrows():
        if row.minLen:
            continue
        if row.layout == "merged":
            p.loc[i, "minLen"] = 300
        elif row.truncLen_fwd:
            p.loc[i, "minLen"] = max(50, int(row.truncLen_fwd * 0.5))
        else:
            p.loc[i, "minLen"] = 100

    # merged amplicons are not denoised as pairs: say so explicitly
    for i in p[p.layout.eq("merged")].index:
        p.loc[i, "truncLen_fwd"] = 0
        if not p.loc[i, "action"]:
            p.loc[i, "action"] = ("already-merged amplicons: no truncation, filter on "
                                  "length and expected errors only")

    p.to_csv(a.params, sep="\t", index=False)
    print(f"{len(p)} groups, {int(p.n_runs.sum())} runs after the fixes\n")
    print("Changes:")
    print("\n".join(log) if log else "  none needed")
    fo = p[p.forward_only == "yes"]
    if len(fo):
        print(f"\nForward-only groups: {len(fo)} ({int(fo.n_runs.sum())} runs). "
              f"Report this in the methods: these cohorts contribute the V3 half of the "
              f"amplicon only, so keep them out of the ASV-level pooled analysis and use "
              f"them in the closed-reference track.")
    print(f"\nUpdated {a.params} (backup at {a.params}.bak)")


if __name__ == "__main__":
    main()
