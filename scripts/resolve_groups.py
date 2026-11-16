#!/usr/bin/env python3
# --------------------------------------------------------------------------
# resolve_groups.py
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
resolve_groups.py - decide the final region and layout for every run, combining
three sources of evidence, and write the processing groups.

  python scripts/resolve_groups.py metadata/read_inventory/read_inventory.tsv \
         --studies metadata/studies.tsv --out metadata/read_inventory

Evidence, in order of trust:
  1. primer found in the reads (strongest, but absent when the submitter trimmed them)
  2. the region declared by the source paper, from studies.tsv
  3. read length plus sequencing platform (weakest, used only to fill gaps)

Layout is corrected too: a 350-450 bp read from a 454 machine is a long single
read, not a merged pair, which the length test alone cannot tell.
"""
import argparse, os, sys
import pandas as pd

# typical single-read lengths, and merged amplicon lengths, per region
MERGED = {"V1-V2": 310, "V1-V3": 490, "V3-V4": 460, "V4": 253, "V4-V5": 375, "V6-V8": 500}
LONG_READ_PLATFORMS = ("454", "ion torrent", "pacbio", "nanopore", "minion", "gridion")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inventory")
    ap.add_argument("--studies", default="metadata/studies.tsv")
    ap.add_argument("--table", help="FINAL table, used to pick up instrument_model")
    ap.add_argument("--out", default="metadata/read_inventory")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="primer hit fraction needed to trust the detected region")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    inv = pd.read_csv(a.inventory, sep="\t", low_memory=False)
    inv["region"] = inv.region.fillna("")

    # instrument, if we can get it
    inv["instrument"] = ""
    if a.table and os.path.exists(a.table):
        t = pd.read_csv(a.table, sep="\t", low_memory=False)
        if "instrument_model" in t:
            inv["instrument"] = inv.run_accession.map(
                dict(zip(t.run_accession, t.instrument_model))).fillna("")

    # region declared by the source paper
    declared = {}
    if os.path.exists(a.studies):
        st = pd.read_csv(a.studies, sep="\t", low_memory=False)
        if "region_16S" in st and "cohort_id" in st:
            for _, r in st.iterrows():
                v = str(r.get("region_16S", "")).strip()
                if v and v.lower() not in ("confirm", "nan", ""):
                    declared[r["cohort_id"]] = v.split("/")[0]

    # per cohort, the region found in reads where the primer WAS detected
    votes = (inv[(inv.region != "") & (inv.region_primer_frac >= a.min_frac)]
             .groupby("cohort_id").region.agg(lambda s: s.value_counts().idxmax()))

    rows = []
    for _, r in inv.iterrows():
        cid = r.cohort_id
        if r.region and r.region_primer_frac >= a.min_frac:
            region, src = r.region, "primer in reads"
        elif cid in votes:
            region, src = votes[cid], "primer in other runs of the same cohort"
        elif cid in declared:
            region, src = declared[cid], "declared in the source paper"
        else:
            region, src = "UNKNOWN", "not resolved"

        # layout: long-read platforms give long single reads, not merged pairs
        layout = r.layout
        inst = str(r.instrument).lower()
        if layout.startswith("merged") and any(p in inst for p in LONG_READ_PLATFORMS):
            layout, lsrc = "single_long", "long-read platform"
        elif layout == "merged?":
            exp = MERGED.get(region, 0)
            if exp >= 350 and r.read_len_median >= exp * 0.85:
                layout, lsrc = "merged", "length matches a merged amplicon"
            else:
                layout, lsrc = "merged", "length too long for one short read"
        else:
            lsrc = "file count / length"

        rows.append(dict(run_accession=r.run_accession, cohort_id=cid,
                         fold_id=r.get("fold_id", ""), layout=layout, layout_evidence=lsrc,
                         region=region, region_evidence=src,
                         read_len_median=r.read_len_median, n_files=r.n_files,
                         files=r.files))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.out, "runs_resolved.tsv"), sep="\t", index=False)

    out["group"] = out.cohort_id + "__" + out.layout + "__" + out.region
    grp = out.groupby("group").agg(n_runs=("run_accession", "size"),
                                   cohort_id=("cohort_id", "first"),
                                   layout=("layout", "first"), region=("region", "first"),
                                   len_median=("read_len_median", "median"),
                                   region_evidence=("region_evidence", "first"))
    grp = grp.sort_values(["region", "layout", "n_runs"], ascending=[True, True, False])
    grp.to_csv(os.path.join(a.out, "processing_groups_resolved.tsv"), sep="\t")

    print(f"{len(out)} runs\n")
    print("Region after resolving:")
    print(out.region.value_counts().to_string())
    print("\nHow each region was decided:")
    print(out.region_evidence.value_counts().to_string())
    print("\nLayout:")
    print(out.layout.value_counts().to_string())
    print(f"\nProcessing groups: {out.group.nunique()}")
    un = out[out.region == "UNKNOWN"]
    if len(un):
        print(f"\nStill unresolved: {len(un)} runs in "
              f"{un.cohort_id.nunique()} cohorts -> "
              f"{', '.join(sorted(un.cohort_id.unique())[:12])}")
        print("For these, read the source paper's methods and put the region in the "
              "region_16S column of studies.tsv, then run this again.")


if __name__ == "__main__":
    main()
