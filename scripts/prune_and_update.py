#!/usr/bin/env python3
# --------------------------------------------------------------------------
# prune_and_update.py
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
prune_and_update.py - remove excluded cohorts from the dataset and refresh every
downstream file, in one step.

  python scripts/prune_and_update.py \
         --drop CHINA_48473 MEXICO_96341 \
         --raw raw --metadata metadata

What it does, in order:
  1. works out which runs belong to the dropped cohorts
  2. moves their FASTQ files to raw/_excluded/ (use --delete to erase instead)
  3. rewrites the sample table and accession lists without them
  4. fills in regions from studies.tsv for any run still marked UNKNOWN
  5. regenerates the processing groups
  6. writes a report of everything it changed

Files are MOVED, not deleted, unless you ask. Every file it rewrites is backed
up first with a .bak extension, so a mistake is recoverable.
"""
import argparse, os, shutil, sys
import pandas as pd

SAMPLE_TABLE = "FINAL_sample_set_M1.tsv"
ACC_16S = "FINAL_accessions_16S.txt"
ACC_WGS = "FINAL_accessions_shotgun.txt"
RESOLVED = os.path.join("read_inventory", "runs_resolved_v2.tsv")
GROUPS = os.path.join("read_inventory", "processing_groups_final.tsv")


def backup(path):
    if os.path.exists(path):
        shutil.copy(path, path + ".bak")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drop", nargs="*", default=[], help="cohort_ids to remove")
    ap.add_argument("--drop-from-studies", action="store_true",
                    help="also drop every cohort marked include=N in studies.tsv")
    ap.add_argument("--raw", default="raw")
    ap.add_argument("--metadata", default="metadata")
    ap.add_argument("--studies", default=None, help="defaults to <metadata>/studies.tsv")
    ap.add_argument("--delete", action="store_true", help="erase files instead of moving them")
    ap.add_argument("--dry-run", action="store_true", help="report only, change nothing")
    ap.add_argument("--sweep-orphans", action="store_true",
                    help="also move any FASTQ in raw/ that is not in the final sample table "
                         "(left over from an earlier version of the list)")
    a = ap.parse_args()

    md = a.metadata
    studies_path = a.studies or os.path.join(md, "studies.tsv")
    table_path = os.path.join(md, SAMPLE_TABLE)
    resolved_path = os.path.join(md, RESOLVED)
    for p in (table_path, studies_path):
        if not os.path.exists(p):
            sys.exit(f"Cannot find {p}")

    t = pd.read_csv(table_path, sep="\t", low_memory=False)
    st = pd.read_csv(studies_path, sep="\t", low_memory=False)
    st["region_16S"] = (st.get("region_16S", "").astype(str)
                        .str.replace(r"\s*\(confirm\)", "", regex=True).str.strip()
                        .replace({"nan": "", "confirm": ""}))

    drop = set(a.drop)
    if a.drop_from_studies:
        drop |= set(st.loc[st.include.astype(str).str.upper().str.startswith("N"), "cohort_id"])
    unknown_cohorts = drop - set(t.cohort_id)
    if unknown_cohorts:
        print(f"Warning: not in the sample table, ignored: {sorted(unknown_cohorts)}")
    drop &= set(t.cohort_id)

    log = [f"Cohorts to remove: {sorted(drop) or 'none'}"]
    gone = t[t.cohort_id.isin(drop)]
    kept = t[~t.cohort_id.isin(drop)]
    log.append(f"Participants: {len(t)} -> {len(kept)}  (removing {len(gone)})")
    log.append(f"Cohorts:      {t.cohort_id.nunique()} -> {kept.cohort_id.nunique()}")
    if "fold_id" in kept:
        log.append(f"Folds:        {t.fold_id.nunique()} -> {kept.fold_id.nunique()}")

    # 2. the FASTQ files of the removed runs
    moved = missing = 0
    dest = os.path.join(a.raw, "_excluded")
    if len(gone):
        if not a.dry_run and not a.delete:
            os.makedirs(dest, exist_ok=True)
        for _, r in gone.iterrows():
            for suffix in ("", "_1", "_2"):
                f = os.path.join(a.raw, f"{r.run_accession}{suffix}.fastq.gz")
                if os.path.exists(f):
                    moved += 1
                    if not a.dry_run:
                        if a.delete:
                            os.remove(f)
                        else:
                            shutil.move(f, os.path.join(dest, os.path.basename(f)))
                elif suffix == "":
                    missing += 1
    log.append(f"FASTQ files {'deleted' if a.delete else 'moved to ' + dest}: {moved}")

    # 2b. files on disk that belong to no run in the final table
    if a.sweep_orphans:
        keep_acc = set(kept.run_accession)
        orphan_dest = os.path.join(a.raw, "_not_in_final_list")
        orphans = []
        for f in sorted(os.listdir(a.raw)):
            if not f.endswith((".fastq.gz", ".fq.gz")):
                continue
            acc = f.split(".")[0]
            for suf in ("_1", "_2"):
                if acc.endswith(suf):
                    acc = acc[:-len(suf)]
            if acc not in keep_acc:
                orphans.append(f)
        if orphans and not a.dry_run:
            os.makedirs(orphan_dest, exist_ok=True)
            for f in orphans:
                shutil.move(os.path.join(a.raw, f), os.path.join(orphan_dest, f))
        log.append(f"Files on disk not in the final list: {len(orphans)}"
                   + (f" -> {orphan_dest}" if orphans else ""))
        if orphans:
            log.append("  e.g. " + ", ".join(orphans[:6]))

    # 2c. runs in the final table with no file on disk
    have = set()
    for f in os.listdir(a.raw):
        if f.endswith((".fastq.gz", ".fq.gz")):
            acc = f.split(".")[0]
            for suf in ("_1", "_2"):
                if acc.endswith(suf):
                    acc = acc[:-len(suf)]
            have.add(acc)
    absent = sorted(set(kept[kept.library_strategy != "WGS"].run_accession) - have)
    if absent:
        log.append(f"WARNING: {len(absent)} runs in the final list have NO file on disk, "
                   f"e.g. {absent[:5]}")
        if not a.dry_run:
            open(os.path.join(md, "runs_missing_files.txt"), "w").write("\n".join(absent) + "\n")
            log.append(f"  listed in {os.path.join(md, 'runs_missing_files.txt')} - "
                       f"re-run the downloader on that file")
    else:
        log.append("Every run in the final list has a file on disk.")

    # 3. rewrite the sample table and accession lists
    if not a.dry_run:
        backup(table_path)
        kept.to_csv(table_path, sep="\t", index=False)
        for name, sub in ((ACC_16S, kept[kept.library_strategy != "WGS"]),
                          (ACC_WGS, kept[kept.library_strategy == "WGS"])):
            p = os.path.join(md, name)
            backup(p)
            open(p, "w").write("\n".join(sub.run_accession) + "\n")
        # the 16S-only table, if it exists
        p16 = os.path.join(md, "FINAL_16S_only.tsv")
        if os.path.exists(p16):
            backup(p16)
            k16 = kept[kept.library_strategy != "WGS"]
            k16.to_csv(p16, sep="\t", index=False)
    log.append(f"Rewritten: {SAMPLE_TABLE}, {ACC_16S}, {ACC_WGS}, FINAL_16S_only.tsv")

    # 4 and 5. regions and processing groups
    if os.path.exists(resolved_path):
        rv = pd.read_csv(resolved_path, sep="\t", low_memory=False)
        rv["region"] = (rv.region.astype(str)
                        .str.replace(r"\s*\(confirm\)", "", regex=True).str.strip())
        rv = rv[~rv.cohort_id.isin(drop)]
        declared = dict(zip(st.cohort_id, st.region_16S))
        filled = 0
        for cid, reg in declared.items():
            if not reg:
                continue
            m = (rv.cohort_id == cid) & (rv.region.isin(["UNKNOWN", "", "nan"]))
            if m.any():
                rv.loc[m, "region"] = reg
                rv.loc[m, "region_evidence"] = "declared in studies.tsv (from the source paper)"
                filled += int(m.sum())
        log.append(f"Regions filled in from studies.tsv: {filled} runs")
        still = rv[rv.region.isin(["UNKNOWN", "", "nan"])]
        if len(still):
            log.append(f"Still without a region: {len(still)} runs in "
                       f"{still.cohort_id.nunique()} cohorts -> "
                       f"{', '.join(sorted(still.cohort_id.unique()))}")
        else:
            log.append("Every remaining run now has a region.")

        rv["group"] = rv.cohort_id + "__" + rv.layout + "__" + rv.region
        grp = rv.groupby("group").agg(
            n_runs=("run_accession", "size"), cohort_id=("cohort_id", "first"),
            layout=("layout", "first"), region=("region", "first"),
            len_median=("read_len_median", "median"),
            region_evidence=("region_evidence", "first")).sort_values(
            ["region", "layout", "n_runs"], ascending=[True, True, False])
        if not a.dry_run:
            backup(resolved_path)
            rv.to_csv(resolved_path, sep="\t", index=False)
            grp.to_csv(os.path.join(md, GROUPS), sep="\t")
        log.append(f"Processing groups: {rv.group.nunique()} "
                   f"-> {os.path.join(md, GROUPS)}")
        log.append("\nRuns per region:\n" + rv.region.value_counts().to_string())
        log.append("\nRuns per region and layout:\n" +
                   pd.crosstab(rv.region, rv.layout, margins=True).to_string())
    else:
        log.append(f"({resolved_path} not found - skipped the region and group update)")

    txt = "\n".join(log)
    if not a.dry_run:
        open(os.path.join(md, "prune_report.txt"), "w").write(txt + "\n")
    print(txt)
    if a.dry_run:
        print("\n(dry run - nothing was changed)")


if __name__ == "__main__":
    main()
