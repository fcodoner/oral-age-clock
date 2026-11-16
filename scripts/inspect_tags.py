#!/usr/bin/env python3
# --------------------------------------------------------------------------
# inspect_tags.py
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
inspect_tags.py - find out which metadata labels a cohort actually uses.
Works from any manifest_runs.tsv (the all_attributes column), so it does not
need attribute_census.tsv and does not download anything.

  # 1. which cohorts have no sample type, and which tags they use instead
  python scripts/inspect_tags.py metadata/harvest_<date>/manifest_runs.tsv --missing site

  # 2. everything one cohort records, with example values
  python scripts/inspect_tags.py metadata/harvest_<date>/manifest_runs.tsv --cohort SPAIN_11189

  # 3. candidate tags across all cohorts for a field
  python scripts/inspect_tags.py metadata/harvest_<date>/manifest_runs.tsv --find sex
"""
import argparse, collections, sys
import pandas as pd

PATTERNS = {
    "site": ["body", "site", "material", "source", "specimen", "medium", "tissue",
             "matrix", "habitat", "sample_type", "env_"],
    "sex":  ["sex", "gender"],
    "country": ["country", "geo", "location", "nation", "region"],
    "age":  ["age", "birth", "year_of"],
    "health": ["disease", "health", "diagnos", "status", "condition", "group", "case"],
    "abx": ["antibiot", "medicat", "drug"],
}


def attrs_of(s):
    d = {}
    if isinstance(s, str):
        for kv in s.split(" | "):
            k, _, v = kv.partition("=")
            if k.strip():
                d[k.strip()] = v.strip()
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--cohort", help="show every tag this cohort uses")
    ap.add_argument("--find", choices=list(PATTERNS), help="candidate tags for a field")
    ap.add_argument("--missing", choices=list(PATTERNS),
                    help="cohorts where our harmonised field is empty, and what they use instead")
    a = ap.parse_args()

    m = pd.read_csv(a.manifest, sep="\t", low_memory=False)
    m["TAGS"] = m.all_attributes.map(attrs_of)

    if a.cohort:
        g = m[m.cohort_id == a.cohort]
        if g.empty:
            sys.exit(f"No rows for cohort {a.cohort}. Cohorts: "
                     + ", ".join(sorted(m.cohort_id.unique())[:30]))
        print(f"{a.cohort}: {len(g)} runs, BioProject {g.bioproject.iloc[0]}\n")
        tally = collections.Counter()
        example = {}
        for d in g["TAGS"]:
            for k, v in d.items():
                tally[k] += 1
                if v and k not in example:
                    example[k] = v
        for k, n in tally.most_common():
            print(f"  {k:<38} {n:>6}   e.g. {str(example.get(k,''))[:45]}")
        return

    if a.find:
        pats = PATTERNS[a.find]
        tally = collections.Counter()
        example = {}
        for d in m["TAGS"]:
            for k, v in d.items():
                if any(p in k.lower() for p in pats):
                    tally[k] += 1
                    if v and k not in example:
                        example[k] = v
        print(f"Tags that might hold '{a.find}':\n")
        for k, n in tally.most_common(40):
            print(f"  {k:<38} {n:>7} samples   e.g. {str(example.get(k,''))[:40]}")
        print("\nAdd the useful ones to the matching *_TAGS list in "
              "harvest_oral_age_metadata.py, then re-run finalize_dataset.sh")
        return

    if a.missing:
        field = {"site": "sample_material", "sex": "sex", "country": "country",
                 "age": "age_years", "health": "health_status_raw",
                 "abx": "antibiotics"}[a.missing]
        empty = m[m[field].fillna("").astype(str).str.strip().eq("")]
        if empty.empty:
            print(f"No cohort is missing '{field}'.")
            return
        print(f"{len(empty)} runs have no '{field}'. By cohort:\n")
        pats = PATTERNS[a.missing]
        for cid, g in sorted(empty.groupby("cohort_id"), key=lambda kv: -len(kv[1]))[:20]:
            cand = collections.Counter()
            example = {}
            for d in g["TAGS"]:
                for k, v in d.items():
                    if any(p in k.lower() for p in pats):
                        cand[k] += 1
                        if v and k not in example:
                            example[k] = v
            bp = g.bioproject.iloc[0]
            print(f"  {cid:<20} n={len(g):<6} {bp}")
            if cand:
                for k, n in cand.most_common(5):
                    print(f"      candidate tag: {k:<32} e.g. {str(example.get(k,''))[:38]}")
            else:
                print("      no candidate tag found - this study deposited nothing usable")
        print("\nAdd the candidate tags to the matching *_TAGS list in "
              "harvest_oral_age_metadata.py, then re-run finalize_dataset.sh")
        return

    ap.print_help()


if __name__ == "__main__":
    main()
