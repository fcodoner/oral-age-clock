#!/usr/bin/env python3
# --------------------------------------------------------------------------
# curate_healthy_adults.py
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
curate_healthy_adults.py
----------------------------------------------------------------------
Phase 1 - Salivary microbiome age model (IATA-CSIC)
Turns manifest_runs.tsv (from harvest_oral_age_metadata.py v1.1) into the
curated healthy-adult analysis set, plus the accession lists to download and
the counts needed for the paper's flow diagram.

  python curate_healthy_adults.py manifest_runs.tsv --out curated \
         [--tier core|lenient|strict] [--rules rules.tsv] [--min-reads 5000]

Tiers (identical code path, different thresholds):
  lenient : antibiotics < 1 month excluded; only explicit disease excluded
  core    : antibiotics < 3 months excluded; disease and pregnancy excluded   <- primary analysis
  strict  : antibiotics < 6 months excluded; unknown health status excluded, but only
            inside cohorts that actually record health (a cohort with no health metadata
            is kept whole and flagged, instead of being deleted)

Outputs (in --out):
  curated_samples.tsv        the analysis set, one row per participant
  excluded_samples.tsv       everything dropped, with the reason
  flow_counts.tsv            step-by-step attrition, per cohort (paper Figure 1)
  cohort_summary.tsv         per-cohort N, age, sex, sample type, platform
  age_distribution.tsv       N per cohort x age decade
  accessions_16S.txt         runs to download (amplicon)
  accessions_shotgun.txt     runs to download (WGS)
  curation_report.txt        human-readable summary of every decision

Per-cohort overrides: --rules rules.tsv with columns
  cohort_id  keep_ids_file  drop_ids_file  notes
keep_ids_file = one sample_alias/subject per line to KEEP (e.g. TRIAD controls).
Standard library + pandas.
"""
import argparse, os, sys
import pandas as pd

TIERS = {"lenient": dict(abx_months=1, drop_unknown_health=False, drop_pregnant=True),
         "core":    dict(abx_months=3, drop_unknown_health=False, drop_pregnant=True),
         "strict":  dict(abx_months=6, drop_unknown_health=True,  drop_pregnant=True)}

AGE_MIN, AGE_MAX = 18, 90
ORAL_WORDS = ("saliva", "oral", "mouth", "rinse", "tongue", "buccal", "spit")
# antibiotic_history vocabularies -> months since last course
ABX_MONTHS = {"week": 0.25, "month": 1, "6 months": 6, "year": 12,
              "i have not taken antibiotics in the past year.": 12}
DIAGNOSED = "diagnosed by a medical professional"
DISEASE_COLS = ["x_ibd", "x_diabetes", "x_cardiovascular_disease"]
EXTRA_DISEASE_ATTRS = ["cancer", "autoimmune", "liver_disease", "kidney_disease", "lung_disease"]


def attr(series, key):
    """Pull one tag out of the all_attributes column."""
    def g(s):
        if not isinstance(s, str):
            return None
        for kv in s.split(" | "):
            k, _, v = kv.partition("=")
            if k.strip().lower() == key:
                return v.strip()
        return None
    return series.map(g)


def abx_months_since(v):
    if not isinstance(v, str) or not v.strip():
        return None
    s = v.strip().lower()
    if s in ABX_MONTHS:
        return ABX_MONTHS[s]
    if "not taken" in s or "never" in s or "no" == s:
        return 12
    for k, m in ABX_MONTHS.items():
        if k in s:
            return m
    return None


def is_diagnosed(v):
    return isinstance(v, str) and DIAGNOSED in v.lower()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--out", default="curated")
    ap.add_argument("--tier", default="core", choices=list(TIERS))
    ap.add_argument("--rules")
    ap.add_argument("--min-reads", type=int, default=5000)
    ap.add_argument("--age-min", type=float, default=AGE_MIN)
    ap.add_argument("--age-max", type=float, default=AGE_MAX)
    a = ap.parse_args()
    cfg = TIERS[a.tier]
    os.makedirs(a.out, exist_ok=True)
    log = []

    m = pd.read_csv(a.manifest, sep="\t", low_memory=False)
    m["age"] = pd.to_numeric(m.age_years, errors="coerce")
    m["reads"] = pd.to_numeric(m.read_count, errors="coerce")
    for k in EXTRA_DISEASE_ATTRS + ["pregnant"]:
        m["y_" + k] = attr(m.all_attributes, k)
    m["abx_months"] = m.antibiotics.map(abx_months_since)
    # sample type: repositories mislabel; body_site wins over env_medium
    site = m.get("x_body_site").fillna("") if "x_body_site" in m else ""
    mat = (m.sample_material.fillna("") + " " + site).str.lower()
    m["sample_type"] = "unknown"
    m.loc[mat.str.contains("saliva|spit|rinse|mouthwash|envo:02000036|oral fluid|oral cavity|^oral$|^oral |^oral,"), "sample_type"] = "saliva/rinse"
    m.loc[mat.str.contains(r"tongue|swab|buccal|cheek|mouth(?!wash)"), "sample_type"] = "swab"
    m.loc[mat.str.contains("oropharyn|nasopharyn|plaque|stool|skin|water|neoplasm|carcinoma|cancer|tumou?r"), "sample_type"] = "unknown"

    keep_rules, drop_rules = {}, {}
    if a.rules:
        for r in pd.read_csv(a.rules, sep="\t").to_dict("records"):
            c = r["cohort_id"]
            for col, store in (("keep_ids_file", keep_rules), ("drop_ids_file", drop_rules)):
                f = str(r.get(col) or "").strip()
                if f and os.path.exists(f):
                    store[c] = {x.strip() for x in open(f) if x.strip()}

    steps, excluded = [], []

    def drop(mask, reason):
        nonlocal m
        if mask.any():
            d = m[mask].copy()
            d["exclusion_reason"] = reason
            excluded.append(d)
        m = m[~mask].copy()
        steps.append((reason, int(mask.sum()), len(m)))
        log.append(f"{reason:<42} removed {int(mask.sum()):>5}   remaining {len(m):>5}")

    log.append(f"Tier: {a.tier}  (antibiotic window {cfg['abx_months']} months, "
               f"age {a.age_min:g}-{a.age_max:g}, min reads {a.min_reads})")
    log.append(f"Input runs: {len(m)}\n")
    steps.append(("0. Input runs", 0, len(m)))

    drop(~m.filter_status.eq("PASS") & m.age.isna(), "1. No usable numeric age")
    drop(m.age.notna() & ~m.age.between(a.age_min, a.age_max), "2. Outside adult age window")
    drop(~m.sample_type.isin(["saliva/rinse", "swab"]), "3. Not an oral sample / type unknown")
    drop(~m.library_strategy.isin(["AMPLICON", "WGS", "OTHER"]), "4. Unsuitable library strategy")
    drop(m.reads.fillna(0) < a.min_reads, "5. Too few reads")

    for c, ids in drop_rules.items():
        drop(m.cohort_id.eq(c) & (m.sample_alias.isin(ids) | m.subject_id.isin(ids)),
             f"6. Per-cohort drop list ({c})")
    for c, ids in keep_rules.items():
        drop(m.cohort_id.eq(c) & ~(m.sample_alias.isin(ids) | m.subject_id.isin(ids)),
             f"6. Not on keep list, e.g. non-controls ({c})")

    drop(m.abx_months.notna() & (m.abx_months < cfg["abx_months"]),
         f"7. Antibiotics within {cfg['abx_months']} months")
    if cfg["drop_unknown_health"]:
        known = m[DISEASE_COLS].notna().any(axis=1) | m.health_status_raw.notna()
        # only enforce in cohorts that record health for at least half their samples
        cov = known.groupby(m.cohort_id).mean()
        enforce = m.cohort_id.map(cov).fillna(0) >= 0.5
        for c in cov[cov < 0.5].index:
            log.append(f"   note: cohort {c} records no health metadata - kept whole, "
                       f"flag as a limitation rather than dropping it")
        drop(~known & enforce, "8. Health status not recorded (strict)")

    dis = pd.Series(False, index=m.index)
    for col in DISEASE_COLS + ["y_" + k for k in EXTRA_DISEASE_ATTRS]:
        if col in m:
            dis |= m[col].map(is_diagnosed)
    if "x_subset_healthy" in m:
        dis |= m.x_subset_healthy.astype(str).str.lower().eq("false") & (a.tier == "strict")
    drop(dis, "9. Diagnosed disease reported")
    if cfg["drop_pregnant"]:
        drop(m.y_pregnant.astype(str).str.lower().isin(["true", "yes"]), "10. Pregnancy reported")

    # one run per sample, one sample per participant, within each sequencing type
    m = m.sort_values("reads", ascending=False)
    before = len(m)
    m = m.drop_duplicates(subset=["library_strategy" if False else "sample_accession"], keep="first")
    m["is_wgs"] = m.library_strategy.eq("WGS")
    m["participant"] = m.subject_id.fillna("").replace("", pd.NA).fillna(m.sample_accession)
    m = m.drop_duplicates(subset=["is_wgs", "cohort_id", "participant"], keep="first")
    steps.append(("11. Duplicate runs / repeat participants", before - len(m), len(m)))
    log.append(f"{'11. Duplicate runs / repeat participants':<42} removed {before-len(m):>5}   remaining {len(m):>5}")

    # ---- outputs ----
    m = m.sort_values(["cohort_id", "age"])
    m.to_csv(os.path.join(a.out, "curated_samples.tsv"), sep="\t", index=False)
    if excluded:
        pd.concat(excluded).to_csv(os.path.join(a.out, "excluded_samples.tsv"), sep="\t", index=False)
    pd.DataFrame(steps, columns=["step", "removed", "remaining"]).to_csv(
        os.path.join(a.out, "flow_counts.tsv"), sep="\t", index=False)

    summ = m.groupby("cohort_id").agg(
        n=("run_accession", "size"), age_min=("age", "min"), age_median=("age", "median"),
        age_max=("age", "max"), pct_female=("sex", lambda s: round(100 * s.eq("female").mean(), 1)),
        sample_types=("sample_type", lambda s: "; ".join(f"{k}:{v}" for k, v in s.value_counts().items())),
        platforms=("instrument_model", lambda s: "; ".join(sorted(set(s.dropna())))),
        strategies=("library_strategy", lambda s: "; ".join(f"{k}:{v}" for k, v in s.value_counts().items())),
        countries=("country", lambda s: "; ".join(sorted({str(x).split(":")[0] for x in s.dropna()})[:4])))
    summ.to_csv(os.path.join(a.out, "cohort_summary.tsv"), sep="\t")

    bins = [18, 30, 40, 50, 60, 70, 200]
    labs = ["18-29", "30-39", "40-49", "50-59", "60-69", "70+"]
    dist = pd.crosstab(m.cohort_id, pd.cut(m.age, bins, right=False, labels=labs), margins=True)
    dist.to_csv(os.path.join(a.out, "age_distribution.tsv"), sep="\t")

    amp = m[~m.is_wgs].run_accession.tolist()
    wgs = m[m.is_wgs].run_accession.tolist()
    open(os.path.join(a.out, "accessions_16S.txt"), "w").write("\n".join(amp) + "\n")
    open(os.path.join(a.out, "accessions_shotgun.txt"), "w").write("\n".join(wgs) + "\n")

    log.append("")
    log.append(f"FINAL: {len(m)} participants  ({len(amp)} 16S, {len(wgs)} shotgun) "
               f"across {m.cohort_id.nunique()} cohorts")
    log.append("")
    log.append(summ.to_string())
    log.append("")
    log.append(dist.to_string())
    n_small = (summ.n < 30).sum()
    if n_small:
        log.append(f"\nWARNING: {n_small} cohort(s) with <30 participants - too small to act as a "
                   f"leave-one-cohort-out fold on their own; pool them by region/protocol.")
    if m.cohort_id.nunique() < 6:
        log.append(f"WARNING: only {m.cohort_id.nunique()} cohorts. Leave-one-cohort-out needs "
                   f"6-8 to be convincing. Add cohorts before modelling.")
    txt = "\n".join(log)
    open(os.path.join(a.out, "curation_report.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
