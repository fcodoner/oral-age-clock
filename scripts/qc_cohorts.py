#!/usr/bin/env python3
# --------------------------------------------------------------------------
# qc_cohorts.py
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
qc_cohorts.py
----------------------------------------------------------------------
Phase 1 - Salivary microbiome age model (IATA-CSIC)
Quality-control pass over the curated set, run AFTER curate_healthy_adults.py
and BEFORE downloading or modelling. It looks for the five problems that would
otherwise corrupt a leave-one-cohort-out design, and writes the fixes.

  python qc_cohorts.py metadata/curated_core \
         --excluded metadata/curated_core/excluded_samples.tsv \
         --manifest metadata/harvest_<date>_v2/manifest_runs.tsv \
         --out metadata/qc

Checks
  1 duplicate cohorts   same participants deposited under two BioProjects
  2 invalid ages        a single value repeated for a whole cohort, or an
                        implausibly narrow spread (placeholder, not real age)
  3 lost samples        why step 3 of curation dropped so many as "not oral"
  4 missing covariates  cohorts where sex (or another field) is absent for all
  5 small cohorts       too small to be their own validation fold

Outputs
  qc_report.txt          everything below, in readable form
  cohorts_to_drop.txt    one cohort_id per line, with the reason as a comment
  pooling_map.tsv        suggested grouping of small cohorts into folds
  tag_candidates.tsv     metadata labels seen in dropped samples that may mean
                         "oral sample" or "sex" and are not being read yet
  curated_samples_qc.tsv the curated set with drops applied and fold_id added

Nothing is deleted automatically: review cohorts_to_drop.txt, then re-run with
--apply to write curated_samples_qc.tsv.
"""
import argparse, os, sys
import pandas as pd

ORAL_WORDS = ("saliva", "oral", "mouth", "rinse", "tongue", "buccal", "spit", "mouthwash")
SEX_WORDS = ("sex", "gender")


def sec(title):
    return f"\n{'='*78}\n{title}\n{'='*78}"


def age_match(c, a, b):
    """Of the sample names two cohorts share, what fraction carry the same age?
    Generic names (S1, sample3) collide by chance and will not match on age."""
    ga = c[c.cohort_id == a].set_index(c[c.cohort_id == a].sample_alias.astype(str))
    gb = c[c.cohort_id == b].set_index(c[c.cohort_id == b].sample_alias.astype(str))
    shared = [x for x in set(ga.index) & set(gb.index)]
    if not shared:
        return 0.0
    same = 0
    for k in shared:
        try:
            if abs(float(pd.Series(ga.loc[k, "age_years"]).iloc[0])
                   - float(pd.Series(gb.loc[k, "age_years"]).iloc[0])) < 0.51:
                same += 1
        except Exception:
            pass
    return same / len(shared)


def check_duplicates(c, log):
    """Two cohorts holding the same participants: by accession overlap, then by profile."""
    log.append(sec("CHECK 1 - duplicate cohorts"))
    drops, cohorts = [], sorted(c.cohort_id.unique())
    accs = {k: set(g.sample_accession.dropna()) for k, g in c.groupby("cohort_id")}
    alias = {k: set(g.sample_alias.dropna().astype(str)) for k, g in c.groupby("cohort_id")}
    found = False
    for i, a in enumerate(cohorts):
        for b in cohorts[i + 1:]:
            ov = len(accs[a] & accs[b])
            ova = len(alias[a] & alias[b])
            if len(accs[a]) != len(accs[b]):
                smaller = a if len(accs[a]) < len(accs[b]) else b
            else:                       # same size: drop the alphabetically later name
                smaller = max(a, b)
            n_small = max(1, min(len(accs[a]), len(accs[b])))
            frac = ova / n_small
            if ov:
                found = True
                log.append(f"  {a} and {b}: {ov} shared sample ACCESSIONS -> "
                           f"certain duplicate, drop {smaller}")
                drops.append((smaller, f"shares deposited samples with {a if smaller==b else b}"))
            elif ova and frac >= 0.5 and ova >= 10 and age_match(c, a, b) >= 0.8:
                found = True
                log.append(f"  {a} and {b}: {ova} shared sample names "
                           f"({frac:.0%} of the smaller cohort) AND the same ages attached "
                           f"to them -> likely duplicate, drop {smaller}")
                drops.append((smaller, f"sample names and ages overlap with "
                                       f"{a if smaller==b else b}"))
            elif ova:
                log.append(f"  note: {a} and {b} share {ova} sample name(s) "
                           f"({frac:.0%} of the smaller cohort) but no deposited samples. "
                           f"Generic names like 'S1' or 'sample3' collide by chance - "
                           f"NOT treated as duplicates.")
    # profile fingerprint: identical n, median age, sex ratio and decade counts
    prof = c.groupby("cohort_id").apply(
        lambda g: (len(g), round(g.age_years.astype(float).median(), 3),
                   round(g.sex.eq("female").mean(), 4),
                   tuple(pd.cut(g.age_years.astype(float), [18, 30, 40, 50, 60, 70, 200],
                                right=False).value_counts().sort_index())), include_groups=False)
    seen, already = {}, {cid for cid, _ in drops}
    for cid, p in prof.items():
        if cid in already:
            continue
        if p in seen and seen[p] not in already:
            found = True
            log.append(f"  {seen[p]} and {cid}: identical size, median age, sex ratio and age "
                       f"distribution -> near-certain duplicate, drop {cid}")
            drops.append((cid, f"identical profile to {seen[p]}"))
        else:
            seen[p] = cid
    if not found:
        log.append("  none found")
    return drops


def check_ages(c, log):
    log.append(sec("CHECK 2 - implausible age distributions"))
    drops = []
    for cid, g in c.groupby("cohort_id"):
        age = g.age_years.astype(float)
        if len(g) < 5:
            continue
        if age.nunique() == 1:
            log.append(f"  {cid}: every one of {len(g)} participants has age {age.iloc[0]:g} "
                       f"-> placeholder, not real age. DROP")
            drops.append((cid, "single age value for the whole cohort"))
        elif age.nunique() <= 3 and len(g) >= 20:
            log.append(f"  {cid}: only {age.nunique()} distinct ages among {len(g)} participants "
                       f"({sorted(age.unique())}) -> check the source before using")
        elif age.std() < 2 and len(g) >= 20:
            log.append(f"  {cid}: age range {age.min():g}-{age.max():g}, sd {age.std():.2f} "
                       f"-> suspiciously narrow, check the source")
    if not drops:
        log.append("  no cohort had a single repeated age value")
    return drops


def check_lost(excluded, log, out):
    log.append(sec("CHECK 3 - samples dropped as 'not an oral sample / type unknown'"))
    if excluded is None:
        log.append("  excluded_samples.tsv not supplied - skipped")
        return
    x = excluded[excluded.exclusion_reason.astype(str).str.contains("oral", case=False)]
    log.append(f"  {len(x)} samples dropped by this rule")
    if x.empty:
        return
    log.append("\n  by cohort (top 20):")
    for cid, n in x.groupby("cohort_id").size().sort_values(ascending=False).head(20).items():
        log.append(f"    {cid:<20} {n:>6}")
    log.append("\n  the material values that caused it (top 25):")
    vc = x.sample_material.fillna("(blank)").value_counts().head(25)
    for v, n in vc.items():
        flag = "  <-- looks oral, tag may be misread" if any(w in str(v).lower() for w in ORAL_WORDS) else ""
        log.append(f"    {str(v)[:60]:<62} {n:>6}{flag}")
    blank = int(x.sample_material.fillna("").astype(str).str.strip().eq("").sum())
    if blank:
        log.append(f"\n  {blank} of these had NO material recorded at all. Those cohorts probably "
                   f"store sample type under a tag the harvester does not read yet.")
    cand = (x[x.sample_material.fillna("").astype(str).str.strip().eq("")]
            .groupby("cohort_id").size().sort_values(ascending=False))
    if len(cand):
        cand.rename("n_samples_with_no_material").to_csv(
            os.path.join(out, "cohorts_missing_sample_type.tsv"), sep="\t")
        log.append(f"  -> written to cohorts_missing_sample_type.tsv "
                   f"({len(cand)} cohorts). Look these up in attribute_census.tsv.")


def check_covariates(c, log):
    log.append(sec("CHECK 4 - cohorts missing a covariate entirely"))
    for field in ["sex", "country"]:
        miss = []
        for cid, g in c.groupby("cohort_id"):
            v = g[field].fillna("").astype(str).str.strip()
            if v.eq("").all():
                miss.append((cid, len(g), "field empty for all participants"))
            elif field == "sex" and v.eq("male").mean() == 1.0 and len(g) >= 20:
                miss.append((cid, len(g), "100% male - verify, may be a misread tag"))
            elif field == "sex" and v.eq("female").mean() == 1.0 and len(g) >= 20:
                miss.append((cid, len(g), "100% female - verify, may be genuine (e.g. a women's cohort)"))
        if miss:
            log.append(f"\n  {field}:")
            for cid, n, why in miss:
                log.append(f"    {cid:<20} n={n:<5} {why}")
        else:
            log.append(f"\n  {field}: recorded in every cohort")
    log.append("\n  For any cohort listed: search attribute_census.tsv for that BioProject and see "
               "which tag holds the value, then add it to the tag list in the harvester and "
               "re-run with --rebuild.")


def check_small(c, log, out, min_n=30):
    log.append(sec(f"CHECK 5 - cohorts with fewer than {min_n} participants"))
    sizes = c.groupby("cohort_id").size().sort_values()
    small = sizes[sizes < min_n]
    if small.empty:
        log.append("  none")
        return
    def mode_or_blank(s):
        m = s.dropna().astype(str)
        m = m[m.str.strip() != ""]
        return m.mode().iloc[0] if len(m) else ""

    meta = c.groupby("cohort_id").agg(country=("country", mode_or_blank),
                                      platform=("instrument_model", mode_or_blank))
    rows = []
    for cid, n in small.items():
        key = f"{meta.loc[cid,'country'] or 'unknown'}|{meta.loc[cid,'platform'] or 'unknown'}"
        rows.append({"cohort_id": cid, "n": n, "country": meta.loc[cid, "country"],
                     "platform": meta.loc[cid, "platform"], "suggested_pool": f"POOL_{key}"})
        log.append(f"  {cid:<20} n={n:<5} {meta.loc[cid,'country'] or '?'} / {meta.loc[cid,'platform'] or '?'}")
    pool = pd.DataFrame(rows)
    live = os.path.join(out, "pooling_map.tsv")
    if os.path.exists(live):
        pool.to_csv(os.path.join(out, "pooling_map_suggested.tsv"), sep="\t", index=False)
        log.append("\n  pooling_map.tsv already exists and was NOT overwritten (your edits are "
                   "kept). Fresh suggestions written to pooling_map_suggested.tsv.")
    else:
        pool.to_csv(live, sep="\t", index=False)
    log.append("\n  A cohort this small cannot be its own validation fold. Either pool it with "
               "cohorts of the same country and platform (see pooling_map.tsv) or exclude it. "
               "Pooling keeps the data; excluding keeps the folds clean. Pool when the pooled "
               "group reaches 30+, exclude otherwise.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("curated_dir")
    ap.add_argument("--excluded")
    ap.add_argument("--manifest")
    ap.add_argument("--out", default="metadata/qc")
    ap.add_argument("--min-n", type=int, default=30)
    ap.add_argument("--apply", action="store_true",
                    help="write curated_samples_qc.tsv with cohorts_to_drop.txt applied")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    cur = os.path.join(a.curated_dir, "curated_samples.tsv")
    c = pd.read_csv(cur, sep="\t", low_memory=False)
    for col in ["sex", "country", "sample_alias", "instrument_model", "sample_accession"]:
        if col not in c:
            c[col] = ""
    excluded = None
    if a.excluded and os.path.exists(a.excluded):
        excluded = pd.read_csv(a.excluded, sep="\t", low_memory=False)

    log = [f"QC pass on {cur}",
           f"{len(c)} participants across {c.cohort_id.nunique()} cohorts"]
    drops = check_duplicates(c, log) + check_ages(c, log)
    check_lost(excluded, log, a.out)
    check_covariates(c, log)
    check_small(c, log, a.out, a.min_n)

    seen, uniq = set(), []
    for cid, why in drops:
        if cid not in seen:
            seen.add(cid); uniq.append((cid, why))
    with open(os.path.join(a.out, "cohorts_to_drop.txt"), "w") as f:
        for cid, why in uniq:
            f.write(f"{cid}\t# {why}\n")
    log.append(sec("SUMMARY"))
    log.append(f"  {len(uniq)} cohort(s) proposed for removal: "
               + (", ".join(cid for cid, _ in uniq) or "none"))
    kept = c[~c.cohort_id.isin(seen)]
    log.append(f"  After removal: {len(kept)} participants across {kept.cohort_id.nunique()} cohorts")
    folds = kept.groupby("cohort_id").size()
    log.append(f"  Cohorts large enough to be their own fold (>={a.min_n}): {(folds>=a.min_n).sum()}")

    if a.apply:
        kept = kept.copy()
        pool_path = os.path.join(a.out, "pooling_map.tsv")
        kept["fold_id"] = kept.cohort_id
        if os.path.exists(pool_path):
            pm = pd.read_csv(pool_path, sep="\t").set_index("cohort_id")["suggested_pool"].to_dict()
            kept["fold_id"] = kept.cohort_id.map(lambda x: pm.get(x, x))
        kept.to_csv(os.path.join(a.out, "curated_samples_qc.tsv"), sep="\t", index=False)
        log.append(f"  Written: curated_samples_qc.tsv with a fold_id column "
                   f"({kept.fold_id.nunique()} folds)")
        for name, sub in [("accessions_16S_qc.txt", kept[kept.library_strategy != "WGS"]),
                          ("accessions_shotgun_qc.txt", kept[kept.library_strategy == "WGS"])]:
            open(os.path.join(a.out, name), "w").write("\n".join(sub.run_accession) + "\n")
        log.append("  Written: accessions_16S_qc.txt and accessions_shotgun_qc.txt")
    else:
        log.append("  Re-run with --apply once you have reviewed cohorts_to_drop.txt")

    txt = "\n".join(log)
    open(os.path.join(a.out, "qc_report.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
