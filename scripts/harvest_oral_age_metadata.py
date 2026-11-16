#!/usr/bin/env python3
# --------------------------------------------------------------------------
# harvest_oral_age_metadata.py
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
harvest_oral_age_metadata.py
----------------------------------------------------------------------
Phase 1 - Salivary microbiome age model (IATA-CSIC)
Builds a per-sample / per-run metadata manifest for a list of INSDC
BioProjects (PRJNA / PRJEB / PRJDB) and writes the run accessions that pass
the basic inclusion filters, ready for a downloader.

Sources (public, no API key needed):
  1. ENA Portal API  -> run list, library info, FASTQ links, MD5
  2. ENA Browser API -> full sample attribute sets (works for SAMN/SAMEA/SAMD)

Usage:
  python harvest_oral_age_metadata.py studies.tsv --out out_dir
  studies.tsv: tab-separated, header required, columns:
      bioproject  cohort_id  include(Y/N)  tax_id(optional)  [notes]
  tax_id example: 447426 = human oral metagenome (use for mixed-site projects like AGP)

Outputs (in --out):
  manifest_runs.tsv       one row per run, raw + harmonised fields
  attribute_census.tsv    which attribute tags each study uses (for manual mapping)
  accessions_included.txt run accessions passing filters (one per line)
  accessions_by_study/    same, split per BioProject

Standard library only (Python >= 3.8).
"""
import argparse, csv, os, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_XML = "https://www.ebi.ac.uk/ena/browser/api/xml/"
RUN_FIELDS = ["run_accession", "sample_accession", "secondary_sample_accession",
              "study_accession", "sample_alias", "sample_title", "library_strategy",
              "library_layout", "library_source", "instrument_platform", "instrument_model",
              "read_count", "base_count", "target_gene", "collection_date",
              "fastq_ftp", "fastq_md5", "fastq_bytes"]
BATCH = 100
PAUSE = 0.4           # be polite to the APIs
ADULT_MIN, ADULT_MAX = 18, 90

# ---------------- network ----------------
def http_get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "iata-csic-oral-age/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))

def ena_runs(bioproject, tax_id=""):
    """All runs of a project; optional tax_id restricts to e.g. 447426 (human oral metagenome)."""
    if tax_id:
        q = urllib.parse.urlencode({"result": "read_run", "format": "tsv", "limit": 0,
                                    "fields": ",".join(RUN_FIELDS),
                                    "query": f'study_accession="{bioproject}" AND tax_eq({tax_id})'})
        txt = http_get(f"{ENA_SEARCH}?{q}")
    else:
        q = urllib.parse.urlencode({"accession": bioproject, "result": "read_run",
                                    "fields": ",".join(RUN_FIELDS), "format": "tsv",
                                    "limit": 0})
        txt = http_get(f"{ENA_FILEREPORT}?{q}")
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    return rows

def ena_sample_attributes(accessions):
    """Return {sample_accession: {tag_lower: value}} for a list of sample accessions."""
    out = {}
    for i in range(0, len(accessions), BATCH):
        chunk = accessions[i:i + BATCH]
        xml = http_get(ENA_XML + ",".join(chunk))
        out.update(parse_sample_xml(xml))
        time.sleep(PAUSE)
    return out

def parse_sample_xml(xml_text):
    res = {}
    root = ET.fromstring(xml_text)
    for s in root.iter("SAMPLE"):
        acc = s.get("accession", "")
        ids = s.find("IDENTIFIERS")
        ext = []
        if ids is not None:
            ext = [e.text for e in ids.iter("EXTERNAL_ID") if e.text]
        attrs = {}
        title = s.findtext("TITLE") or ""
        attrs["_title"] = title
        for a in s.iter("SAMPLE_ATTRIBUTE"):
            tag = (a.findtext("TAG") or "").strip().lower()
            val = (a.findtext("VALUE") or "").strip()
            unit = (a.findtext("UNITS") or "").strip()
            if tag:
                attrs[tag] = f"{val} {unit}".strip() if unit else val
        for key in [acc] + ext:
            if key:
                res[key] = attrs
    return res

# ---------------- harmonisation ----------------
AGE_TAGS = ["host age", "host_age", "age", "age_years", "age (years)", "subject age",
            "host_age_years", "patient age"]
AGE_BIN_TAGS = ["age_cat", "age category", "age_group", "age group", "age range"]
BIRTH_TAGS = ["birth_year", "year of birth", "host birth year"]
SEX_TAGS = ["host sex", "host_sex", "sex", "gender", "host_gender", "sex_s"]
COUNTRY_TAGS = ["geographic location (country and/or sea)", "geographic location",
                "geo_loc_name", "country_residence", "country"]
SITE_TAGS = ["environmental medium", "env_material", "host body product", "host tissue sampled", "body_site", "isolation source", "isolation_source", "sample_type", "host_body_site", "local-scale environmental context", "broad-scale environmental context", "env_medium", "env_local_scale", "env_broad_scale", "environment (material)", "tissue", "tissue_type", "sample_site", "env_biome", "host body site"]
BMI_TAGS = ["host body mass index", "host_body_mass_index", "bmi", "host_bmi"]
SMOKE_TAGS = ["smoking_frequency", "smoker", "smoking", "tobacco use", "host smoking status"]
ABX_TAGS = ["antibiotic_history", "antibiotics", "antibiotic use", "host medication"]
HEALTH_TAGS = ["host disease", "host_disease", "host_disease_stat", "disease", "disease status",
               "host health state", "host_health_state", "health_status", "host_status",
               "diagnosis", "group", "case_control"]
SUBJECT_TAGS = ["host_subject_id", "host subject id", "anonymized_name", "subject_id"]
AGP_EXTRA = ["subset_healthy", "ibd", "diabetes", "cardiovascular_disease", "body_site"]
DATE_TAGS = ["collection date", "collection_date", "collection_timestamp"]

def first(attrs, tags):
    for t in tags:
        v = attrs.get(t)
        if v and v.lower() not in {"missing", "not collected", "not applicable", "na",
                                   "n/a", "unknown", "unspecified", "not provided", "none", ""}:
            return v, t
    return "", ""

def parse_age(value):
    """'15.499 years' -> 15.499 ; '51' -> 51 ; '18 months' -> 1.5 ; '31-40' -> None (binned)."""
    v = value.lower().strip()
    if re.fullmatch(r"\d+\s*-\s*\d+.*", v) or "+" in v or "<" in v or ">" in v:
        return None
    m = re.search(r"(\d+(?:\.\d+)?)", v)
    if not m:
        return None
    x = float(m.group(1))
    if "month" in v:
        x /= 12.0
    elif "week" in v:
        x /= 52.18
    elif "day" in v:
        x /= 365.25
    return round(x, 2)

def year_of(date_str):
    m = re.search(r"(19|20)\d{2}", date_str or "")
    return int(m.group(0)) if m else None

def harmonise(run, attrs):
    h = {}
    age_raw, age_tag = first(attrs, AGE_TAGS)
    age = parse_age(age_raw) if age_raw else None
    precision, source = "", ""
    if age is not None:
        precision, source = "exact", f"attribute:{age_tag}"
    else:
        by, by_tag = first(attrs, BIRTH_TAGS)
        cd, _ = first(attrs, DATE_TAGS)
        cy = year_of(cd) or year_of(run.get("collection_date", ""))
        if by and year_of(by) and cy:
            age = cy - year_of(by)
            precision, source = "derived (+/-1 y)", f"{by_tag} & collection year"
    age_bin, bin_tag = first(attrs, AGE_BIN_TAGS)
    if age is None and not age_bin and age_raw:
        age_bin = age_raw            # e.g. '31-40' stored under 'host age'
    h["age_years"] = "" if age is None else age
    h["age_precision"] = precision or ("binned" if age_bin else "missing")
    h["age_source"] = source or (f"attribute:{bin_tag or age_tag}" if age_bin else "")
    h["age_bin"] = age_bin
    h["sex"] = norm_sex(first(attrs, SEX_TAGS)[0])
    h["country"] = first(attrs, COUNTRY_TAGS)[0]
    h["sample_material"] = first(attrs, SITE_TAGS)[0]
    h["bmi"] = first(attrs, BMI_TAGS)[0]
    h["smoking"] = first(attrs, SMOKE_TAGS)[0]
    h["antibiotics"] = first(attrs, ABX_TAGS)[0]
    h["health_status_raw"] = first(attrs, HEALTH_TAGS)[0]
    h["collection_date"] = first(attrs, DATE_TAGS)[0] or run.get("collection_date", "")
    h["sample_title"] = attrs.get("_title", "")
    h["subject_id"] = first(attrs, SUBJECT_TAGS)[0]
    for k in AGP_EXTRA:
        h["x_" + k] = attrs.get(k, "")
    return h

def norm_sex(v):
    v = v.lower()
    if v.startswith(("f", "wom")):
        return "female"
    if v.startswith(("m", "man")) and not v.startswith("mi"):
        return "male"
    return v

def basic_filter(run, h):
    reasons = []
    if h["age_years"] == "":
        reasons.append("no_numeric_age")
    elif not (ADULT_MIN <= float(h["age_years"]) <= ADULT_MAX):
        reasons.append("outside_18_90")
    # v1.1: only judge sample type when a site/material attribute exists
    mat = h["sample_material"].lower()
    if mat and not any(k in mat for k in ["saliva", "oral", "mouth", "rinse", "tongue", "buccal"]):
        reasons.append("check_sample_type")
    if run.get("library_strategy", "") not in {"AMPLICON", "WGS", "OTHER"}:
        reasons.append("library_strategy")
    try:
        if int(run.get("read_count") or 0) < 5000:
            reasons.append("low_reads")
    except ValueError:
        pass
    return "PASS" if not reasons else ";".join(reasons)

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="?")
    ap.add_argument("--out", default="harvest_out")
    ap.add_argument("--rebuild", help="existing manifest_runs.tsv to re-harmonise offline")
    a = ap.parse_args()
    if a.rebuild:
        rebuild(a.rebuild, a.out)
        return
    os.makedirs(os.path.join(a.out, "accessions_by_study"), exist_ok=True)

    studies = [r for r in csv.DictReader(open(a.studies), delimiter="\t")
               if r.get("include", "Y").upper().startswith("Y")]
    all_rows, census = [], defaultdict(Counter)
    for st in studies:
        bp = st["bioproject"].strip()
        print(f"[{bp}] fetching runs ...", file=sys.stderr)
        runs = ena_runs(bp, st.get("tax_id", "").strip())
        samples = sorted({r["sample_accession"] for r in runs if r.get("sample_accession")})
        print(f"[{bp}] {len(runs)} runs / {len(samples)} samples; fetching attributes ...",
              file=sys.stderr)
        attrs = ena_sample_attributes(samples)
        for r in runs:
            at = attrs.get(r["sample_accession"]) or attrs.get(r.get("secondary_sample_accession", "")) or {}
            for tag in at:
                census[bp][tag] += 1
            h = harmonise(r, at)
            row = {"cohort_id": st.get("cohort_id", bp), "bioproject": bp, **r, **h}
            row["filter_status"] = basic_filter(r, h)
            row["all_attributes"] = " | ".join(f"{k}={v}" for k, v in at.items() if k != "_title")
            all_rows.append(row)
        passed = [r["run_accession"] for r in all_rows if r["bioproject"] == bp and r["filter_status"] == "PASS"]
        open(os.path.join(a.out, "accessions_by_study", f"{bp}.txt"), "w").write("\n".join(passed) + "\n")
        print(f"[{bp}] {len(passed)} runs pass basic filters", file=sys.stderr)

    mark_primary(all_rows)
    write_outputs(all_rows, a.out)

def mark_primary(rows):
    """One run per sample (most reads), one sample per subject (most reads); 16S and WGS kept apart."""
    best_run, best_subj = {}, {}
    for r in rows:
        if r["filter_status"] != "PASS":
            continue
        key = (r["library_strategy"] == "WGS", r["sample_accession"])
        if key not in best_run or rd(r) > rd(best_run[key]):
            best_run[key] = r
    for r in best_run.values():
        subj = (r["library_strategy"] == "WGS", r["cohort_id"], r.get("subject_id") or r["sample_accession"])
        if subj not in best_subj or rd(r) > rd(best_subj[subj]):
            best_subj[subj] = r
    keep = {id(r) for r in best_subj.values()}
    for r in rows:
        r["primary"] = "Y" if id(r) in keep else "N"

def rd(r):
    try:
        return int(r.get("read_count") or 0)
    except ValueError:
        return 0

def write_outputs(all_rows, out):
    os.makedirs(os.path.join(out, "accessions_by_study"), exist_ok=True)
    cols = ["cohort_id", "bioproject", "run_accession", "sample_accession", "sample_alias", "subject_id",
            "age_years", "age_precision", "age_bin", "age_source", "sex", "country",
            "sample_material", "bmi", "smoking", "antibiotics", "health_status_raw",
            "collection_date", "library_strategy", "library_layout", "instrument_model",
            "target_gene", "read_count", "base_count", "fastq_ftp", "fastq_md5",
            "filter_status", "primary", "x_subset_healthy", "x_ibd", "x_diabetes",
            "x_cardiovascular_disease", "x_body_site", "sample_title", "all_attributes"]
    with open(os.path.join(out, "manifest_runs.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    amp = [r["run_accession"] for r in all_rows if r["primary"] == "Y" and r["library_strategy"] != "WGS"]
    wgs = [r["run_accession"] for r in all_rows if r["primary"] == "Y" and r["library_strategy"] == "WGS"]
    open(os.path.join(out, "accessions_included.txt"), "w").write("\n".join(amp) + "\n")
    open(os.path.join(out, "accessions_included_wgs.txt"), "w").write("\n".join(wgs) + "\n")
    for bp in sorted({r["bioproject"] for r in all_rows}):
        ids = [r["run_accession"] for r in all_rows if r["bioproject"] == bp and r["primary"] == "Y"]
        open(os.path.join(out, "accessions_by_study", f"{bp}.txt"), "w").write("\n".join(ids) + "\n")
    print(f"Done: {len(all_rows)} runs; {len(amp)} primary 16S + {len(wgs)} primary shotgun "
          f"(one run per sample, one sample per participant) -> {out}", file=sys.stderr)

def rebuild(manifest, out):
    """Re-harmonise an existing manifest offline from its all_attributes column (no downloads)."""
    rows = list(csv.DictReader(open(manifest), delimiter="\t"))
    for r in rows:
        at = {}
        for kv in (r.get("all_attributes") or "").split(" | "):
            if "=" in kv:
                k, v = kv.split("=", 1)
                at[k.strip().lower()] = v.strip()
        at["_title"] = r.get("sample_title", "")
        r.update(harmonise(r, at))
        r["filter_status"] = basic_filter(r, r)
    mark_primary(rows)
    write_outputs(rows, out)

if __name__ == "__main__":
    main()
