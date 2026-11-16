#!/usr/bin/env python3
# --------------------------------------------------------------------------
# expand_studies.py
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
expand_studies.py
----------------------------------------------------------------------
Phase 1 - Salivary microbiome age model (IATA-CSIC)
Grows studies.tsv two ways:

  A) --from-supplement : mine accessions out of the source-paper supplementary
     tables (Ruan 2022 = 47 studies, Saliva MicroAge 2025 = 56 datasets,
     Huang 2020 = Qiita studies). Reads .xlsx, .csv, .tsv, .txt, .docx or .pdf,
     pulls every INSDC / Qiita / SRA identifier, and tags each with its source.

  B) --discover : ask ENA directly for every study containing human oral
     metagenome runs (taxon 447426), keep the ones whose title looks like adult
     saliva work, then screen a sample of each study's BioSamples to see whether
     age is actually deposited. This is reproducible and does not depend on any
     supplement being available.

Both modes write/extend studies.tsv (never overwriting your include= decisions)
and a candidate_studies.tsv review sheet with the evidence behind each row.

  python expand_studies.py --from-supplement ruan2022_tableS1.xlsx \
         --source RUAN2022 --studies metadata/studies.tsv
  python expand_studies.py --discover --studies metadata/studies.tsv \
         --screen 8 --min-runs 30
  python expand_studies.py --discover --no-screen     # fast, title filter only

Then run harvest_oral_age_metadata.py on the updated studies.tsv.
Requires: pandas; openpyxl for .xlsx; pdfplumber or pdftotext for .pdf (optional).
"""
import argparse, csv, os, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_XML = "https://www.ebi.ac.uk/ena/browser/api/xml/"
ORAL_TAXON = 447426                      # human oral metagenome
COLS = ["bioproject", "cohort_id", "include", "tax_id", "source", "notes"]

ACC_RE = re.compile(r"\b(PRJ(?:NA|EB|DB)\d{4,9}|(?:SRP|ERP|DRP)\d{5,9}|"
                    r"(?:SRA|ERA|DRA)\d{5,9}|(?:SRS|ERS|DRS)\d{5,9})\b")
QIITA_RE = re.compile(r"\b(?:qiita[ _:-]*(?:study[ _:-]*)?(?:id[ _:-]*)?)(\d{3,6})\b", re.I)

# Titles we want: adult saliva / oral studies. Titles we don't: disease-only,
# animal, in vitro, plaque-only. Screening still decides; this is a cheap prefilter.
GOOD = ["saliva", "salivary", "oral microbiome", "oral microbiota", "mouth",
        "tongue", "oral cavity", "spit", "mouthwash", "oral rinse"]
BAD = ["mouse", "mice", "rat ", "canine", "dog", "lion", "primate", "macaque",
       "in vitro", "denture", "root canal", "endodontic", "calculus", "ancient",
       "cell line", "biofilm model", "pig", "bovine", "cat "]


def http_get(url, tries=4, timeout=180):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "iata-csic-oral-age/1.1"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


# ---------------------------------------------------------------- A) supplements
def read_text_any(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".csv", ".tsv", ".txt", ".md"):
        return open(path, errors="replace").read()
    if ext == ".xlsx":
        import pandas as pd
        return "\n".join(
            df.astype(str).to_csv(index=False)
            for df in pd.read_excel(path, sheet_name=None, header=None, dtype=str).values())
    if ext == ".docx":
        import zipfile
        with zipfile.ZipFile(path) as z:
            return re.sub(r"<[^>]+>", " ", z.read("word/document.xml").decode("utf-8", "replace"))
    if ext == ".pdf":
        try:
            import pdfplumber
            with pdfplumber.open(path) as pdf:
                return "\n".join((p.extract_text() or "") for p in pdf.pages)
        except ImportError:
            import subprocess
            return subprocess.run(["pdftotext", path, "-"], capture_output=True,
                                  text=True).stdout
    raise SystemExit(f"Unsupported supplement format: {ext}")


def mine_supplement(path, source):
    txt = read_text_any(path)
    insdc = sorted(set(ACC_RE.findall(txt)))
    qiita = sorted({f"qiita:{q}" for q in QIITA_RE.findall(txt)})
    print(f"{source}: {len(insdc)} INSDC accessions, {len(qiita)} Qiita study ids "
          f"found in {os.path.basename(path)}", file=sys.stderr)
    return [(a, source) for a in insdc + qiita]


def resolve_to_bioproject(acc):
    """SRP/ERP/SRS style accessions -> their parent BioProject, via ENA."""
    if acc.startswith("PRJ") or acc.startswith("qiita:"):
        return acc
    q = urllib.parse.urlencode({"accession": acc, "result": "read_run",
                                "fields": "study_accession", "format": "tsv", "limit": 1})
    try:
        rows = list(csv.DictReader(http_get(f"{ENA_FILEREPORT}?{q}").splitlines(), delimiter="\t"))
        return rows[0]["study_accession"] if rows else acc
    except Exception:
        return acc


# ---------------------------------------------------------------- B) discovery
def discover_studies():
    q = urllib.parse.urlencode({"result": "study", "query": f"tax_eq({ORAL_TAXON})",
                                "fields": "study_accession,study_title",
                                "format": "tsv", "limit": 0})
    rows = list(csv.DictReader(http_get(f"{ENA_SEARCH}?{q}").splitlines(), delimiter="\t"))
    keep = []
    for r in rows:
        title = (r.get("study_title") or "").lower()
        if any(b in title for b in BAD):
            continue
        if any(g in title for g in GOOD):
            keep.append((r["study_accession"], r.get("study_title", "")))
    print(f"ENA: {len(rows)} studies with human oral metagenome runs; "
          f"{len(keep)} pass the title filter", file=sys.stderr)
    return keep


def run_summary(bp):
    q = urllib.parse.urlencode({"result": "read_run", "format": "tsv", "limit": 0,
                                "fields": "run_accession,sample_accession,library_strategy,"
                                          "instrument_platform,country",
                                "query": f'study_accession="{bp}" AND tax_eq({ORAL_TAXON})'})
    rows = list(csv.DictReader(http_get(f"{ENA_SEARCH}?{q}").splitlines(), delimiter="\t"))
    samples = sorted({r["sample_accession"] for r in rows if r.get("sample_accession")})
    strat = {}
    for r in rows:
        strat[r.get("library_strategy", "")] = strat.get(r.get("library_strategy", ""), 0) + 1
    countries = sorted({r.get("country", "").split(":")[0] for r in rows if r.get("country")})
    return rows, samples, strat, countries


AGE_KEYS = ["host age", "host_age", "age", "age_years", "birth_year", "age_cat",
            "age_group", "host_age_years"]


def screen_age(sample_accs, n=8):
    """Fetch a few BioSamples and report which age-like tags exist, and whether numeric."""
    if not sample_accs:
        return "", "", 0
    chunk = sample_accs[:n]
    try:
        xml = http_get(ENA_XML + ",".join(chunk))
        root = ET.fromstring(xml)
    except Exception:
        return "fetch failed", "", 0
    tags_seen, numeric, checked = set(), 0, 0
    for s in root.iter("SAMPLE"):
        checked += 1
        for a in s.iter("SAMPLE_ATTRIBUTE"):
            tag = (a.findtext("TAG") or "").strip().lower()
            val = (a.findtext("VALUE") or "").strip()
            if tag in AGE_KEYS and val and val.lower() not in (
                    "missing", "not collected", "not applicable", "na", "unknown", "not provided"):
                tags_seen.add(tag)
                if re.fullmatch(r"\d+(\.\d+)?\s*\w*", val) and not re.fullmatch(r"\d+\s*-\s*\d+", val):
                    numeric += 1
    status = ("exact age deposited" if numeric else
              "age-like tag but not numeric" if tags_seen else "no age in repository")
    return status, ";".join(sorted(tags_seen)), checked


# ---------------------------------------------------------------- studies.tsv
def load_studies(path):
    if not os.path.exists(path):
        return []
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    for r in rows:
        for c in COLS:
            r.setdefault(c, "")
    return rows


def save_studies(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def cohort_id_for(bp, title, existing):
    base = re.sub(r"[^A-Za-z0-9]+", "_", (title or bp)).strip("_").upper()[:16] or bp
    cid, i = base, 2
    while cid in existing:
        cid, i = f"{base}_{i}", i + 1
    return cid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", default="metadata/studies.tsv")
    ap.add_argument("--from-supplement")
    ap.add_argument("--source", default="SUPPLEMENT")
    ap.add_argument("--discover", action="store_true")
    ap.add_argument("--no-screen", action="store_true")
    ap.add_argument("--screen", type=int, default=8, help="BioSamples to check per study")
    ap.add_argument("--min-runs", type=int, default=30)
    ap.add_argument("--out", default="metadata/candidate_studies.tsv")
    a = ap.parse_args()

    rows = load_studies(a.studies)
    have = {r["bioproject"] for r in rows}
    ids = {r["cohort_id"] for r in rows}
    candidates, added = [], 0

    if a.from_supplement:
        for acc, src in mine_supplement(a.from_supplement, a.source):
            bp = resolve_to_bioproject(acc)
            if bp in have:
                for r in rows:                       # record that a 2nd source lists it
                    if r["bioproject"] == bp and src not in r["source"]:
                        r["source"] = (r["source"] + "," + src).strip(",")
                continue
            cid = cohort_id_for(bp, "", ids); ids.add(cid); have.add(bp)
            rows.append({"bioproject": bp, "cohort_id": cid, "include": "N", "tax_id": "",
                         "source": src,
                         "notes": "from supplement; verify sample type, age and N before include=Y"})
            added += 1

    if a.discover:
        for bp, title in discover_studies():
            if bp in have:
                continue
            try:
                runs, samples, strat, countries = run_summary(bp)
            except Exception as e:
                print(f"  {bp}: run query failed ({e})", file=sys.stderr)
                continue
            if len(runs) < a.min_runs:
                continue
            if a.no_screen:
                status, tags, checked = "not screened", "", 0
            else:
                status, tags, checked = screen_age(samples, a.screen)
                time.sleep(0.4)
            candidates.append({"bioproject": bp, "study_title": title, "n_runs": len(runs),
                               "n_samples": len(samples), "age_status": status,
                               "age_tags": tags, "screened": checked,
                               "strategies": ";".join(f"{k}:{v}" for k, v in strat.items()),
                               "countries": ";".join(countries[:5])})
            inc = "Y" if status == "exact age deposited" and len(samples) >= a.min_runs else "N"
            cid = cohort_id_for(bp, title, ids); ids.add(cid); have.add(bp)
            rows.append({"bioproject": bp, "cohort_id": cid, "include": inc, "tax_id": str(ORAL_TAXON),
                         "source": "ENA_DISCOVERY",
                         "notes": f"{status}; {len(samples)} oral samples; {title[:70]}"})
            added += 1
            print(f"  {bp:<12} {len(samples):>5} samples  {status:<28} {title[:60]}", file=sys.stderr)

    save_studies(a.studies, rows)
    if candidates:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(candidates[0]), delimiter="\t")
            w.writeheader(); w.writerows(candidates)
    ready = sum(1 for r in rows if r["include"].upper().startswith("Y"))
    print(f"\n{added} study rows added; {len(rows)} total, {ready} marked include=Y.\n"
          f"Review {a.out} and the include column, then run harvest_oral_age_metadata.py.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
