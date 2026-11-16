#!/usr/bin/env bash
# --------------------------------------------------------------------------
# check_group.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# check_group.sh - how did one DADA2 group go?
#   bash scripts/check_group.sh SOUTHKORE_04331__paired__V3-V4
G="${1:?Usage: bash scripts/check_group.sh <group name>}"
T="processed/${G}_track.tsv"
[ -f "$T" ] || { echo "No track file for $G - did the run finish?"; exit 1; }
python3 - "$T" <<'PY'
import sys, pandas as pd
t = pd.read_csv(sys.argv[1], sep="\t")
print(f"samples: {len(t)}")
for c in ["input", "filtered", "denoisedF", "merged", "nonchim"]:
    if c in t and t[c].notna().any():
        print(f"  {c:<10} total {int(t[c].sum()):>12,}   median per sample {int(t[c].median()):>8,}")
pk = t.pct_kept.dropna()
print(f"\npct_kept: median {pk.median():.1f}%  range {pk.min():.1f}-{pk.max():.1f}%")
lost = t[t.nonchim < 1000]
print(f"samples under 1,000 reads: {len(lost)}" + (f"  -> {list(lost.run_accession[:8])}" if len(lost) else ""))
if "merged" in t and t.merged.notna().any():
    drop = 1 - t.merged.sum() / t.denoisedF.sum()
    print(f"lost at the merge step: {drop:.0%}" +
          ("   <-- too high, the pair barely overlaps" if drop > 0.4 else ""))
v = "OK" if pk.median() >= 40 else "LOW - revisit truncLen for this group"
print(f"\nverdict: {v}")
PY
