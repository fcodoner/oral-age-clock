#!/usr/bin/env bash
# --------------------------------------------------------------------------
# finalize_dataset.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# finalize_dataset.sh - Phase 1 salivary age clock (IATA-CSIC)
# Runs the fix -> re-curate -> re-QC loop. Safe to run repeatedly.
# Usage:  bash scripts/finalize_dataset.sh <harvest_dir>
# e.g.    bash scripts/finalize_dataset.sh metadata/harvest_2026-09-24_v2
set -euo pipefail

HARVEST="${1:?Usage: bash scripts/finalize_dataset.sh <harvest_dir>}"
STAMP=$(date +%F_%H%M)
NEW="metadata/harvest_${STAMP}"

echo ">>> 1/3  re-harmonising metadata offline (no downloads)"
python scripts/harvest_oral_age_metadata.py \
       --rebuild "${HARVEST}/manifest_runs.tsv" \
       --out "${NEW}"

echo ">>> 2/3  re-curating healthy adults (core tier)"
python scripts/curate_healthy_adults.py \
       "${NEW}/manifest_runs.tsv" \
       --out metadata/curated_core --tier core --rules metadata/rules.tsv

echo ">>> 3/3  re-running QC"
python scripts/qc_cohorts.py metadata/curated_core \
       --excluded metadata/curated_core/excluded_samples.tsv \
       --out metadata/qc

echo
echo "New harvest dir: ${NEW}"
echo "Pass this path to the next run of this script."
grep -E "FINAL:|participants across" metadata/curated_core/curation_report.txt || true
