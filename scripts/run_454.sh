#!/usr/bin/env bash
# --------------------------------------------------------------------------
# run_454.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# run_454.sh - 454 pyrosequencing groups (Japanese cohorts)
#
# Reads vary in length, so they are NOT truncated: minLen does the filtering.
# DADA2 is given its long-read settings (HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32).
# 454 has a higher error rate than Illumina: 30-50% of reads kept is normal here,
# lower than the 50-70% you should see elsewhere. Under 20% means something is wrong.
#
# All output goes to processed/, the same folder as every other group.
# Already-finished groups are skipped, so this is safe to re-run.

set -u
P="metadata/dada2/dada2_params.tsv"
O="processed"

for G in \
  JAPAN_10985__single__V1-V2 \
  JAPAN_10651__single__V1-V2 \
  JAPAN_10985__single_long__V1-V2 \
  JAPAN_10651__single_long__V1-V2
do
  echo "=============================================================="
  echo "  $G"
  echo "=============================================================="
  Rscript scripts/run_dada2.R "$G" --params "$P" --out "$O" || \
    { echo "FAILED: $G"; continue; }
  bash scripts/check_group.sh "$G" || true
done

echo "454 groups finished."
