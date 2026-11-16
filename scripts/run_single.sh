#!/usr/bin/env bash
# --------------------------------------------------------------------------
# run_single.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# run_single.sh - single-end Illumina groups. Output goes to processed/.
# Finished groups are skipped, so this is safe to re-run.
set -u
P="metadata/dada2/dada2_params.tsv"
O="processed"

for G in \
  JAPAN_10985__single__V1-V2 \
  USA_29169__single__V4 \
  CHINA_48248__single__V4 \
  JAPAN_10651__single__V1-V2 \
  YU20_CN__single__V4 \
  USA_48040__single__V4 \
  CHEN18_US__single__V4 \
  USA_74379__single__V1-V3 \
  USA_39492__single__V1-V3 \
  ISRAEL_38062__single__V4 \
  LEE17_CN__single__V4 \
  AGP_ORAL__single__V4
do
  echo "=== $G ==="
  Rscript scripts/run_dada2.R "$G" --params "$P" --out "$O" || { echo "FAILED: $G"; continue; }
  bash scripts/check_group.sh "$G" || true
done
echo "single-end groups finished."
