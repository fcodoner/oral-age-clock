#!/usr/bin/env bash
# --------------------------------------------------------------------------
# run_paired.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# run_paired.sh - Paired-end Illumina groups
#
# Forward and reverse reads are denoised separately, then merged.
# Expect 50-70% kept. Watch the merge step for groups with tight overlap:
# SLL2_ES (32 bp), JAPAN_25734 (32), ZHU20_CN_2 (38), WANG20_CN (39).
# If more than 40% is lost at merging there, relax truncLen before keeping the result.
# The last three are large (309-470 runs) - run them in tmux.
#
# All output goes to processed/, the same folder as every other group.
# Already-finished groups are skipped, so this is safe to re-run.

set -u
P="metadata/dada2/dada2_params.tsv"
O="processed"

for G in \
  JAPAN_11746__paired__V1-V2 \
  CABRAL17_US__paired__V1-V2 \
  CHINA_48248__paired__V4 \
  WOLF17_EU__paired__V4 \
  INDIA_36346__paired__V4 \
  ANTARCTIC_19274__paired__V4 \
  CHINA_43872__paired__V4 \
  HIJAZI20_EU__paired__V1-V2 \
  LINDHEIM16_EU__paired__V1-V2 \
  SPAIN_83278__paired__V3-V4 \
  WANG16_CN__paired__V4 \
  PORTUGAL_74247__paired__V1-V2 \
  WANG20_CN__paired__V3-V4 \
  ZHU20_CN__paired__V3-V4 \
  SPAIN_74796__paired__V1-V2 \
  SOHAIL19_QA__paired__V3-V4 \
  TAIWAN_37501__paired__V3-V4 \
  TAIWAN_39064__paired__V3-V4 \
  RELVAS21_EU__paired__V3-V4 \
  USA_12320__paired__V1-V2 \
  ZHU20_CN_2__paired__V3-V4 \
  JAPAN_25734__paired__V3-V4 \
  JAPAN_17623__paired__V1-V2 \
  USA_30897__paired__V4 \
  KUMPITSCH20_EU__paired__V4 \
  USA_27304__paired__V4 \
  SOUTHKORE_13236__paired__V3-V4 \
  SLL2_ES__paired__V3-V4 \
  TRIAD_ES__paired__V3-V4
do
  echo "=============================================================="
  echo "  $G"
  echo "=============================================================="
  Rscript scripts/run_dada2.R "$G" --params "$P" --out "$O" || \
    { echo "FAILED: $G"; continue; }
  bash scripts/check_group.sh "$G" || true
done

echo "paired groups finished."
