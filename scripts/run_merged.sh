#!/usr/bin/env bash
# --------------------------------------------------------------------------
# run_merged.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# run_merged.sh - Already-merged amplicons
#
# The submitter merged the read pairs before depositing, so there is no pair-merging
# step. DADA2 still denoises the merged sequences, which is correct.
# maxEE is 5 rather than 2 because these reads are about 460 bp: maxEE 2 over that
# length demands twice the per-base accuracy of a 250 bp read and fails most reads.
#
# All output goes to processed/, the same folder as every other group.
# Already-finished groups are skipped, so this is safe to re-run.

set -u
P="metadata/dada2/dada2_params.tsv"
O="processed"

for G in \
  CHAUDHARI20_IN__merged__V3-V4 \
  CHINA_41273__merged__V3-V4 \
  TIAN17_US__merged__V3-V4 \
  NETHERLAN_46654__merged__V3-V4 \
  SUN20_CN__merged__V3-V4
do
  echo "=============================================================="
  echo "  $G"
  echo "=============================================================="
  Rscript scripts/run_dada2.R "$G" --params "$P" --out "$O" || \
    { echo "FAILED: $G"; continue; }
  bash scripts/check_group.sh "$G" || true
done

echo "merged groups finished."
