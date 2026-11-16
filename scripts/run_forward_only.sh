#!/usr/bin/env bash
# --------------------------------------------------------------------------
# run_forward_only.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# run_forward_only.sh - Forward reads only
#
# After quality truncation the pair could not overlap across a 465 bp amplicon,
# so merging would have returned nothing. These samples contribute the V3 half only.
# Keep them OUT of the ASV-level pooled analysis; use them in the closed-reference track.
#
# All output goes to processed/, the same folder as every other group.
# Already-finished groups are skipped, so this is safe to re-run.

set -u
P="metadata/dada2/dada2_params.tsv"
O="processed"

for G in \
  CHAUDHARI20_IN__paired__V3-V4
do
  echo "=============================================================="
  echo "  $G"
  echo "=============================================================="
  Rscript scripts/run_dada2.R "$G" --params "$P" --out "$O" || \
    { echo "FAILED: $G"; continue; }
  bash scripts/check_group.sh "$G" || true
done

echo "forward_only groups finished."
