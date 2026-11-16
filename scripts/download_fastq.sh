#!/usr/bin/env bash
# --------------------------------------------------------------------------
# download_fastq.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# download_fastq.sh - resumable download of the final sample set (IATA-CSIC)
#
#   bash scripts/download_fastq.sh metadata/FINAL_accessions_16S.txt raw
#
# Safe to stop and restart at any time: runs already downloaded are skipped.
# Writes a log and a list of failures. Nothing is ever deleted.
set -uo pipefail

LIST="${1:?Usage: bash scripts/download_fastq.sh <accession list> <output folder>}"
OUT="${2:-raw}"
mkdir -p "$OUT"
LOG="$OUT/download.log"
FAILED="$OUT/failed_accessions.txt"
: > "$FAILED"

TOTAL=$(grep -c . "$LIST")
DONE=0; SKIP=0; FAIL=0; I=0

echo "=== started $(date) : $TOTAL accessions ===" | tee -a "$LOG"

while read -r ACC; do
  [ -z "$ACC" ] && continue
  I=$((I+1))

  # already there? (either gzipped or not)
  if ls "$OUT/${ACC}"*.fastq.gz >/dev/null 2>&1; then
    SKIP=$((SKIP+1))
    printf "[%d/%d] %s already downloaded, skipping\n" "$I" "$TOTAL" "$ACC" | tee -a "$LOG"
    continue
  fi

  printf "[%d/%d] %s downloading...\n" "$I" "$TOTAL" "$ACC" | tee -a "$LOG"
  OK=0
  for TRY in 1 2 3; do
    if fasterq-dump --split-files --skip-technical --outdir "$OUT" "$ACC" >>"$LOG" 2>&1; then
      OK=1; break
    fi
    echo "   attempt $TRY failed, waiting 30s" | tee -a "$LOG"
    rm -f "$OUT/${ACC}"*.fastq          # remove any half-written file
    sleep 30
  done

  if [ "$OK" -eq 1 ]; then
    gzip -f "$OUT/${ACC}"*.fastq 2>/dev/null
    DONE=$((DONE+1))
  else
    FAIL=$((FAIL+1))
    echo "$ACC" >> "$FAILED"
    echo "   GIVING UP on $ACC (recorded in failed_accessions.txt)" | tee -a "$LOG"
  fi
done < "$LIST"

echo "=== finished $(date) ===" | tee -a "$LOG"
echo "downloaded: $DONE   already had: $SKIP   failed: $FAIL   of $TOTAL" | tee -a "$LOG"
[ "$FAIL" -gt 0 ] && echo "Retry the failures with:  bash scripts/download_fastq.sh $FAILED $OUT" | tee -a "$LOG"
exit 0
