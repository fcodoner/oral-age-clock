#!/usr/bin/env bash
# --------------------------------------------------------------------------
# download_from_ena.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# download_from_ena.sh - resumable FASTQ download straight from the European
# Nucleotide Archive, using the URLs and checksums already stored in the
# final sample table. Needs only curl, which macOS and Linux already have.
#
#   bash scripts/download_from_ena.sh metadata/FINAL_sample_set_M1.tsv raw
#
# Safe to stop and restart: finished files are verified and skipped.
set -uo pipefail

TABLE="${1:?Usage: bash scripts/download_from_ena.sh <sample table.tsv> <output folder>}"
OUT="${2:-raw}"
mkdir -p "$OUT"
LOG="$OUT/download.log"
FAILED="$OUT/failed_accessions.txt"
: > "$FAILED"

# md5 tool differs between macOS and Linux
if command -v md5sum >/dev/null 2>&1; then MD5() { md5sum "$1" | cut -d' ' -f1; }
elif command -v md5  >/dev/null 2>&1; then MD5() { md5 -q "$1"; }
else echo "No md5 tool found."; exit 1; fi

# extract just the three columns we need, with awk, which preserves empty cells
# (bash's read collapses repeated tabs, which silently shifts every column)
TMP="$OUT/.fields.txt"
awk -F'\t' -v OFS='|' '
  NR==1 {
    for (i=1; i<=NF; i++) h[$i]=i
    if (!("run_accession" in h) || !("fastq_ftp" in h)) {
      print "MISSINGCOLS" > "/dev/stderr"; exit 1
    }
    a=h["run_accession"]; f=h["fastq_ftp"]; m=h["fastq_md5"]
    next
  }
  { print $a, $f, (m ? $m : "") }
' "$TABLE" > "$TMP" || { echo "Could not find run_accession / fastq_ftp columns in $TABLE"; exit 1; }

TOTAL=$(wc -l < "$TMP")
I=0; GOT=0; SKIP=0; FAIL=0
echo "=== started $(date) : $TOTAL runs ===" | tee -a "$LOG"

# IFS='|' is a non-whitespace separator, so empty fields are preserved
while IFS='|' read -r ACC FTP SUMS; do
  I=$((I+1))
  [ -z "$ACC" ] && continue
  if [ -z "$FTP" ] || [ "$FTP" = "nan" ]; then
    echo "[$I/$TOTAL] $ACC has no FTP link recorded" | tee -a "$LOG"
    echo "$ACC" >> "$FAILED"; FAIL=$((FAIL+1)); continue
  fi

  IFS=';' read -ra URLS <<< "$FTP"
  IFS=';' read -ra MD5S <<< "$SUMS"
  OK=1
  for k in "${!URLS[@]}"; do
    URL="${URLS[$k]}"; [ -z "$URL" ] && continue
    # only a real scheme counts: a bare "ftp.sra.ebi.ac.uk/..." must become https,
    # otherwise curl guesses the FTP protocol and connects on port 21
    case "$URL" in
      http://*|https://*|ftp://*) ;;
      *) URL="https://$URL" ;;
    esac
    FILE="$OUT/$(basename "$URL")"
    WANT="${MD5S[$k]:-}"

    if [ -s "$FILE" ] && [ -n "$WANT" ] && [ "$(MD5 "$FILE")" = "$WANT" ]; then
      SKIP=$((SKIP+1)); continue
    fi

    # This network refuses plain GETs and open-ended ranges, but serves an explicit
    # byte range fine. A range ending far beyond the file simply returns the whole
    # file, so no size lookup (and no extra connection) is needed.
    GOOD=0
    for TRY in 1 2 3 4 5; do
      rm -f "$FILE"
      curl -sS -4 --connect-timeout 45 -r 0-99999999999 \
           -o "$FILE" "$URL" >>"$LOG" 2>&1
      RC=$?
      [ "$RC" -ne 0 ] && echo "   curl exit $RC on attempt $TRY for $(basename "$FILE")" >>"$LOG"
      if [ -s "$FILE" ]; then
        if [ -z "$WANT" ] || [ "$(MD5 "$FILE")" = "$WANT" ]; then GOOD=1; break; fi
        echo "   checksum mismatch on $(basename "$FILE"), retrying" | tee -a "$LOG"
        rm -f "$FILE"
      fi
      sleep $((TRY * 10))
    done
    [ "$GOOD" -eq 1 ] || OK=0
  done

  sleep 1
  if [ "$OK" -eq 1 ]; then
    GOT=$((GOT+1))
    printf "[%d/%d] %s ok\n" "$I" "$TOTAL" "$ACC" | tee -a "$LOG"
  else
    FAIL=$((FAIL+1)); echo "$ACC" >> "$FAILED"
    printf "[%d/%d] %s FAILED (recorded in failed_accessions.txt)\n" "$I" "$TOTAL" "$ACC" | tee -a "$LOG"
  fi
done < "$TMP"

rm -f "$TMP"
echo "=== finished $(date) ===" | tee -a "$LOG"
echo "files already verified and skipped, plus newly downloaded runs, are in $OUT" | tee -a "$LOG"
echo "any failures are listed in $FAILED" | tee -a "$LOG"
exit 0
