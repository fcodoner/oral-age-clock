#!/usr/bin/env bash
# --------------------------------------------------------------------------
# run_cutadapt.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
set -euo pipefail
PARAMS=metadata/dada2/dada2_params.tsv
LIST=${1:-metadata/reprocess_list.tsv}
OUT=raw_trimmed
mkdir -p "$OUT" logs/cutadapt metadata/dada2/group_manifests_cut

fwd_primer () { case "$1" in
  V1-V2|V1-V3) echo "AGAGTTTGATCMTGGCTCAG";;
  V3-V4)       echo "CCTACGGGNGGCWGCAG";;
  V4)          echo "GTGYCAGCMGCCGCGGTAA";; esac; }
rev_opts () { case "$1" in
  V1-V2) echo "-G TGCTGCCTCCCGTAGGAGT";;
  V1-V3) echo "-G ATTACCGCGGCTGCTGG";;
  V3-V4) echo "-G GACTACHVGGGTATCTAATCC -G GGACTACNVGGGTWTCTAAT";;
  V4)    echo "-G GGACTACNVGGGTWTCTAAT";; esac; }

tail -n +2 "$LIST" | while IFS=$'\t' read -r G FM RM; do
  [ -z "$G" ] && continue
  REG=$(awk -F'\t' -v g="$G" '$1==g{print $4}'  "$PARAMS")
  LAY=$(awk -F'\t' -v g="$G" '$1==g{print $3}'  "$PARAMS")
  FO=$( awk -F'\t' -v g="$G" '$1==g{print $20}' "$PARAMS")
  MAN="metadata/dada2/group_manifests/${G}.tsv"
  NEW="metadata/dada2/group_manifests_cut/${G}.tsv"
  FP=$(fwd_primer "$REG"); RO=$(rev_opts "$REG")
  mkdir -p "$OUT/$G"; : > "logs/cutadapt/${G}.log"
  printf 'run_accession\tfwd\trev\n' > "$NEW"
  n=0
  tail -n +2 "$MAN" | while IFS=$'\t' read -r ACC F R; do
    [ -z "$ACC" ] && continue
    OF="$OUT/$G/${ACC}_1.fastq.gz"; OR="$OUT/$G/${ACC}_2.fastq.gz"
    if [ "$LAY" = "paired" ] && [ "$FO" != "yes" ] && [ -n "${R:-}" ] && [ -f "$R" ]; then
      cutadapt -g "$FP" $RO -e 0.15 -j 4 \
        -o "$OF" -p "$OR" "$F" "$R" >> "logs/cutadapt/${G}.log" 2>&1
      printf '%s\t%s\t%s\n' "$ACC" "$OF" "$OR" >> "$NEW"
    else
      cutadapt -g "$FP" -e 0.15 -j 4 \
        -o "$OF" "$F" >> "logs/cutadapt/${G}.log" 2>&1
      printf '%s\t%s\t\n' "$ACC" "$OF" >> "$NEW"
    fi
    n=$((n+1))
  done
  echo "finished $G  ($(($(wc -l < "$NEW")-1)) samples)"
done
