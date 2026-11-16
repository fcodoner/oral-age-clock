#!/usr/bin/env bash
# --------------------------------------------------------------------------
# run_all_groups.sh
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# run_all_groups.sh - one DADA2 command per group, smallest and most unusual first.
# Run the lines ONE AT A TIME, checking pct_kept after each. Or run the whole file
# with:  bash scripts/run_all_groups.sh   (it skips groups already done).
set -u
P="metadata/dada2/dada2_params.tsv"
O="processed"


# ---- 454 short-read groups: no truncation, minLen 150. CHECK THESE FIRST ----
Rscript scripts/run_dada2.R JAPAN_10985__single__V1-V2 --params "$P" --out "$O"   # 8 runs, V1-V2
Rscript scripts/run_dada2.R JAPAN_10651__single__V1-V2 --params "$P" --out "$O"   # 21 runs, V1-V2

# ---- 454 long-read groups: no truncation, minLen 200 ----
Rscript scripts/run_dada2.R JAPAN_10985__single_long__V1-V2 --params "$P" --out "$O"   # 13 runs, V1-V2
Rscript scripts/run_dada2.R JAPAN_10651__single_long__V1-V2 --params "$P" --out "$O"   # 71 runs, V1-V2

# ---- forward reads only (the pair cannot overlap) ----
Rscript scripts/run_dada2.R CHAUDHARI20_IN__paired__V3-V4 --params "$P" --out "$O"   # 24 runs, V3-V4

# ---- already-merged amplicons: no denoising of pairs ----
Rscript scripts/run_dada2.R CHAUDHARI20_IN__merged__V3-V4 --params "$P" --out "$O"   # 6 runs, V3-V4
Rscript scripts/run_dada2.R CHINA_41273__merged__V3-V4 --params "$P" --out "$O"   # 36 runs, V3-V4
Rscript scripts/run_dada2.R TIAN17_US__merged__V3-V4 --params "$P" --out "$O"   # 62 runs, V3-V4
Rscript scripts/run_dada2.R NETHERLAN_46654__merged__V3-V4 --params "$P" --out "$O"   # 73 runs, V3-V4
Rscript scripts/run_dada2.R SUN20_CN__merged__V3-V4 --params "$P" --out "$O"   # 130 runs, V3-V4

# ---- single-end ----
Rscript scripts/run_dada2.R ANTARCTIC_19274__single__V4 --params "$P" --out "$O"   # 1 runs, V4
Rscript scripts/run_dada2.R KUMPITSCH20_EU__single__V4 --params "$P" --out "$O"   # 1 runs, V4
Rscript scripts/run_dada2.R USA_29169__single__V4 --params "$P" --out "$O"   # 9 runs, V4
Rscript scripts/run_dada2.R CHINA_48248__single__V4 --params "$P" --out "$O"   # 10 runs, V4
Rscript scripts/run_dada2.R YU20_CN__single__V4 --params "$P" --out "$O"   # 40 runs, V4
Rscript scripts/run_dada2.R USA_48040__single__V4 --params "$P" --out "$O"   # 56 runs, V4
Rscript scripts/run_dada2.R CHEN18_US__single__V4 --params "$P" --out "$O"   # 81 runs, V4
Rscript scripts/run_dada2.R USA_74379__single__V1-V3 --params "$P" --out "$O"   # 127 runs, V1-V3
Rscript scripts/run_dada2.R USA_39492__single__V1-V3 --params "$P" --out "$O"   # 143 runs, V1-V3
Rscript scripts/run_dada2.R ISRAEL_38062__single__V4 --params "$P" --out "$O"   # 181 runs, V4
Rscript scripts/run_dada2.R LEE17_CN__single__V4 --params "$P" --out "$O"   # 376 runs, V4
Rscript scripts/run_dada2.R AGP_ORAL__single__V4 --params "$P" --out "$O"   # 838 runs, V4

# ---- paired-end ----
Rscript scripts/run_dada2.R JAPAN_11746__paired__V1-V2 --params "$P" --out "$O"   # 5 runs, V1-V2
Rscript scripts/run_dada2.R CABRAL17_US__paired__V1-V2 --params "$P" --out "$O"   # 12 runs, V1-V2
Rscript scripts/run_dada2.R CHINA_48248__paired__V4 --params "$P" --out "$O"   # 20 runs, V4
Rscript scripts/run_dada2.R WOLF17_EU__paired__V4 --params "$P" --out "$O"   # 20 runs, V4
Rscript scripts/run_dada2.R INDIA_36346__paired__V4 --params "$P" --out "$O"   # 31 runs, V4
Rscript scripts/run_dada2.R ANTARCTIC_19274__paired__V4 --params "$P" --out "$O"   # 35 runs, V4
Rscript scripts/run_dada2.R CHINA_43872__paired__V4 --params "$P" --out "$O"   # 37 runs, V4
Rscript scripts/run_dada2.R HIJAZI20_EU__paired__V1-V2 --params "$P" --out "$O"   # 43 runs, V1-V2
Rscript scripts/run_dada2.R LINDHEIM16_EU__paired__V1-V2 --params "$P" --out "$O"   # 44 runs, V1-V2
Rscript scripts/run_dada2.R SPAIN_83278__paired__V3-V4 --params "$P" --out "$O"   # 50 runs, V3-V4
Rscript scripts/run_dada2.R WANG16_CN__paired__V4 --params "$P" --out "$O"   # 55 runs, V4
Rscript scripts/run_dada2.R PORTUGAL_74247__paired__V1-V2 --params "$P" --out "$O"   # 58 runs, V1-V2
Rscript scripts/run_dada2.R WANG20_CN__paired__V3-V4 --params "$P" --out "$O"   # 60 runs, V3-V4
Rscript scripts/run_dada2.R ZHU20_CN__paired__V3-V4 --params "$P" --out "$O"   # 60 runs, V3-V4
Rscript scripts/run_dada2.R SPAIN_74796__paired__V1-V2 --params "$P" --out "$O"   # 66 runs, V1-V2
Rscript scripts/run_dada2.R SOHAIL19_QA__paired__V3-V4 --params "$P" --out "$O"   # 68 runs, V3-V4
Rscript scripts/run_dada2.R TAIWAN_37501__paired__V3-V4 --params "$P" --out "$O"   # 73 runs, V3-V4
Rscript scripts/run_dada2.R TAIWAN_39064__paired__V3-V4 --params "$P" --out "$O"   # 75 runs, V3-V4
Rscript scripts/run_dada2.R RELVAS21_EU__paired__V3-V4 --params "$P" --out "$O"   # 77 runs, V3-V4
Rscript scripts/run_dada2.R USA_12320__paired__V1-V2 --params "$P" --out "$O"   # 78 runs, V1-V2
Rscript scripts/run_dada2.R ZHU20_CN_2__paired__V3-V4 --params "$P" --out "$O"   # 120 runs, V3-V4
Rscript scripts/run_dada2.R JAPAN_25734__paired__V3-V4 --params "$P" --out "$O"   # 139 runs, V3-V4
Rscript scripts/run_dada2.R JAPAN_17623__paired__V1-V2 --params "$P" --out "$O"   # 140 runs, V1-V2
Rscript scripts/run_dada2.R USA_30897__paired__V4 --params "$P" --out "$O"   # 140 runs, V4
Rscript scripts/run_dada2.R KUMPITSCH20_EU__paired__V4 --params "$P" --out "$O"   # 155 runs, V4
Rscript scripts/run_dada2.R USA_27304__paired__V4 --params "$P" --out "$O"   # 215 runs, V4
Rscript scripts/run_dada2.R SOUTHKORE_13236__paired__V3-V4 --params "$P" --out "$O"   # 309 runs, V3-V4
Rscript scripts/run_dada2.R SLL2_ES__paired__V3-V4 --params "$P" --out "$O"   # 446 runs, V3-V4
Rscript scripts/run_dada2.R TRIAD_ES__paired__V3-V4 --params "$P" --out "$O"   # 470 runs, V3-V4

# ---- already processed ----
# Rscript scripts/run_dada2.R SOUTHKORE_04331__paired__V3-V4 --params "$P" --out "$O"   # 32 runs, V3-V4

# after each group, check how many reads survived:
#   awk -F'\t' 'NR>1{s+=$7; n++} END{print "median-ish pct_kept:", s/n}' processed/<GROUP>_track.tsv
