# --------------------------------------------------------------------------
# picrust_qc.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
tre <- readLines("picrust/placed.tre", warn=FALSE)
tips <- unique(regmatches(tre, gregexpr("ASV[0-9]+", tre))[[1]])
map  <- read.delim("picrust/asv_map.tsv", stringsAsFactors=FALSE)
lost <- setdiff(map$asv_id, tips)
cat(sprintf("ASVs submitted %d | placed %d | excluded %d (%.1f%%)\n",
    nrow(map), length(intersect(map$asv_id, tips)), length(lost),
    100*length(lost)/nrow(map)))
cat(sprintf("excluded carry %.2f%% of the reads in the table\n",
    100*sum(map$total_reads[map$asv_id %in% lost])/sum(map$total_reads)))

A  <- read.delim("picrust/asv_table.tsv", row.names=1, check.names=FALSE)
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(colnames(A), sm$run_accession), ]
keep <- rownames(A) %in% tips
fr <- 100*colSums(A[keep,,drop=FALSE])/colSums(A)
q  <- tapply(fr, sm$cohort_id, median)
cat(sprintf("\nreads surviving placement: median %.1f%% | range %.1f%% - %.1f%%\n",
    median(fr), min(q), max(q)))
cat("lowest cohorts:\n"); print(round(sort(q)[1:6], 1))
