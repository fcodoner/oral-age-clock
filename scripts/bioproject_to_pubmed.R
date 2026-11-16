# --------------------------------------------------------------------------
# bioproject_to_pubmed.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
A <- read.delim("metadata/studies.tsv", stringsAsFactors=FALSE)
S <- read.csv("results/supp_samples.csv", stringsAsFactors=FALSE)
use <- sort(unique(S$cohort_id))
A <- A[A$cohort_id %in% use & nzchar(A$bioproject), ]
A <- A[!grepl("^qiita", A$bioproject), ]
E <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
grab <- function(u) tryCatch(paste(readLines(u, warn=FALSE), collapse=""), error=function(e) "")
out <- list()
for (i in seq_len(nrow(A))) {
  acc <- A$bioproject[i]; cid <- A$cohort_id[i]; pm <- NA
  x <- grab(sprintf("%s/esearch.fcgi?db=bioproject&term=%s&retmode=xml", E, acc))
  uid <- sub(".*<IdList>.*?<Id>([0-9]+)</Id>.*", "\\1", x)
  if (grepl("^[0-9]+$", uid)) {
    Sys.sleep(0.4)
    y <- grab(sprintf("%s/elink.fcgi?dbfrom=bioproject&db=pubmed&id=%s&retmode=xml", E, uid))
    if (grepl("<LinkSetDb>", y)) {
      seg <- sub(".*<LinkSetDb>(.*?)</LinkSetDb>.*", "\\1", y)
      ids <- regmatches(seg, gregexpr("<Id>([0-9]+)</Id>", seg))[[1]]
      ids <- gsub("[^0-9]", "", ids)
      if (length(ids)) pm <- paste(ids, collapse=";")
    }
  }
  out[[i]] <- data.frame(cohort_id=cid, bioproject=acc, pmids=pm, stringsAsFactors=FALSE)
  cat(sprintf("[%2d/%d] %-20s %-16s %s\n", i, nrow(A), cid, acc,
      ifelse(is.na(pm), "no linked publication", pm)))
  Sys.sleep(0.4)
}
B <- do.call(rbind, out)
write.table(B, "results/cohort_pubmed.tsv", sep="\t", row.names=FALSE, quote=FALSE, na="")
cat("\nwrote results/cohort_pubmed.tsv |", sum(!is.na(B$pmids)), "of", nrow(B), "linked to a publication\n")
