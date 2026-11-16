# --------------------------------------------------------------------------
# check_lengths.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
f <- list.files(c("final_split/train","final_split/valid"), pattern="rds$", full.names=TRUE)
reg <- sub(".*__", "", sub("_seqtab.rds$", "", basename(f)))
L <- list(); A <- list()
for (i in seq_along(f)) {
  st <- readRDS(f[i]); r <- reg[i]
  L[[r]] <- c(L[[r]], nchar(colnames(st)))
  A[[r]] <- c(A[[r]], colSums(st))
}
for (r in names(L)) {
  n <- L[[r]]; a <- A[[r]]
  md <- as.integer(names(which.max(tapply(a, n, sum))))
  ok <- abs(n - md) <= 30
  cat(sprintf("%-6s ASVs %7d | modal len %3d | outside +/-30bp: %6d ASVs (%5.2f%% of ASVs, %6.3f%% of reads) | range %d-%d\n",
      r, length(n), md, sum(!ok), 100*mean(!ok), 100*sum(a[!ok])/sum(a), min(n), max(n)))
}
