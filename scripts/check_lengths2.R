# --------------------------------------------------------------------------
# check_lengths2.R
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
b <- sub("_seqtab.rds$", "", basename(f))
reg <- sub(".*__", "", b)
rte <- sapply(strsplit(b, "__"), function(x) x[2])
key <- paste(reg, rte)
L <- list(); A <- list(); G <- list()
for (i in seq_along(f)) {
  st <- readRDS(f[i]); k <- key[i]
  L[[k]] <- c(L[[k]], nchar(colnames(st)))
  A[[k]] <- c(A[[k]], colSums(st))
  G[[k]] <- c(G[[k]], b[i])
}
cat(sprintf("%-20s %5s %7s %6s %8s %9s %12s\n",
    "region x route","grps","ASVs","modal","out ASV%","out read%","range"))
for (k in names(L)[order(names(L))]) {
  n <- L[[k]]; a <- A[[k]]
  md <- as.integer(names(which.max(tapply(a, n, sum))))
  ok <- abs(n - md) <= 30
  cat(sprintf("%-20s %5d %7d %6d %7.2f%% %8.3f%% %6d-%d\n",
      k, length(G[[k]]), length(n), md, 100*mean(!ok), 100*sum(a[!ok])/sum(a), min(n), max(n)))
}
