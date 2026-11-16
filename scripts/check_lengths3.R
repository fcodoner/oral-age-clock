# --------------------------------------------------------------------------
# check_lengths3.R
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
cl <- ifelse(grepl("/train/", f), "T", "V")
o <- order(sub(".*__","",b), sapply(strsplit(b,"__"), `[`, 2), b)
cat(sprintf("%-34s %2s %7s %6s %6s %9s %10s\n",
    "group","cl","ASVs","modal","unifm","out read%","range"))
for (i in o) {
  st <- readRDS(f[i]); n <- nchar(colnames(st)); a <- colSums(st)
  md <- as.integer(names(which.max(tapply(a, n, sum))))
  ok <- abs(n - md) <= 30
  cat(sprintf("%-34s %2s %7d %6d %6s %8.3f%% %6d-%d%s\n",
      b[i], cl[i], length(n), md, if (length(unique(n))==1) "yes" else "no",
      100*sum(a[!ok])/sum(a), min(n), max(n),
      if (100*sum(a[!ok])/sum(a) > 5) "   <-- CHECK" else ""))
}
