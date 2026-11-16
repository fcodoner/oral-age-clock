# --------------------------------------------------------------------------
# host_and_depth.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
tx <- readRDS("taxonomy/taxa_oriented.rds")$tax
bact <- !is.na(tx[,"Phylum"]) & !is.na(tx[,"Kingdom"]) & tx[,"Kingdom"]=="Bacteria"
names(bact) <- rownames(tx)
MIN <- 3000

f  <- list.files(c("final_split/train","final_split/valid"), pattern="rds$", full.names=TRUE)
b  <- sub("_seqtab.rds$","",basename(f)); cl <- ifelse(grepl("/train/",f),"T","V")
o  <- order(-sapply(seq_along(f), function(i) 1))  # keep file order
cat(sprintf("%-34s %2s %7s %9s %9s %9s %7s\n",
    "group","cl","n","host%","medBactDep","minBact","lost"))
tot <- c(n=0,lost=0); hostf <- c()
for (i in seq_along(f)) {
  st <- readRDS(f[i]); s <- colnames(st)
  bd <- rowSums(st[, bact[s], drop=FALSE]); td <- rowSums(st)
  h <- 100*(1 - sum(bd)/sum(td)); hostf <- c(hostf, setNames(h, b[i]))
  lo <- bd < MIN
  cat(sprintf("%-34s %2s %7d %8.2f%% %9.0f %9.0f %7d%s\n",
      b[i], cl[i], nrow(st), h, median(bd), min(bd), sum(lo),
      if (sum(lo)>0) "  <--" else ""))
  tot <- tot + c(nrow(st), sum(lo))
}
cat(sprintf("\nsamples %d | would be dropped on bacterial depth <%d: %d\n",
    tot[["n"]], MIN, tot[["lost"]]))
cat(sprintf("host fraction: median %.2f%% | worst: %s\n", median(hostf),
    paste(sprintf("%s %.1f%%", names(sort(hostf,decreasing=TRUE))[1:5],
          sort(hostf,decreasing=TRUE)[1:5]), collapse=", ")))
