# --------------------------------------------------------------------------
# diag_build.R
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
MIN <- 3000
lab <- tx[,"Genus"]
for (r in c("Family","Order","Class","Phylum")) {
  i <- is.na(lab) & !is.na(tx[,r]); lab[i] <- paste0(tx[i,r], " (unclassified)")
}
keep <- !is.na(tx[,"Phylum"]) & !is.na(tx[,"Kingdom"]) & tx[,"Kingdom"]=="Bacteria"
names(lab) <- names(keep) <- rownames(tx)

f <- list.files(c("final_split/train","final_split/valid"), pattern="rds$", full.names=TRUE)
tot <- 0L
cat(sprintf("%-34s %6s %6s %6s %6s %6s\n","group","in","bactASV","afterD","out","LOST"))
for (i in seq_along(f)) {
  g <- sub("_seqtab.rds$","",basename(f[i]))
  st <- readRDS(f[i]); n0 <- nrow(st); s <- colnames(st)
  nk <- sum(keep[s])
  st2 <- st[, keep[s], drop=FALSE]
  bd <- rowSums(st2); n1 <- sum(bd >= MIN)
  st3 <- st2[bd >= MIN, , drop=FALSE]
  n2 <- if (nrow(st3)) nrow(t(rowsum(t(st3), lab[colnames(st3)]))) else 0L
  tot <- tot + n2
  cat(sprintf("%-34s %6d %6d %6d %6d %6d%s\n", g, n0, nk, n1, n2, n0-n2,
      if (n0-n2 > 0) "  <--" else ""))
}
cat("\ntotal out:", tot, "\n")
