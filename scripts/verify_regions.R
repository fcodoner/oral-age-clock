#!/usr/bin/env Rscript
# --------------------------------------------------------------------------
# verify_regions.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# verify_regions.R - check each group's observed amplicon length against the
# other groups annotated with the same hypervariable region.
args <- commandArgs(trailingOnly = TRUE)
params_file <- if (length(args) >= 1) args[1] else "metadata/dada2/dada2_params.tsv"
outdir      <- if (length(args) >= 2) args[2] else "processed"
p <- read.delim(params_file, stringsAsFactors = FALSE)

wq <- function(L, w, q) { o <- order(L); L <- L[o]; w <- w[o]
                          L[which(cumsum(w)/sum(w) >= q)[1]] }

res <- do.call(rbind, lapply(seq_len(nrow(p)), function(i) {
  g <- p$group[i]; f <- file.path(outdir, paste0(g, "_seqtab.rds"))
  if (!file.exists(f)) return(NULL)
  s <- readRDS(f); L <- nchar(colnames(s)); w <- colSums(s)
  tl <- suppressWarnings(as.integer(p$trimLeft[i]));  if (is.na(tl)) tl <- 0
  tr <- suppressWarnings(as.integer(p$trimRight[i])); if (is.na(tr)) tr <- 0
  fo  <- !is.na(p$forward_only[i]) && p$forward_only[i] == "yes"
  lay <- p$layout[i]
  verif <- (lay == "paired" && !fo) || lay == "merged"
  adj <- if (lay == "merged") tl + tr else 2 * tl
  data.frame(group=g, region=p$region[i], layout=lay, fwd_only=ifelse(fo,"yes",""),
             len_p10=wq(L,w,.10), len_med=wq(L,w,.50), len_p90=wq(L,w,.90),
             norm_len=ifelse(verif, wq(L,w,.50)+adj, NA_integer_),
             verifiable=verif, stringsAsFactors=FALSE)
}))

res$flag <- ""
for (rg in unique(res$region)) {
  k <- res$region == rg & res$verifiable
  if (sum(k) < 2) { res$flag[k] <- "sole testable group in region"; next }
  ref <- median(res$norm_len[k]); dev <- res$norm_len[k] - ref
  res$flag[k] <- ifelse(abs(dev) > 25,
    sprintf("CHECK %+d bp vs region median %d", round(dev), round(ref)), "ok")
}
res$flag[!res$verifiable] <- "not testable by length"

res <- res[order(res$region, -res$verifiable, res$norm_len), ]
print(res[, c("group","region","layout","fwd_only","len_p10","len_med","len_p90","norm_len","flag")],
      row.names = FALSE)
cat("\nRegion summary (primer-inclusive equivalent length):\n")
for (rg in unique(res$region)) {
  k <- res$region == rg & res$verifiable
  if (!any(k)) { cat(sprintf("  %-7s  no testable group\n", rg)); next }
  cat(sprintf("  %-7s  n=%2d   median %4d bp   range %4d-%4d\n", rg, sum(k),
      round(median(res$norm_len[k])), min(res$norm_len[k]), max(res$norm_len[k])))
}
cat(sprintf("\n%d of %d groups testable by length; %d flagged.\n",
    sum(res$verifiable), nrow(res), sum(grepl("^CHECK", res$flag))))
write.table(res, "metadata/region_check.tsv", sep="\t", quote=FALSE, row.names=FALSE)
