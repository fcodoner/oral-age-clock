# --------------------------------------------------------------------------
# salivotypes.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages({library(vegan); library(cluster)}); set.seed(42)
R  <- readRDS("features/genus_relab.rds"); rownames(R) <- unname(rownames(R))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(R), sm$run_accession), ]
R  <- R[sm$class=="train", ]; sm <- sm[sm$class=="train", ]
R  <- R[, colMeans(R>0) > 0.10, drop=FALSE]; R <- R/rowSums(R)
cat(sprintf("healthy samples %d | genera %d | cohorts %d\n",
    nrow(R), ncol(R), length(unique(sm$cohort_id))))

ari <- function(a, b) {
  t <- table(a,b); n <- sum(t)
  s <- sum(choose(t,2)); sa <- sum(choose(rowSums(t),2)); sb <- sum(choose(colSums(t),2))
  e <- sa*sb/choose(n,2); (s-e)/((sa+sb)/2 - e) }

cat("\ncomputing Bray-Curtis...\n"); D <- vegdist(R, method="bray")
cat(sprintf("\n%4s %10s %10s %10s %10s\n","k","avgSil","ARI cohort","ARI region","ARI route"))
best <- NULL
for (k in 2:8) {
  p <- pam(D, k, diss=TRUE)
  s <- p$silinfo$avg.width
  cat(sprintf("%4d %10.3f %10.3f %10.3f %10.3f\n", k, s,
      ari(p$clustering, sm$cohort_id), ari(p$clustering, sm$region), ari(p$clustering, sm$route)))
  if (is.null(best) || s > best$s) best <- list(k=k, s=s, p=p)
}
K <- best$k; cl <- best$p$clustering
cat(sprintf("\nbest k = %d (avg silhouette %.3f)\n", K, best$s))

cat("\n--- is each type present across many cohorts? ---\n")
tb <- table(sm$cohort_id, cl)
pr <- sweep(tb, 1, rowSums(tb), "/")
cat(sprintf("%6s %7s %9s %12s %10s\n","type","n","% samples","cohorts>=10%","medianAge"))
for (j in 1:K)
  cat(sprintf("%6d %7d %8.1f%% %12d %10.0f\n", j, sum(cl==j), 100*mean(cl==j),
      sum(pr[,j] >= 0.10), median(sm$age[cl==j])))
cat("\n(a biological type appears in most cohorts; a batch artefact appears in few)\n")

cat("\n--- driver genera per type (mean % abundance) ---\n")
mu <- t(sapply(1:K, function(j) colMeans(R[cl==j,,drop=FALSE])))
for (j in 1:K) {
  enr <- sort(mu[j,] / (colMeans(R) + 1e-6), decreasing=TRUE)[1:5]
  cat(sprintf("type %d: %s\n", j,
      paste(sprintf("%s %.1f%% (x%.1f)", names(enr), 100*mu[j,names(enr)], enr), collapse=", ")))
}
saveRDS(list(k=K, clustering=cl, pam=best$p, samples=rownames(R)), "features/salivotypes.rds")
cat("\nsaved: features/salivotypes.rds\n")
