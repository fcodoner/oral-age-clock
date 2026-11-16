# --------------------------------------------------------------------------
# sizing_estimator.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages(library(glmnet)); set.seed(42)
M  <- readRDS("features/genus_counts.rds"); rownames(M) <- unname(rownames(M))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(M), sm$run_accession), ]
B  <- readRDS("models/oralclock_final.rds")
L  <- log(M[,B$genera,drop=FALSE]+0.5); X <- L - rowMeans(L)
sc <- as.numeric(predict(B$model, X, s="lambda.min"))
Z <- qnorm(0.90)
est <- list(
  empirical  = function(s) quantile(s, 0.90),
  parametric = function(s) mean(s) + Z*sd(s),
  robust     = function(s) median(s) + Z*mad(s))

for (co in c("TRIAD_ES","SLL2_ES","SOUTHKORE_13236")) {
  H <- sc[which(sm$class=="train" & sm$cohort_id==co)]
  cat(sprintf("\n=== %s (%d healthy; normality p = %.3g) ===\n", co, length(H),
      shapiro.test(H[1:min(5000,length(H))])$p.value))
  cat(sprintf("%6s %26s %26s %26s\n","n_ref","empirical","parametric","robust"))
  for (n in c(10,15,20,30,50)) {
    out <- sapply(names(est), function(e) {
      f <- replicate(1000, { i <- sample(length(H), n)
             mean(H[-i] >= est[[e]](H[i])) })
      sprintf("%.1f%% (sd %.1f%%)", 100*mean(f), 100*sd(f)) })
    cat(sprintf("%6d %26s %26s %26s\n", n, out[1], out[2], out[3]))
  }
}
cat("\nTarget: 10% with the smallest spread. CLSI EP28-A3 requires n>=120 for a direct\n")
cat("empirical reference interval; parametric estimation is standard at smaller n.\n")
