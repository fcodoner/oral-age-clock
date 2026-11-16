# --------------------------------------------------------------------------
# predict_sample.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# Salivary oral-health scorer -- standalone. Usage:
#   source("scripts/predict_sample.R")
#   score_sample(counts)   # named numeric vector of genus counts for ONE sample
suppressMessages(library(glmnet))
.B <- readRDS("models/oralclock_v1.rds")
score_sample <- function(counts, verbose=TRUE) {
  g <- .B$genera
  v <- setNames(rep(0, length(g)), g)
  hit <- intersect(names(counts), g); v[hit] <- as.numeric(counts[hit])
  cov <- sum(counts[hit]) / max(sum(counts), 1)
  L <- log(v + 0.5); x <- matrix(L - mean(L), nrow=1, dimnames=list(NULL, g))
  lo <- as.numeric(predict(.B$oral_health$model, x, s="lambda.min"))
  p  <- as.numeric(plogis(predict(.B$oral_health$calib, data.frame(s=lo))))
  oc <- as.numeric(predict(.B$oscc$model, x, s="lambda.min"))
  pc <- as.numeric(plogis(predict(.B$oscc$calib, data.frame(s=oc))))
  ag <- as.numeric(predict(.B$age$model, x, s="lambda.min"))
  band <- if (p < 0.30) "LOW" else if (p < 0.70) "INTERMEDIATE" else "ELEVATED"
  cf <- .B$oral_health$coefs; cf <- cf[names(cf)!="(Intercept)"]
  z  <- (x[1, names(cf)] - .B$ref$mean[names(cf)]) / .B$ref$sd[names(cf)]
  contrib <- sort(cf * z, decreasing=TRUE)
  out <- list(oral_health_prob=p, band=band, oscc_prob=pc,
              age_estimate=ag, age_mae=.B$age$loco_mae, genus_coverage=cov,
              pushing_up=round(head(contrib[contrib>0], 5), 3),
              pushing_down=round(head(rev(contrib[contrib<0]), 5), 3))
  if (verbose) {
    cat(sprintf("Oral-health score : %.3f  [%s]\n", p, band))
    cat(sprintf("OSCC score        : %.3f  (interpret only in high-risk referral)\n", pc))
    cat(sprintf("Age estimate      : %.0f years (+/- %.0f, low precision)\n", ag, .B$age$loco_mae))
    cat(sprintf("Genus coverage    : %.1f%% of reads mapped to model features\n", 100*cov))
    cat("\nRaising the score:\n"); print(out$pushing_up)
    cat("\nLowering the score:\n"); print(out$pushing_down)
  }
  invisible(out)
}
