# --------------------------------------------------------------------------
# score.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# Salivary oral-health scorer (final). Usage:
#   source("scripts/score.R")
#   ref <- lab_reference(control_counts)   # matrix: healthy controls x genus counts
#   report(patient_counts, ref)
suppressMessages(library(glmnet))
.B <- readRDS("models/oralclock_final.rds")
.clr <- function(cts) { g <- .B$genera
  v <- setNames(rep(0, length(g)), g)
  h <- intersect(names(cts), g); v[h] <- as.numeric(cts[h])
  L <- log(v + 0.5); list(x=matrix(L-mean(L), 1, length(g), dimnames=list(NULL,g)),
                          cov=sum(cts[h])/max(sum(cts),1)) }
.raw <- function(cts) as.numeric(predict(.B$model, .clr(cts)$x, s="lambda.min"))
lab_reference <- function(control_counts) {
  s <- apply(control_counts, 1, .raw)
  if (length(s) < 15) warning("fewer than 15 controls: reference unreliable; 30+ advised")
  list(scores=s, n=length(s), q50=quantile(s,.50), q90=quantile(s,.90)) }
report <- function(cts, ref=NULL) {
  o <- .clr(cts); s <- as.numeric(predict(.B$model, o$x, s="lambda.min"))
  R <- if (is.null(ref)) list(scores=.B$global_ref, n=.B$n_train,
         q50=quantile(.B$global_ref,.50), q90=quantile(.B$global_ref,.90)) else ref
  pct  <- mean(R$scores <= s)
  band <- if (s >= R$q90) "ELEVATED" else if (s >= R$q50) "INTERMEDIATE" else "LOW"
  z  <- (o$x[1,] - .B$ref$mean) / .B$ref$sd
  ax <- c(PC2=sum(z * .B$loadings[,1]), PC3=sum(z * .B$loadings[,2]))
  cf <- .B$coefs[names(.B$coefs)!="(Intercept)"]
  dr <- sort(cf * z[names(cf)], decreasing=TRUE)
  cat(sprintf("Oral-health band   : %s  (%.0fth percentile of %d reference controls)\n",
      band, 100*pct, R$n))
  if (is.null(ref)) cat("  NOTE: no local reference supplied -- using the global healthy pool.\n",
                        "       Scores are not comparable across laboratories; supply local controls.\n")
  cat(sprintf("Genus coverage     : %.1f%% of reads\n", 100*o$cov))
  if (o$cov < 0.80) cat("  WARNING: low coverage -- score unreliable, do not report.\n")
  cat(sprintf("\nPeriodontal pathogen axis : %+.2f\nCommensal / aerobe axis   : %+.2f\n",
      ax[1], ax[2]))
  cat("\nRaising this score:\n"); print(round(head(dr[dr>0],5),3))
  cat("\nLowering this score:\n"); print(round(head(rev(dr[dr<0]),5),3))
  invisible(list(score=s, percentile=pct, band=band, axes=ax, coverage=o$cov,
                 drivers_up=head(dr[dr>0],5), drivers_down=head(rev(dr[dr<0]),5))) }
