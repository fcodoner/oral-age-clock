# --------------------------------------------------------------------------
# matched_comparison.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages(library(glmnet))
src <- readLines("scripts/reported_values.R")
i <- grep("^folds *<-", src)[1]; stopifnot(!is.na(i))
eval(parse(text=paste(src[1:i], collapse="\n")), envir=globalenv())
cat("eligible folds:", length(folds), "\n")

keep <- co %in% folds
Xs <- XG[tr,][keep,,drop=FALSE]; ys <- y[keep]; cs <- co[keep]
cat("matched set:", nrow(Xs), "samples |", length(unique(cs)), "cohorts\n\n")

mae <- function(p,a) mean(abs(p-a))
PL <- rep(NA_real_, length(ys))
for (f in folds) { te <- cs==f
  cv <- fit(Xs[!te,,drop=FALSE], ys[!te])
  PL[te] <- as.numeric(predict(cv, Xs[te,,drop=FALSE], s="lambda.min")) }
cat(sprintf("LOCO       MAE %.3f  r %.3f   (should match 12.003 / 0.469)\n", mae(PL,ys), cor(PL,ys)))

set.seed(1); kk <- sample(FF(length(ys))); PR <- rep(NA_real_, length(ys))
for (j in 1:5) { cv <- fit(Xs[kk!=j,,drop=FALSE], ys[kk!=j])
  PR[kk==j] <- as.numeric(predict(cv, Xs[kk==j,,drop=FALSE], s="lambda.min")) }
cat(sprintf("randomCV   MAE %.3f  r %.3f   (SAME samples, same cohorts)\n", mae(PR,ys), cor(PR,ys)))
cat(sprintf("difference %.3f years\n\n", mae(PL,ys)-mae(PR,ys)))

CU <- unique(cs); set.seed(7); N <- 2000; bg <- bl <- bn <- numeric(N)
for (b in 1:N) { s <- sample(CU, length(CU), replace=TRUE)
  ix <- unlist(lapply(s, function(k) which(cs==k)))
  bl[b] <- mae(PL[ix], ys[ix]); bg[b] <- mae(PL[ix],ys[ix]) - mae(PR[ix],ys[ix])
  bn[b] <- mean(sapply(s, function(k) mean(abs(ys[cs==k]-mean(ys[cs!=k]))))) }
q <- function(v) sprintf("%.2f [%.2f, %.2f]", mean(v), quantile(v,.025), quantile(v,.975))
cat("-- paired cluster bootstrap, identical samples both arms --\n")
cat("LOCO MAE             ", q(bl), "\n")
cat("improvement vs null  ", q(bn-bl), "\n")
cat("LOCO - randomCV      ", q(bg), "\n")
saveRDS(list(loco=PL, rand=PR, y=ys, cohort=cs, boot=list(gap=bg,loco=bl,null=bn)),
        "results/matched_comparison.rds")
cat("\nsaved results/matched_comparison.rds\n")
