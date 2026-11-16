# --------------------------------------------------------------------------
# loco_benchmark.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages({library(glmnet); library(ranger)})
set.seed(42)

M  <- readRDS("features/genus_counts.rds"); rownames(M) <- unname(rownames(M))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(M), sm$run_accession), ]
M  <- M[sm$class=="train", ]; sm <- sm[sm$class=="train", ]
M  <- M[, colMeans(M>0) > 0.10, drop=FALSE]

L <- log(M + 0.5)
X <- L - rowMeans(L)                      # centred log-ratio
y <- sm$age; co <- sm$cohort_id
cat(sprintf("train %d samples | %d genera | age mean %.1f sd %.1f\n",
    nrow(X), ncol(X), mean(y), sd(y)))

ag <- tapply(y, co, function(v) diff(range(v))); nn <- table(co)
folds <- names(ag)[ag >= 40 & nn[names(ag)] >= 40]
cat("LOCO folds:", length(folds), "\n\n")

cat(sprintf("%-18s %5s %7s %7s %7s %7s %7s\n",
    "held-out cohort","n","nullMAE","enMAE","rfMAE","en_r","rf_r"))
rows <- list()
for (f in folds) {
  te <- co == f; tr <- !te
  mu <- mean(y[tr])
  nm <- mean(abs(y[te] - mu))
  cvf <- cv.glmnet(X[tr,], y[tr], alpha=0.5, nfolds=5, standardize=TRUE)
  pe  <- as.numeric(predict(cvf, X[te,], s="lambda.min"))
  rf  <- ranger(x=X[tr,], y=y[tr], num.trees=500, num.threads=3, seed=42)
  pr  <- predict(rf, X[te,])$predictions
  r1 <- suppressWarnings(cor(pe, y[te])); r2 <- suppressWarnings(cor(pr, y[te]))
  cat(sprintf("%-18s %5d %7.2f %7.2f %7.2f %7.3f %7.3f\n",
      f, sum(te), nm, mean(abs(pe-y[te])), mean(abs(pr-y[te])), r1, r2))
  rows[[f]] <- data.frame(cohort=f, n=sum(te), age=y[te], null=mu, en=pe, rf=pr)
}
D <- do.call(rbind, rows)
saveRDS(D, "features/loco_predictions.rds")

cat("\n================ POOLED ACROSS FOLDS ================\n")
cat(sprintf("samples evaluated: %d (%d cohorts)\n", nrow(D), length(folds)))
for (m in c("null","en","rf")) {
  p <- D[[m]]
  cat(sprintf("%-5s  MAE %6.2f  RMSE %6.2f  Pearson %6.3f  Spearman %6.3f  R2 %6.3f\n",
      m, mean(abs(p-D$age)), sqrt(mean((p-D$age)^2)),
      cor(p,D$age), cor(p,D$age,method="spearman"),
      1 - sum((p-D$age)^2)/sum((D$age-mean(D$age))^2)))
}
cat("\nimprovement over null: elastic net %.2f yr (%.1f%%) | RF %.2f yr (%.1f%%)\n")
en <- mean(abs(D$en-D$age)); rf <- mean(abs(D$rf-D$age)); nu <- mean(abs(D$null-D$age))
cat(sprintf("  elastic net: %.2f yr better (%.1f%%)\n  random forest: %.2f yr better (%.1f%%)\n",
    nu-en, 100*(nu-en)/nu, nu-rf, 100*(nu-rf)/nu))
cat("\nslope of observed~predicted (1.0 = unbiased, <1 = compression):\n")
cat(sprintf("  elastic net %.3f | random forest %.3f\n",
    coef(lm(age ~ en, D))[2], coef(lm(age ~ rf, D))[2]))
