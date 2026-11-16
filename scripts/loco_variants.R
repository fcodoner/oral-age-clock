# --------------------------------------------------------------------------
# loco_variants.R
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
M  <- M[sm$class=="train", ]; sm <- sm[sm$class=="train", ]
M  <- M[, colMeans(M>0) > 0.10, drop=FALSE]
L <- log(M + 0.5); X0 <- L - rowMeans(L)
y <- sm$age; co <- sm$cohort_id
ag <- tapply(y, co, function(v) diff(range(v))); nn <- table(co)
folds <- names(ag)[ag >= 40 & nn[names(ag)] >= 40]

zc <- function(X, co) {            # z-score each feature within its own cohort
  out <- X
  for (c in unique(co)) { i <- co==c
    mu <- colMeans(X[i,,drop=FALSE]); sd <- apply(X[i,,drop=FALSE],2,sd)
    sd[sd==0] <- 1
    out[i,] <- sweep(sweep(X[i,,drop=FALSE],2,mu,"-"),2,sd,"/") }
  out }

reg <- model.matrix(~ region - 1, data=sm)
variants <- list(
  baseline      = list(X=X0,                  drop=NULL),
  cohort_zscore = list(X=zc(X0, co),          drop=NULL),
  plus_region   = list(X=cbind(X0, reg),      drop=NULL),
  no_AGP        = list(X=X0,                  drop="AGP_ORAL"))

cat(sprintf("%-14s %7s %7s %7s %7s %9s\n",
    "variant","MAE","RMSE","r","R2","medWithinR"))
for (v in names(variants)) {
  X <- variants[[v]]$X; dr <- variants[[v]]$drop
  P <- c(); A <- c(); CO <- c()
  for (f in folds) {
    te <- co==f; tr <- !te & !(co %in% dr)
    if (!any(tr)) next
    cv <- cv.glmnet(X[tr,], y[tr], alpha=0.5, nfolds=5, standardize=TRUE)
    P <- c(P, as.numeric(predict(cv, X[te,], s="lambda.min")))
    A <- c(A, y[te]); CO <- c(CO, rep(f, sum(te)))
  }
  wr <- sapply(unique(CO), function(f) { i <- CO==f
    if (sd(P[i])==0) NA_real_ else suppressWarnings(cor(P[i], A[i])) })
  cat(sprintf("%-14s %7.2f %7.2f %7.3f %7.3f %9.3f\n", v,
      mean(abs(P-A)), sqrt(mean((P-A)^2)), cor(P,A),
      1-sum((P-A)^2)/sum((A-mean(A))^2), median(wr, na.rm=TRUE)))
}
