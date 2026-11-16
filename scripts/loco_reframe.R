# --------------------------------------------------------------------------
# loco_reframe.R
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
L <- log(M + 0.5); X <- L - rowMeans(L)
y <- sm$age; co <- sm$cohort_id; rg <- sm$region
ag <- tapply(y, co, function(v) diff(range(v))); nn <- table(co)
folds <- names(ag)[ag >= 40 & nn[names(ag)] >= 40]

cat("=== 1. RANDOM 5-FOLD CV (the optimistic, publication-style estimate) ===\n")
k <- sample(rep(1:5, length.out=nrow(X)))
P <- numeric(nrow(X))
for (i in 1:5) {
  cv <- cv.glmnet(X[k!=i,], y[k!=i], alpha=0.5, nfolds=5, standardize=TRUE)
  P[k==i] <- as.numeric(predict(cv, X[k==i,], s="lambda.min"))
}
cat(sprintf("random CV : MAE %.2f | r %.3f | R2 %.3f\n", mean(abs(P-y)), cor(P,y),
    1-sum((P-y)^2)/sum((y-mean(y))^2)))
cat(sprintf("LOCO      : MAE 12.17 | r 0.449 | R2 0.177   (from previous run)\n"))

cat("\n=== 2. LOCO, V3-V4 COHORTS ONLY (longest amplicon, one region) ===\n")
i34 <- rg == "V3-V4"
f34 <- intersect(folds, unique(co[i34]))
cat("folds:", length(f34), "| training pool:", sum(i34), "samples\n")
Pv <- c(); Av <- c(); Cv <- c()
for (f in f34) {
  te <- i34 & co==f; tr <- i34 & co!=f
  cv <- cv.glmnet(X[tr,], y[tr], alpha=0.5, nfolds=5, standardize=TRUE)
  Pv <- c(Pv, as.numeric(predict(cv, X[te,], s="lambda.min")))
  Av <- c(Av, y[te]); Cv <- c(Cv, rep(f, sum(te)))
}
wr <- sapply(unique(Cv), function(f){i<-Cv==f; cor(Pv[i],Av[i])})
cat(sprintf("V3-V4 LOCO: MAE %.2f | r %.3f | medWithinR %.3f\n",
    mean(abs(Pv-Av)), cor(Pv,Av), median(wr)))
print(round(sort(wr), 3))

cat("\n=== 3. LOCO AGE-BAND CLASSIFICATION (3 bands) ===\n")
br <- quantile(y, c(0,1/3,2/3,1)); br[1] <- br[1]-1
g  <- cut(y, br, labels=c("young","mid","older"))
cat("bands:", paste(sprintf("%s (%.0f-%.0f, n=%d)", levels(g),
    br[-4]+c(1,0,0), br[-1], table(g)), collapse=" | "), "\n")
TR <- c(); PR <- c()
for (f in folds) {
  te <- co==f; tr <- !te
  cv <- cv.glmnet(X[tr,], g[tr], family="multinomial", alpha=0.5,
                  nfolds=5, standardize=TRUE)
  pc <- predict(cv, X[te,], s="lambda.min", type="class")
  TR <- c(TR, as.character(g[te])); PR <- c(PR, as.character(pc))
}
lv <- c("young","mid","older")
acc <- mean(TR==PR)
adj <- mean(abs(match(TR,lv)-match(PR,lv)) <= 1)
cat(sprintf("exact accuracy %.3f (chance 0.333) | within-one-band %.3f\n", acc, adj))
print(table(true=factor(TR,lv), pred=factor(PR,lv)))
