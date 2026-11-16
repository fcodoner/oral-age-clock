# --------------------------------------------------------------------------
# age_calibration.R
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
CONT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
          "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
          "Burkholderia","Cutibacterium","Escherichia")
K <- setdiff(colnames(M)[colMeans(M>0)>0.10], CONT)
L <- log(M[,K,drop=FALSE]+0.5); X <- L - rowMeans(L)
y <- sm$age; co <- sm$cohort_id
nn <- table(co); ag <- tapply(y, co, function(v) diff(range(v)))
big <- names(nn)[nn >= 200 & ag[names(nn)] >= 40]

cat("=== A. HOW GENTLE SHOULD LOCAL CALIBRATION BE? ===\n")
cat("A local panel of 50 age-labelled samples is used only to CALIBRATE a shipped model.\n")
cat(sprintf("\n%-18s %8s %9s %9s %9s %9s\n",
    "held-out","raw MAE","shift","shift+scale","quantile","r"))
for (f in big) {
  te <- which(co==f); tr <- which(co!=f)
  cv <- cv.glmnet(X[tr,], y[tr], alpha=0.5, nfolds=5, standardize=TRUE)
  p  <- as.numeric(predict(cv, X[te,], s="lambda.min"))
  res <- replicate(200, {
    pan <- sample(length(te), 50); rest <- setdiff(seq_along(te), pan)
    b   <- mean(y[te][pan] - p[pan])                       # intercept only
    lm1 <- lm(y[te][pan] ~ p[pan])                         # intercept + slope
    s2  <- coef(lm1)[1] + coef(lm1)[2]*p[rest]
    qm  <- quantile(y[te][pan], probs=rank(p[rest])/(length(rest)+1))
    c(raw=mean(abs(p[rest]-y[te][rest])), shift=mean(abs(p[rest]+b-y[te][rest])),
      scale=mean(abs(s2-y[te][rest])),   quant=mean(abs(qm-y[te][rest])))
  })
  m <- rowMeans(res)
  cat(sprintf("%-18s %8.2f %9.2f %9.2f %9.2f %9.3f\n", f,
      m["raw"], m["shift"], m["scale"], m["quant"], cor(p, y[te])))
}

cat("\n=== B. HOW MANY LOCAL TRAINING SAMPLES BEFORE LOCAL BEATS SHIPPED? ===\n")
for (f in big) {
  i <- which(co==f); tr <- which(co!=f)
  cv <- cv.glmnet(X[tr,], y[tr], alpha=0.5, nfolds=5, standardize=TRUE)
  cat(sprintf("\n-- %s (n=%d, shipped MAE %.2f, null %.2f) --\n", f, length(i),
      mean(abs(as.numeric(predict(cv, X[i,], s="lambda.min")) - y[i])),
      mean(abs(y[i]-mean(y[i])))))
  cat(sprintf("%8s %9s %8s\n","n_local","MAE","sd"))
  for (n in c(50,100,150,200,300)) {
    if (n > length(i)-60) next
    v <- replicate(20, {
      k <- sample(length(i), n); ho <- setdiff(seq_along(i), k)
      m <- cv.glmnet(X[i[k],], y[i[k]], alpha=0.5, nfolds=5, standardize=TRUE)
      mean(abs(as.numeric(predict(m, X[i[ho],], s="lambda.min")) - y[i[ho]]))
    })
    cat(sprintf("%8d %9.2f %8.2f\n", n, mean(v), sd(v)))
  }
}
