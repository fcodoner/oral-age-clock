# --------------------------------------------------------------------------
# age_local.R
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

cat("=== A. TRAINED AND TESTED WITHIN ONE LABORATORY (5-fold CV) ===\n")
cat(sprintf("%-18s %5s %7s %8s %8s %7s %7s\n","cohort","n","ageSD","nullMAE","MAE","r","R2"))
for (f in big) {
  i <- which(co==f); k <- sample(rep(1:5, length.out=length(i))); p <- numeric(length(i))
  for (j in 1:5) {
    cv <- cv.glmnet(X[i[k!=j],], y[i[k!=j]], alpha=0.5, nfolds=5, standardize=TRUE)
    p[k==j] <- as.numeric(predict(cv, X[i[k==j],], s="lambda.min"))
  }
  cat(sprintf("%-18s %5d %7.1f %8.2f %8.2f %7.3f %7.3f\n", f, length(i), sd(y[i]),
      mean(abs(y[i]-mean(y[i]))), mean(abs(p-y[i])), cor(p,y[i]),
      1-sum((p-y[i])^2)/sum((y[i]-mean(y[i]))^2)))
}
cat("\n(cross-cohort reference: MAE 12.18, r 0.449, R2 0.177)\n")

cat("\n=== B. SHIPPED MODEL + LOCAL CALIBRATION ===\n")
cat("Train on all other cohorts, then map predictions onto the target cohort's own\n")
cat("age distribution by quantile matching -- what a clinic could do with a reference panel.\n")
cat(sprintf("\n%-18s %5s %8s %10s %12s %7s\n","held-out","n","nullMAE","raw MAE","calibMAE","r"))
for (f in big) {
  te <- which(co==f); tr <- which(co!=f)
  cv <- cv.glmnet(X[tr,], y[tr], alpha=0.5, nfolds=5, standardize=TRUE)
  p  <- as.numeric(predict(cv, X[te,], s="lambda.min"))
  cal <- quantile(y[te], probs=rank(p)/(length(p)+1))   # quantile mapping
  cat(sprintf("%-18s %5d %8.2f %10.2f %12.2f %7.3f\n", f, length(te),
      mean(abs(y[te]-mean(y[tr]))), mean(abs(p-y[te])), mean(abs(cal-y[te])), cor(p,y[te])))
}
