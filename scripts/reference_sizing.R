# --------------------------------------------------------------------------
# reference_sizing.R
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
rv  <- which(sm$cohort_id=="RELVAS21_EU")
lab <- ifelse(sm$class[rv]=="train","0_healthy", sub(".*Grade_","grade_",sm$label[rv]))
H   <- sc[rv][lab=="0_healthy"]; D <- sc[rv][lab %in% c("grade_2","grade_3")]

cat(sprintf("%6s %12s %14s %18s\n","n_ref","q90 (mean)","q90 (sd)","%grade2-3 flagged"))
for (n in c(5,8,10,12,15,17)) {
  q <- c(); f <- c()
  for (i in 1:2000) {
    s <- sample(H, n, replace = n > length(H))
    t <- quantile(s, 0.90); q <- c(q, t); f <- c(f, mean(D >= t))
  }
  cat(sprintf("%6d %12.3f %14.3f %10.0f%% (sd %.0f%%)\n", n, mean(q), sd(q), 100*mean(f), 100*sd(f)))
}
cat("\n-- with replacement beyond n=17, so larger n shows the asymptote, not new data --\n")
for (n in c(25,40,60)) {
  q <- c(); f <- c()
  for (i in 1:2000) { s <- sample(H, n, replace=TRUE)
    t <- quantile(s, 0.90); q <- c(q, t); f <- c(f, mean(D >= t)) }
  cat(sprintf("%6d %12.3f %14.3f %10.0f%% (sd %.0f%%)\n", n, mean(q), sd(q), 100*mean(f), 100*sd(f)))
}
