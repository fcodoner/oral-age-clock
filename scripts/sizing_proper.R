# --------------------------------------------------------------------------
# sizing_proper.R
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

big <- c("TRIAD_ES","SLL2_ES","SOUTHKORE_13236")
cat("Each cohort is one laboratory. Draw n controls, set the 90th-percentile threshold,\n")
cat("then measure the false-positive rate on the HELD-OUT healthy from the same lab.\n")
cat("A well-sized panel gives 10% with little spread.\n\n")
cat(sprintf("%-18s %6s %12s %14s\n","cohort","n_ref","FPR mean","FPR sd"))
for (co in big) {
  H <- sc[which(sm$class=="train" & sm$cohort_id==co)]
  for (n in c(10,15,20,30,50,100)) {
    if (n >= length(H)-20) next
    f <- replicate(1000, { i <- sample(length(H), n)
      mean(H[-i] >= quantile(H[i], 0.90)) })
    cat(sprintf("%-18s %6d %11.1f%% %13.1f%%\n", co, n, 100*mean(f), 100*sd(f)))
  }
  cat("\n")
}
