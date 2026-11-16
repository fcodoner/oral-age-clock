# --------------------------------------------------------------------------
# gradients.R
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
L  <- log(M + 0.5); X <- L - rowMeans(L)

CONT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
          "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
          "Burkholderia","Cutibacterium","Escherichia")
cc <- intersect(CONT, colnames(M))
R  <- M/rowSums(M)
cat("=== suspected reagent/environmental contaminants ===\n")
cat(sprintf("%-20s %9s %10s %12s\n","genus","mean %","prevalence","cor w/ depth"))
for (g in cc) cat(sprintf("%-20s %8.3f%% %9.1f%% %12.3f\n", g, 100*mean(R[,g]),
    100*mean(R[,g]>0), cor(R[,g], rowSums(M), method="spearman")))
cat(sprintf("\ntotal abundance of these taxa: %.2f%% of reads\n", 100*mean(rowSums(R[,cc,drop=FALSE]))))

X <- X[, !colnames(X) %in% cc, drop=FALSE]
cat(sprintf("\ngradient analysis on %d genera (contaminants removed)\n", ncol(X)))
pc <- prcomp(X, center=TRUE, scale.=FALSE)
ve <- 100*pc$sdev^2/sum(pc$sdev^2)
cat(sprintf("variance explained: PC1 %.1f%%  PC2 %.1f%%  PC3 %.1f%%  PC4 %.1f%%  PC5 %.1f%%\n",
    ve[1],ve[2],ve[3],ve[4],ve[5]))

bc <- function(z) { n<-length(z); s<-sum((z-mean(z))^3)/n/sd(z)^3
  k<-sum((z-mean(z))^4)/n/sd(z)^4 - 3; (s^2+1)/(k + 3*(n-1)^2/((n-2)*(n-3))) }
cat("\n--- continuous or discrete? (bimodality coefficient; >0.555 suggests two modes) ---\n")
for (i in 1:4) cat(sprintf("PC%d: BC %.3f  %s\n", i, bc(pc$x[,i]),
    ifelse(bc(pc$x[,i])>0.555, "possible two modes", "unimodal / continuous")))

cat("\n--- what each gradient tracks ---\n")
r2 <- function(y, f) summary(lm(y ~ factor(f)))$r.squared
cat(sprintf("%5s %9s %9s %9s %9s %9s\n","PC","age r","cohort R2","region R2","route R2","depth r"))
for (i in 1:4) cat(sprintf("%5d %9.3f %9.3f %9.3f %9.3f %9.3f\n", i,
    cor(pc$x[,i], sm$age), r2(pc$x[,i], sm$cohort_id), r2(pc$x[,i], sm$region),
    r2(pc$x[,i], sm$route), cor(pc$x[,i], rowSums(M), method="spearman")))

cat("\n--- top loadings ---\n")
for (i in 1:3) {
  l <- sort(pc$rotation[,i])
  cat(sprintf("\nPC%d (%.1f%%)\n  negative: %s\n  positive: %s\n", i, ve[i],
      paste(names(head(l,6)), collapse=", "),
      paste(rev(names(tail(l,6))), collapse=", ")))
}

source("scripts/predict_sample.R")
oh <- sapply(1:nrow(M), function(k) score_sample(M[k,], verbose=FALSE)$oral_health_prob)
cat("\n--- do gradients relate to the oral-health score? ---\n")
for (i in 1:4) cat(sprintf("PC%d vs oral-health score: Spearman %+.3f\n", i,
    cor(pc$x[,i], oh, method="spearman")))
saveRDS(list(pca=pc, var=ve, oh=oh, samples=rownames(M), dropped=cc), "features/gradients.rds")
cat("\nsaved: features/gradients.rds\n")
