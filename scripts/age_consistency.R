# --------------------------------------------------------------------------
# age_consistency.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
R  <- readRDS("features/genus_relab.rds"); rownames(R) <- unname(rownames(R))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(R), sm$run_accession), ]
R  <- R[sm$class=="train", ]; sm <- sm[sm$class=="train", ]
R  <- R[, colMeans(R>0) > 0.10, drop=FALSE]
cat(sprintf("training samples %d | genera %d\n", nrow(R), ncol(R)))

ag <- tapply(sm$age, sm$cohort_id, function(x) diff(range(x)))
nn <- table(sm$cohort_id)
use <- names(ag)[ag >= 40 & nn[names(ag)] >= 40]
cat("informative cohorts (span>=40y, n>=40):", length(use), "\n")
cat(paste(sprintf("%s(n=%d,span=%.0f)", use, nn[use], ag[use]), collapse=", "), "\n\n")

rho <- sapply(use, function(c) {
  i <- sm$cohort_id == c
  apply(R[i, , drop=FALSE], 2, function(y)
    if (sd(y)==0) NA_real_ else suppressWarnings(cor(y, sm$age[i], method="spearman")))
})
ok <- rowSums(!is.na(rho)) == length(use)
rho <- rho[ok, , drop=FALSE]
cat("genera testable in all cohorts:", nrow(rho), "\n\n")

pos <- rowSums(rho > 0); k <- ncol(rho)
agree <- pmax(pos, k - pos)
p <- sapply(agree, function(a) binom.test(a, k, 0.5)$p.value)
q <- p.adjust(p, "BH")
res <- data.frame(genus=rownames(rho), mean_rho=round(rowMeans(rho),3),
                  agree=paste0(agree,"/",k), q=signif(q,3), stringsAsFactors=FALSE)
res <- res[order(-abs(res$mean_rho)), ]

cat("--- sign agreement across cohorts ---\n")
cat(sprintf("genera agreeing in >=%d/%d cohorts: %d (q<0.05: %d)\n",
    ceiling(0.8*k), k, sum(agree >= ceiling(0.8*k)), sum(q < 0.05)))
cat("\n--- top 20 by mean Spearman rho with age ---\n")
print(head(res, 20), row.names=FALSE)
cat("\n--- strongest consistent (q<0.05, ranked by |mean rho|) ---\n")
print(head(res[res$q < 0.05, ], 20), row.names=FALSE)
saveRDS(rho, "features/per_cohort_age_rho.rds")
