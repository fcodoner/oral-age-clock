# --------------------------------------------------------------------------
# gradient_tests.R
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
M  <- M[, colMeans(M>0) > 0.10, drop=FALSE]
L  <- log(M + 0.5); X <- L - rowMeans(L)
sm$dis <- ifelse(sm$class=="train", 0L, 1L)
auc <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }

cat("=== A. DEDICATED PRECANCEROUS MODEL (is 0.59 the ceiling?) ===\n")
pc <- sm$cohort_id=="LEE17_CN" & grepl("precancer", sm$label, ignore.case=TRUE)
hl <- sm$cohort_id=="LEE17_CN" & sm$dis==0
cat(sprintf("Lee17: healthy n=%d (mean age %.1f) | precancerous n=%d (mean age %.1f)\n",
    sum(hl), mean(sm$age[hl]), sum(pc), mean(sm$age[pc])))
i <- which(hl | pc); d <- as.integer(pc[i])
k <- sample(rep(1:5, length.out=length(i))); p <- numeric(length(i))
for (j in 1:5) {
  cv <- cv.glmnet(X[i[k!=j],], d[k!=j], family="binomial", alpha=0.5, nfolds=5)
  p[k==j] <- as.numeric(predict(cv, X[i[k==j],], s="lambda.min"))
}
cat(sprintf("  within-Lee17 5-fold CV, healthy vs precancerous AUC: %.3f  (optimistic ceiling)\n", auc(p,d)))
ovh <- sm$cohort_id=="TAIWAN_39064" & grepl("verrucous", sm$label, ignore.case=TRUE)
th  <- sm$cohort_id=="TAIWAN_39064" & sm$dis==0
if (sum(ovh)>=5) {
  cv <- cv.glmnet(X[i,], d, family="binomial", alpha=0.5, nfolds=5)
  jj <- which(ovh | th)
  cat(sprintf("  Lee17-trained -> Taiwan OVH vs healthy AUC: %.3f (n=%d)\n",
      auc(as.numeric(predict(cv, X[jj,], s="lambda.min")), as.integer(ovh[jj])), length(jj)))
}

cat("\n=== B. PERIODONTITIS MODEL -> RELVAS ORAL-HEALTH GRADIENT ===\n")
tr <- which(sm$cohort_id %in% c("PORTUGAL_74247","SPAIN_74796") &
            (sm$dis==0 | grepl("periodontitis", sm$label, ignore.case=TRUE)))
cv <- cv.glmnet(X[tr,], sm$dis[tr], family="binomial", alpha=0.5, nfolds=5)
rv <- which(sm$cohort_id=="RELVAS21_EU")
s  <- as.numeric(predict(cv, X[rv,], s="lambda.min"))
gr <- ifelse(sm$dis[rv]==0, "0_healthy", sub(".*Grade_", "grade_", sm$label[rv]))
cat(sprintf("trained on %d Iberian samples; Relvas never seen\n\n", length(tr)))
for (g in sort(unique(gr)))
  cat(sprintf("  %-12s n=%3d  mean score %+6.3f  median %+6.3f\n",
      g, sum(gr==g), mean(s[gr==g]), median(s[gr==g])))
o <- as.integer(factor(gr))
ct <- suppressWarnings(cor.test(s, o, method="spearman"))
cat(sprintf("\n  ordered trend across grades: rho %.3f, p %.3g\n", ct$estimate, ct$p.value))
cat(sprintf("  healthy vs any-disease AUC : %.3f\n", auc(s, as.integer(gr!="0_healthy"))))
cat(sprintf("  healthy vs grade_3     AUC : %.3f (n=%d)\n",
    auc(s[gr %in% c("0_healthy","grade_3")], as.integer(gr[gr %in% c("0_healthy","grade_3")]=="grade_3")),
    sum(gr %in% c("0_healthy","grade_3"))))
cat(sprintf("  age of Relvas groups: %s\n",
    paste(sprintf("%s %.0fy", sort(unique(gr)), tapply(sm$age[rv], gr, mean)), collapse=" | ")))
