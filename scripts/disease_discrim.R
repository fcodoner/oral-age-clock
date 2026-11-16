# --------------------------------------------------------------------------
# disease_discrim.R
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

auc <- function(s, d) { r <- rank(s); n1 <- sum(d==1); n0 <- sum(d==0)
  if (n1==0 || n0==0) return(NA_real_)
  (sum(r[d==1]) - n1*(n1+1)/2) / (n1*n0) }

cat("=== A. WITHIN-COHORT 5-fold CV (does a disease signal exist?) ===\n")
cat(sprintf("%-18s %5s %5s %28s %7s\n","cohort","nHlth","nDis","condition(s)","AUC"))
tab <- list()
for (c in unique(sm$cohort_id)) {
  i <- sm$cohort_id==c; d <- sm$dis[i]
  if (sum(d==1) < 8 || sum(d==0) < 8) next
  Xi <- X[i,,drop=FALSE]; k <- sample(rep(1:5, length.out=sum(i)))
  p <- numeric(sum(i))
  for (j in 1:5) {
    if (length(unique(d[k!=j])) < 2) next
    cv <- cv.glmnet(Xi[k!=j,], d[k!=j], family="binomial", alpha=0.5, nfolds=5)
    p[k==j] <- as.numeric(predict(cv, Xi[k==j,], s="lambda.min"))
  }
  a <- auc(p, d)
  cond <- paste(names(sort(table(sm$label[i & sm$dis==1]), decreasing=TRUE))[1:2], collapse="/")
  cat(sprintf("%-18s %5d %5d %28s %7.3f\n", c, sum(d==0), sum(d==1),
      substr(gsub("NA","",cond),1,28), a))
  tab[[c]] <- a
}
cat(sprintf("\nmedian within-cohort AUC: %.3f (n=%d cohorts)\n",
    median(unlist(tab), na.rm=TRUE), length(tab)))

cat("\n=== B. CROSS-COHORT (train on other cohorts, test on held-out) ===\n")
grp <- list(
  OSCC          = c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
  periodontitis = c("PORTUGAL_74247","SPAIN_74796"))
for (g in names(grp)) {
  cs <- intersect(grp[[g]], unique(sm$cohort_id))
  keep <- sm$cohort_id %in% cs &
          (sm$dis==0 | grepl(g, sm$label, ignore.case=TRUE) |
           (g=="OSCC" & grepl("OSCC|carcinoma", sm$label, ignore.case=TRUE)))
  cat(sprintf("\n-- %s: %d cohorts, %d healthy, %d cases --\n", g, length(cs),
      sum(keep & sm$dis==0), sum(keep & sm$dis==1)))
  for (f in cs) {
    te <- keep & sm$cohort_id==f; tr <- keep & sm$cohort_id!=f
    if (sum(sm$dis[te]==1) < 5 || sum(sm$dis[te]==0) < 5) {
      cat(sprintf("   %-18s skipped (too few on one side)\n", f)); next }
    cv <- cv.glmnet(X[tr,], sm$dis[tr], family="binomial", alpha=0.5, nfolds=5)
    p  <- as.numeric(predict(cv, X[te,], s="lambda.min"))
    cat(sprintf("   %-18s n=%3d (%d case) AUC %.3f\n", f, sum(te),
        sum(sm$dis[te]==1), auc(p, sm$dis[te])))
  }
}
