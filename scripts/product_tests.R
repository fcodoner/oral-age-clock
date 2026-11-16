# --------------------------------------------------------------------------
# product_tests.R
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

cat("=== A. PRECANCEROUS DOSE-RESPONSE (model never saw Lee17) ===\n")
trc <- c("KUMPITSCH20_EU","TAIWAN_39064","WOLF17_EU")
tr  <- which(sm$cohort_id %in% trc &
             (sm$dis==0 | grepl("OSCC|carcinoma", sm$label, ignore.case=TRUE)))
cv  <- cv.glmnet(X[tr,], sm$dis[tr], family="binomial", alpha=0.5, nfolds=5)
le  <- which(sm$cohort_id=="LEE17_CN")
s   <- as.numeric(predict(cv, X[le,], s="lambda.min"))
grp <- ifelse(sm$dis[le]==0, "1_healthy",
       ifelse(grepl("precancer", sm$label[le], ignore.case=TRUE), "2_precancerous", "3_OSCC"))
cat(sprintf("trained on %d samples from %s\n\n", length(tr), paste(trc, collapse=", ")))
for (g in sort(unique(grp)))
  cat(sprintf("  %-16s n=%3d  mean score %+6.3f  median %+6.3f\n",
      g, sum(grp==g), mean(s[grp==g]), median(s[grp==g])))
cat(sprintf("\n  healthy vs precancerous AUC : %.3f\n", auc(s[grp!="3_OSCC"], as.integer(grp[grp!="3_OSCC"]=="2_precancerous"))))
cat(sprintf("  healthy vs OSCC          AUC : %.3f\n", auc(s[grp!="2_precancerous"], as.integer(grp[grp!="2_precancerous"]=="3_OSCC"))))
cat(sprintf("  precancerous vs OSCC     AUC : %.3f\n", auc(s[grp!="1_healthy"], as.integer(grp[grp!="1_healthy"]=="3_OSCC"))))
cat(sprintf("  Spearman trend across 3 ordered groups: rho %.3f, p %.3g\n",
    cor(s, as.integer(factor(grp)), method="spearman"),
    cor.test(s, as.integer(factor(grp)), method="spearman")$p.value))

cat("\n=== B. REDUCED PANEL (how few genera suffice?) ===\n")
sets <- list(OSCC          = c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
             periodontitis = c("PORTUGAL_74247","SPAIN_74796"))
for (g in names(sets)) {
  cs <- sets[[g]]
  keep <- sm$cohort_id %in% cs & (sm$dis==0 |
          grepl(if (g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  cat(sprintf("\n-- %s --\n", g))
  cat(sprintf("%-18s %8s %8s %8s %8s %8s\n","held-out","top10","top20","top30","top50","all"))
  for (f in cs) {
    te <- which(keep & sm$cohort_id==f); trr <- which(keep & sm$cohort_id!=f)
    if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
    cvf <- cv.glmnet(X[trr,], sm$dis[trr], family="binomial", alpha=0.5, nfolds=5)
    b   <- as.matrix(coef(cvf, s="lambda.min"))[-1,1]
    ord <- names(sort(abs(b), decreasing=TRUE))
    out <- c()
    for (k in c(10,20,30,50)) {
      sel <- ord[1:min(k,length(ord))]
      c2  <- cv.glmnet(X[trr,sel,drop=FALSE], sm$dis[trr], family="binomial", alpha=0.5, nfolds=5)
      out <- c(out, auc(as.numeric(predict(c2, X[te,sel,drop=FALSE], s="lambda.min")), sm$dis[te]))
    }
    out <- c(out, auc(as.numeric(predict(cvf, X[te,], s="lambda.min")), sm$dis[te]))
    cat(sprintf("%-18s %8.3f %8.3f %8.3f %8.3f %8.3f\n", f, out[1],out[2],out[3],out[4],out[5]))
  }
}
cat("\n=== C. GENERA SELECTED MOST OFTEN ACROSS FOLDS ===\n")
for (g in names(sets)) {
  cs <- sets[[g]]
  keep <- sm$cohort_id %in% cs & (sm$dis==0 |
          grepl(if (g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  acc <- list()
  for (f in cs) {
    trr <- which(keep & sm$cohort_id!=f)
    cvf <- cv.glmnet(X[trr,], sm$dis[trr], family="binomial", alpha=0.5, nfolds=5)
    b <- as.matrix(coef(cvf, s="lambda.min"))[-1,1]
    acc[[f]] <- b[b!=0]
  }
  v <- unlist(acc); nm <- sub("^.*?\\.","",names(v))
  cnt <- tapply(v, nm, function(z) c(n=length(z), mean=mean(z)))
  dd <- data.frame(genus=names(cnt),
                   folds=sapply(cnt,`[`,1), meanCoef=round(sapply(cnt,`[`,2),3))
  dd <- dd[order(-dd$folds, -abs(dd$meanCoef)), ]
  cat(sprintf("\n-- %s (sign: + = higher in disease) --\n", g))
  print(head(dd, 15), row.names=FALSE)
}
