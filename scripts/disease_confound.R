# --------------------------------------------------------------------------
# disease_confound.R
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
  if(n1==0||n0==0) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }

sets <- list(OSCC          = c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
             periodontitis = c("PORTUGAL_74247","SPAIN_74796"))
for (g in names(sets)) {
  cs <- sets[[g]]
  keep <- sm$cohort_id %in% cs & (sm$dis==0 |
          grepl(if (g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  cat(sprintf("\n########## %s ##########\n", g))
  cat(sprintf("%-18s %10s %10s %8s\n","cohort","age ctrl","age case","diff"))
  for (f in cs) {
    i <- keep & sm$cohort_id==f
    a0 <- sm$age[i & sm$dis==0]; a1 <- sm$age[i & sm$dis==1]
    if (!length(a0) || !length(a1)) next
    cat(sprintf("%-18s %4.1f (n=%2d) %4.1f (n=%3d) %+8.1f\n", f,
        mean(a0), length(a0), mean(a1), length(a1), mean(a1)-mean(a0)))
  }
  cat(sprintf("\n%-18s %8s %8s %8s %8s\n","held-out","ageOnly","microb","micr+age","ageMatch"))
  for (f in cs) {
    te <- keep & sm$cohort_id==f; tr <- keep & sm$cohort_id!=f
    if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
    d  <- sm$dis[te]
    a1 <- auc(sm$age[te], d)
    cv <- cv.glmnet(X[tr,], sm$dis[tr], family="binomial", alpha=0.5, nfolds=5)
    a2 <- auc(as.numeric(predict(cv, X[te,], s="lambda.min")), d)
    XA <- cbind(X, age=sm$age)
    cv3<- cv.glmnet(XA[tr,], sm$dis[tr], family="binomial", alpha=0.5, nfolds=5)
    a3 <- auc(as.numeric(predict(cv3, XA[te,], s="lambda.min")), d)
    # age-matched held-out subset: restrict to the overlapping age range
    lo <- max(min(sm$age[te & d==0]), min(sm$age[te & d==1]))
    hi <- min(max(sm$age[te & d==0]), max(sm$age[te & d==1]))
    mm <- te & sm$age >= lo & sm$age <= hi
    a4 <- if (sum(sm$dis[mm]==1)>=5 && sum(sm$dis[mm]==0)>=5)
            auc(as.numeric(predict(cv, X[mm,], s="lambda.min")), sm$dis[mm]) else NA_real_
    cat(sprintf("%-18s %8.3f %8.3f %8.3f %8.3f  [matched n=%d, %.0f-%.0f y]\n",
        f, a1, a2, a3, a4, sum(mm), lo, hi))
  }
}
