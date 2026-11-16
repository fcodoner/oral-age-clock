# --------------------------------------------------------------------------
# disease_ageadj.R
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

sets <- list(OSCC          = c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
             periodontitis = c("PORTUGAL_74247","SPAIN_74796"))
NREP <- 200; TOL <- 5
for (g in names(sets)) {
  cs <- sets[[g]]
  keep <- sm$cohort_id %in% cs & (sm$dis==0 |
          grepl(if (g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  cat(sprintf("\n########## %s ##########\n", g))
  cat(sprintf("%-18s %8s %9s %10s %9s %9s\n",
      "held-out","AUCraw","OR/SD","p(score)","matchAUC","matchN"))
  for (f in cs) {
    te <- which(keep & sm$cohort_id==f); tr <- which(keep & sm$cohort_id!=f)
    if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
    cv <- cv.glmnet(X[tr,], sm$dis[tr], family="binomial", alpha=0.5, nfolds=5)
    s  <- as.numeric(predict(cv, X[te,], s="lambda.min"))
    d  <- sm$dis[te]; a <- sm$age[te]
    sz <- (s - mean(s))/sd(s)
    fit <- glm(d ~ a + sz, family=binomial)
    cf  <- summary(fit)$coefficients
    # 1:1 age-matched case-control resampling, |age difference| <= TOL
    ca <- which(d==1); ct <- which(d==0); accs <- c()
    for (rep in 1:NREP) {
      pool <- ct; pick_c <- c(); pick_t <- c()
      for (i in sample(ca)) {
        cand <- pool[abs(a[pool] - a[i]) <= TOL]
        if (!length(cand)) next
        j <- if (length(cand)==1) cand else sample(cand, 1)
        pick_c <- c(pick_c, i); pick_t <- c(pick_t, j); pool <- setdiff(pool, j)
      }
      if (length(pick_c) >= 5) accs <- c(accs, auc(s[c(pick_c,pick_t)], d[c(pick_c,pick_t)]))
    }
    cat(sprintf("%-18s %8.3f %9.2f %10.3g %9.3f %9d\n", f, auc(s,d),
        exp(cf["sz","Estimate"]), cf["sz","Pr(>|z|)"],
        mean(accs, na.rm=TRUE), if (length(accs)) 2*length(pick_c) else 0))
  }
}
cat("\nOR/SD = odds ratio per SD of the microbiome score, adjusted for age\n")
cat("p(score) = does the microbiome add information beyond age alone?\n")
cat(sprintf("matchAUC = mean AUC over %d 1:1 age-matched resamples (+/-%d years)\n", NREP, TOL))
