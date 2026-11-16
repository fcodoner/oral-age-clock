# --------------------------------------------------------------------------
# depth_check.R
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
sm$dis <- ifelse(sm$class=="train", 0L, 1L); sm$depth <- rowSums(M)
CONT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
          "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
          "Burkholderia","Cutibacterium","Escherichia")
auc <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }
mk <- function(drop) { K <- colnames(M)[colMeans(M>0)>0.10]; K <- setdiff(K, drop)
  L <- log(M[,K,drop=FALSE]+0.5); L - rowMeans(L) }
X1 <- mk(character(0)); X2 <- mk(CONT)

sets <- list(OSCC=c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
             periodontitis=c("PORTUGAL_74247","SPAIN_74796"))
for (g in names(sets)) {
  cs <- sets[[g]]
  keep <- sm$cohort_id %in% cs & (sm$dis==0 |
          grepl(if (g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  cat(sprintf("\n########## %s ##########\n", g))
  cat(sprintf("%-18s %10s %10s %9s %9s %9s\n",
      "held-out","depth ctrl","depth case","AUCdepth","AUCwith","AUCnoCont"))
  for (f in cs) {
    te <- which(keep & sm$cohort_id==f); tr <- which(keep & sm$cohort_id!=f)
    if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
    d <- sm$dis[te]
    a0 <- auc(sm$depth[te], d)
    c1 <- cv.glmnet(X1[tr,], sm$dis[tr], family="binomial", alpha=0.5, nfolds=5)
    c2 <- cv.glmnet(X2[tr,], sm$dis[tr], family="binomial", alpha=0.5, nfolds=5)
    cat(sprintf("%-18s %10.0f %10.0f %9.3f %9.3f %9.3f\n", f,
        median(sm$depth[te][d==0]), median(sm$depth[te][d==1]), a0,
        auc(as.numeric(predict(c1, X1[te,], s="lambda.min")), d),
        auc(as.numeric(predict(c2, X2[te,], s="lambda.min")), d)))
  }
}

cat("\n########## can 2 axes replace 113 genera? ##########\n")
gr <- readRDS("features/gradients.rds")
tr <- which(sm$class=="train")
cat(sprintf("healthy samples: PC2 vs oral-health score rho %.3f | PC3 vs age rho %.3f\n",
    cor(gr$pca$x[,2], gr$oh, method="spearman"),
    cor(gr$pca$x[,3], sm$age[tr], method="spearman")))
P <- predict(gr$pca, mk(gr$dropped))
for (g in names(sets)) {
  cs <- sets[[g]]
  keep <- sm$cohort_id %in% cs & (sm$dis==0 |
          grepl(if (g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  cat(sprintf("\n-- %s: AUC from PC2+PC3 alone --\n", g))
  for (f in cs) {
    te <- which(keep & sm$cohort_id==f); tr2 <- which(keep & sm$cohort_id!=f)
    if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
    fit <- glm(sm$dis[tr2] ~ P[tr2,2] + P[tr2,3], family=binomial)
    p <- cbind(1, P[te,2], P[te,3]) %*% coef(fit)
    cat(sprintf("   %-18s AUC %.3f\n", f, auc(as.numeric(p), sm$dis[te])))
  }
}
