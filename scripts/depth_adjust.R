# --------------------------------------------------------------------------
# depth_adjust.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages({library(glmnet); library(vegan)}); set.seed(42)
M  <- readRDS("features/genus_counts.rds"); rownames(M) <- unname(rownames(M))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(M), sm$run_accession), ]
sm$dis <- ifelse(sm$class=="train", 0L, 1L); sm$depth <- rowSums(M)
CONT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
          "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
          "Burkholderia","Cutibacterium","Escherichia")
auc <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }
K  <- setdiff(colnames(M)[colMeans(M>0)>0.10], CONT)
clrm <- function(C){ L <- log(C[,K,drop=FALSE]+0.5); L - rowMeans(L) }
X  <- clrm(M)
RAR <- 10000
Mr <- M[sm$depth >= RAR, , drop=FALSE]
Mr <- rrarefy(Mr, RAR); Xr <- clrm(Mr)
keepr <- rownames(M) %in% rownames(Mr)
cat(sprintf("rarefied to %d reads: %d of %d samples retained\n", RAR, nrow(Mr), nrow(M)))

sets <- list(OSCC=c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
             periodontitis=c("PORTUGAL_74247","SPAIN_74796"))
for (g in names(sets)) {
  cs <- sets[[g]]
  base <- sm$cohort_id %in% cs & (sm$dis==0 |
          grepl(if (g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  cat(sprintf("\n########## %s ##########\n", g))
  cat(sprintf("%-18s %10s %10s %9s %9s %8s\n",
      "held-out","|AUCdepth|","OR/SD adj","p(score)","AUCrarefied","n_rar"))
  for (f in cs) {
    te <- which(base & sm$cohort_id==f); tr <- which(base & sm$cohort_id!=f)
    if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
    d <- sm$dis[te]
    ad <- auc(sm$depth[te], d); ad <- max(ad, 1-ad)
    cv <- cv.glmnet(X[tr,], sm$dis[tr], family="binomial", alpha=0.5, nfolds=5)
    s  <- as.numeric(predict(cv, X[te,], s="lambda.min")); sz <- (s-mean(s))/sd(s)
    fit <- glm(d ~ log(sm$depth[te]) + sz, family=binomial)
    cf  <- summary(fit)$coefficients
    ter <- which(base & sm$cohort_id==f & keepr); trr <- which(base & sm$cohort_id!=f & keepr)
    ar <- NA_real_
    if (length(ter)>=10 && sum(sm$dis[ter]==1)>=5 && sum(sm$dis[ter]==0)>=5) {
      i1 <- match(rownames(M)[trr], rownames(Mr)); i2 <- match(rownames(M)[ter], rownames(Mr))
      c2 <- cv.glmnet(Xr[i1,], sm$dis[trr], family="binomial", alpha=0.5, nfolds=5)
      ar <- auc(as.numeric(predict(c2, Xr[i2,], s="lambda.min")), sm$dis[ter])
    }
    cat(sprintf("%-18s %10.3f %10.2f %9.3g %9.3f %8d\n", f, ad,
        exp(cf["sz","Estimate"]), cf["sz","Pr(>|z|)"], ar, length(ter)))
  }
}

cat("\n########## PC2 + PC3 as a compact, explainable model ##########\n")
tr <- which(sm$class=="train")
pc <- prcomp(X[tr,], center=TRUE, scale.=FALSE)
P  <- scale(X, center=pc$center, scale=FALSE) %*% pc$rotation
ve <- 100*pc$sdev^2/sum(pc$sdev^2)
cat(sprintf("PC1 %.1f%% (depth r %.2f) | PC2 %.1f%% | PC3 %.1f%% (age r %.2f)\n",
    ve[1], cor(P[tr,1], sm$depth[tr]), ve[2], ve[3], cor(P[tr,3], sm$age[tr])))
for (g in names(sets)) {
  cs <- sets[[g]]
  base <- sm$cohort_id %in% cs & (sm$dis==0 |
          grepl(if (g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  cat(sprintf("\n-- %s --\n%-18s %10s %10s\n", g, "held-out","PC2+PC3","PC1+PC2+PC3"))
  for (f in cs) {
    te <- which(base & sm$cohort_id==f); tr2 <- which(base & sm$cohort_id!=f)
    if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
    f1 <- glm(sm$dis[tr2] ~ P[tr2,2] + P[tr2,3], family=binomial)
    f2 <- glm(sm$dis[tr2] ~ P[tr2,1] + P[tr2,2] + P[tr2,3], family=binomial)
    cat(sprintf("%-18s %10.3f %10.3f\n", f,
      auc(as.numeric(cbind(1,P[te,2],P[te,3]) %*% coef(f1)), sm$dis[te]),
      auc(as.numeric(cbind(1,P[te,1],P[te,2],P[te,3]) %*% coef(f2)), sm$dis[te])))
  }
}
