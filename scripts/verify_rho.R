# --------------------------------------------------------------------------
# verify_rho.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages(library(glmnet))
CT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
        "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
        "Burkholderia","Cutibacterium","Escherichia")
G  <- readRDS("features/genus_counts.rds"); rownames(G) <- unname(rownames(G))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(G), sm$run_accession), ]
sm$dis <- ifelse(sm$class=="train", 0L, 1L)
rv  <- which(sm$cohort_id=="RELVAS21_EU")
lab <- ifelse(sm$class[rv]=="train","grade 0", sub(".*Grade_","grade ", sm$label[rv]))
ord <- as.integer(factor(lab, levels=c("grade 0","grade 1","grade 2","grade 3")))
dis <- as.integer(lab!="grade 0")
auc <- function(s,d){r<-rank(s);n1<-sum(d==1);n0<-sum(d==0)
  (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0)}

fitrho <- function(drop, seed=NULL, foldid=NULL) {
  K <- setdiff(colnames(G)[colMeans(G>0)>0.10], drop)
  L <- log(G[,K,drop=FALSE]+0.5); X <- L - rowMeans(L)
  ib <- which(sm$cohort_id %in% c("PORTUGAL_74247","SPAIN_74796") &
        (sm$dis==0 | grepl("periodontitis", sm$label, ignore.case=TRUE)))
  if (!is.null(seed)) set.seed(seed)
  cv <- if (is.null(foldid)) cv.glmnet(X[ib,], sm$dis[ib], family="binomial", alpha=.5, nfolds=5)
        else cv.glmnet(X[ib,], sm$dis[ib], family="binomial", alpha=.5, foldid=foldid)
  s <- as.numeric(predict(cv, X[rv,], s="lambda.min"))
  c(rho=cor(s, ord, method="spearman"), auc=auc(s, dis), nfeat=sum(as.matrix(coef(cv,s="lambda.min"))[-1,1]!=0))
}

for (nm in c("with contaminants","contaminants removed")) {
  drop <- if (nm=="with contaminants") character(0) else CT
  r <- t(sapply(1:200, function(i) fitrho(drop, seed=i)))
  cat(sprintf("\n=== %s ===\n", nm))
  cat(sprintf("rho : median %.3f | mean %.3f | sd %.3f | range %.3f - %.3f\n",
      median(r[,"rho"]), mean(r[,"rho"]), sd(r[,"rho"]), min(r[,"rho"]), max(r[,"rho"])))
  cat(sprintf("AUC : median %.3f | range %.3f - %.3f\n",
      median(r[,"auc"]), min(r[,"auc"]), max(r[,"auc"])))
  cat(sprintf("features selected: median %d | range %d - %d\n",
      median(r[,"nfeat"]), min(r[,"nfeat"]), max(r[,"nfeat"])))
}

cat("\n=== DETERMINISTIC PROCEDURE (fixed fold assignment) ===\n")
ib <- which(sm$cohort_id %in% c("PORTUGAL_74247","SPAIN_74796") &
      (sm$dis==0 | grepl("periodontitis", sm$label, ignore.case=TRUE)))
fid <- rep(1:5, length.out=length(ib))      # fixed, no RNG
d <- fitrho(CT, foldid=fid)
ct <- suppressWarnings(cor.test(
  { K <- setdiff(colnames(G)[colMeans(G>0)>0.10], CT)
    L <- log(G[,K,drop=FALSE]+0.5); X <- L - rowMeans(L)
    cv <- cv.glmnet(X[ib,], sm$dis[ib], family="binomial", alpha=.5, foldid=fid)
    as.numeric(predict(cv, X[rv,], s="lambda.min")) }, ord, method="spearman"))
cat(sprintf("rho = %.3f  (P = %.2e)  AUC = %.3f  features = %d\n",
    d["rho"], ct$p.value, d["auc"], d["nfeat"]))
cat("\nThis is reproducible on any machine: report it, and use the same foldid in the figure.\n")
