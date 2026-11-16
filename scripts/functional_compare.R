# --------------------------------------------------------------------------
# functional_compare.R
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
P  <- read.delim(gzfile("picrust/pathways/path_abun_unstrat.tsv.gz"), row.names=1, check.names=FALSE)
P  <- t(as.matrix(P))
G  <- readRDS("features/genus_counts.rds"); rownames(G) <- unname(rownames(G))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
ids <- intersect(rownames(P), rownames(G))
P <- P[ids,]; G <- G[ids,]; sm <- sm[match(ids, sm$run_accession), ]
sm$dis <- ifelse(sm$class=="train", 0L, 1L)
CONT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
          "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
          "Burkholderia","Cutibacterium","Escherichia")
clr <- function(M){ L <- log(M + 0.5); L - rowMeans(L) }
XG <- clr(G[, setdiff(colnames(G)[colMeans(G>0)>0.10], CONT), drop=FALSE])
XP <- clr(P[, colMeans(P>0) > 0.10, drop=FALSE])
XB <- cbind(XG, XP)
cat(sprintf("samples %d | genera %d | pathways %d\n", length(ids), ncol(XG), ncol(XP)))
auc <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }

tr <- which(sm$class=="train")
y  <- sm$age[tr]; co <- sm$cohort_id[tr]
ag <- tapply(y, co, function(v) diff(range(v))); nn <- table(co)
folds <- names(ag)[ag >= 40 & nn[names(ag)] >= 40]
cat("\n=== A. AGE, leave-one-cohort-out (", length(folds), "folds ) ===\n")
cat(sprintf("%-12s %8s %8s %8s %9s %12s\n","features","MAE","nullMAE","r","R2","medWithinR"))
for (nm in c("genus","pathway","both")) {
  X <- switch(nm, genus=XG[tr,], pathway=XP[tr,], both=XB[tr,])
  p <- c(); a <- c(); cc <- c()
  for (f in folds) {
    te <- co==f
    cv <- cv.glmnet(X[!te,], y[!te], alpha=0.5, nfolds=5, standardize=TRUE)
    p <- c(p, as.numeric(predict(cv, X[te,], s="lambda.min")))
    a <- c(a, y[te]); cc <- c(cc, rep(f, sum(te)))
  }
  wr <- sapply(unique(cc), function(f){i<-cc==f; suppressWarnings(cor(p[i],a[i]))})
  nu <- sapply(folds, function(f) mean(abs(y[co==f] - mean(y[co!=f]))))
  cat(sprintf("%-12s %8.2f %8.2f %8.3f %9.3f %12.3f\n", nm, mean(abs(p-a)),
      weighted.mean(nu, table(cc)[folds]), cor(p,a),
      1-sum((p-a)^2)/sum((a-mean(a))^2), median(wr, na.rm=TRUE)))
}

cat("\n=== B. DISEASE, cross-cohort ===\n")
sets <- list(OSCC=c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
             periodontitis=c("PORTUGAL_74247","SPAIN_74796"))
for (g in names(sets)) {
  cs <- sets[[g]]
  keep <- sm$cohort_id %in% cs & (sm$dis==0 |
          grepl(if (g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  cat(sprintf("\n-- %s --\n%-18s %9s %9s %9s\n", g, "held-out","genus","pathway","both"))
  for (f in cs) {
    te <- which(keep & sm$cohort_id==f); t2 <- which(keep & sm$cohort_id!=f)
    if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
    o <- sapply(list(XG,XP,XB), function(X) {
      cv <- cv.glmnet(X[t2,], sm$dis[t2], family="binomial", alpha=0.5, nfolds=5)
      auc(as.numeric(predict(cv, X[te,], s="lambda.min")), sm$dis[te]) })
    cat(sprintf("%-18s %9.3f %9.3f %9.3f\n", f, o[1], o[2], o[3]))
  }
}

cat("\n=== C. WITHIN-LABORATORY AGE (5-fold CV) ===\n")
big <- names(nn)[nn >= 200 & ag[names(nn)] >= 40]
cat(sprintf("%-18s %5s %8s %9s %9s\n","cohort","n","nullMAE","genusMAE","pathMAE"))
for (f in big) {
  i <- which(co==f); k <- sample(rep(1:5, length.out=length(i)))
  o <- sapply(list(XG[tr,], XP[tr,]), function(X) {
    p <- numeric(length(i))
    for (j in 1:5) { cv <- cv.glmnet(X[i[k!=j],], y[i[k!=j]], alpha=0.5, nfolds=5, standardize=TRUE)
      p[k==j] <- as.numeric(predict(cv, X[i[k==j],], s="lambda.min")) }
    mean(abs(p - y[i])) })
  cat(sprintf("%-18s %5d %8.2f %9.2f %9.2f\n", f, length(i),
      mean(abs(y[i]-mean(y[i]))), o[1], o[2]))
}
