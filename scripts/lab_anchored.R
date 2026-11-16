# --------------------------------------------------------------------------
# lab_anchored.R
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
sm$dis <- ifelse(sm$class=="train", 0L, 1L)
CONT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
          "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
          "Burkholderia","Cutibacterium","Escherichia")
K <- setdiff(colnames(M)[colMeans(M>0)>0.10], CONT)
L <- log(M[,K,drop=FALSE]+0.5); X <- L - rowMeans(L)
auc <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }

ib <- which(sm$cohort_id %in% c("PORTUGAL_74247","SPAIN_74796") &
            (sm$dis==0 | grepl("periodontitis", sm$label, ignore.case=TRUE)))
cv <- cv.glmnet(X[ib,], sm$dis[ib], family="binomial", alpha=0.5, nfolds=5)
tr <- which(sm$class=="train")
pc <- prcomp(X[tr,], center=TRUE, scale.=FALSE)
P  <- scale(X, center=pc$center, scale=FALSE) %*% pc$rotation
ax <- glm(sm$dis[ib] ~ P[ib,2] + P[ib,3], family=binomial)

rv  <- which(sm$cohort_id=="RELVAS21_EU")
lab <- ifelse(sm$class[rv]=="train","0_healthy", sub(".*Grade_","grade_",sm$label[rv]))
cands <- list(
  `elastic net` = as.numeric(predict(cv, X[rv,], s="lambda.min")),
  `PC2+PC3`     = as.numeric(cbind(1,P[rv,2],P[rv,3]) %*% coef(ax)))

for (nm in names(cands)) {
  s <- cands[[nm]]; ref <- s[lab=="0_healthy"]
  q50 <- quantile(ref,0.50); q90 <- quantile(ref,0.90)
  band <- ifelse(s>=q90,"ELEV", ifelse(s>=q50,"INT","LOW"))
  pct  <- sapply(s, function(z) mean(ref <= z))
  cat(sprintf("\n=== %s, anchored on this lab's own %d healthy controls ===\n", nm, length(ref)))
  cat(sprintf("%-12s %4s %11s %7s %7s %7s\n","grade","n","med pctile","%LOW","%INT","%ELEV"))
  for (g in sort(unique(lab))) { i <- lab==g
    cat(sprintf("%-12s %4d %10.0f%% %6.0f%% %6.0f%% %6.0f%%\n", g, sum(i),
        100*median(pct[i]), 100*mean(band[i]=="LOW"),
        100*mean(band[i]=="INT"), 100*mean(band[i]=="ELEV"))) }
  d2 <- as.integer(lab %in% c("grade_2","grade_3"))
  k  <- lab %in% c("0_healthy","grade_2","grade_3")
  cat(sprintf("  AUC healthy vs grade2+3: %.3f | %% of grade2+3 flagged ELEV: %.0f%%\n",
      auc(s[k], d2[k]), 100*mean(band[lab %in% c("grade_2","grade_3")]=="ELEV")))
}
cat("\n(reference n=17 is small; this simulates a clinic running its own control panel)\n")
