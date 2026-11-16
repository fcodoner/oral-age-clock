# --------------------------------------------------------------------------
# freeze_v3.R
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
auc <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }
mkX <- function(drop){ K <- setdiff(colnames(M)[colMeans(M>0)>0.10], drop)
  L <- log(M[,K,drop=FALSE]+0.5); list(X=L-rowMeans(L), K=K) }

ibc <- c("PORTUGAL_74247","SPAIN_74796")
ib  <- which(sm$cohort_id %in% ibc & (sm$dis==0 |
             grepl("periodontitis", sm$label, ignore.case=TRUE)))
rv  <- which(sm$cohort_id=="RELVAS21_EU")
lab <- ifelse(sm$class[rv]=="train","0_healthy", sub(".*Grade_","grade_",sm$label[rv]))
ord <- as.integer(factor(lab)); dis <- as.integer(lab!="0_healthy")

evaluate <- function(nm, raw_train, raw_rv, healthy_scores) {
  cutoff <- quantile(healthy_scores, 0.90)
  band   <- ifelse(raw_rv >= cutoff, "ELEV",
             ifelse(raw_rv >= quantile(healthy_scores, 0.50), "INT", "LOW"))
  pct <- sapply(raw_rv, function(z) mean(healthy_scores <= z))
  cat(sprintf("\n--- %s ---\n  Relvas AUC %.3f | Spearman %.3f\n", nm, auc(raw_rv, dis),
      cor(raw_rv, ord, method="spearman")))
  cat(sprintf("  %-12s %4s %10s %7s %7s %7s\n","grade","n","med pctile","%LOW","%INT","%ELEV"))
  for (g in sort(unique(lab))) { i <- lab==g
    cat(sprintf("  %-12s %4d %9.0f%% %6.0f%% %6.0f%% %6.0f%%\n", g, sum(i),
        100*median(pct[i]), 100*mean(band[i]=="LOW"),
        100*mean(band[i]=="INT"), 100*mean(band[i]=="ELEV"))) }
}

for (nm in c("with contaminants","contaminants removed")) {
  o <- mkX(if (nm=="with contaminants") character(0) else CONT)
  X <- o$X
  cv <- cv.glmnet(X[ib,], sm$dis[ib], family="binomial", alpha=0.5, nfolds=5)
  hs <- as.numeric(predict(cv, X[which(sm$class=="train"),], s="lambda.min"))
  evaluate(sprintf("elastic net, %s (%d genera, %d selected)", nm, length(o$K),
                   sum(as.matrix(coef(cv,s="lambda.min"))[-1,1]!=0)),
           NULL, as.numeric(predict(cv, X[rv,], s="lambda.min")), hs)
  if (nm=="contaminants removed") BEST <- list(model=cv, K=o$K, healthy=hs, X=X)
}

o  <- mkX(CONT); X <- o$X; tr <- which(sm$class=="train")
pc <- prcomp(X[tr,], center=TRUE, scale.=FALSE)
P  <- scale(X, center=pc$center, scale=FALSE) %*% pc$rotation
f  <- glm(sm$dis[ib] ~ P[ib,2] + P[ib,3], family=binomial)
evaluate("PC2+PC3 axes", NULL,
         as.numeric(cbind(1,P[rv,2],P[rv,3]) %*% coef(f)),
         as.numeric(cbind(1,P[tr,2],P[tr,3]) %*% coef(f)))

saveRDS(list(genera=BEST$K, contaminants=CONT, model=BEST$model,
             healthy_ref=BEST$healthy, pca=pc, axis_model=f,
             q50=quantile(BEST$healthy,0.50), q90=quantile(BEST$healthy,0.90),
             loadings=pc$rotation[,2:3], ref=list(mean=colMeans(X[tr,]), sd=apply(X[tr,],2,sd)),
             built=Sys.Date()), "models/oralclock_v3.rds")
cat("\nsaved: models/oralclock_v3.rds (percentile-anchored banding)\n")
