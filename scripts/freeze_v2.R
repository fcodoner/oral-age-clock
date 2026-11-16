# --------------------------------------------------------------------------
# freeze_v2.R
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
clrm <- function(C){ L <- log(C[,K,drop=FALSE]+0.5); L - rowMeans(L) }
X  <- clrm(M)
auc <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }

tr <- which(sm$class=="train")
pc <- prcomp(X[tr,], center=TRUE, scale.=FALSE)
P  <- scale(X, center=pc$center, scale=FALSE) %*% pc$rotation
cat(sprintf("axes from %d healthy samples, %d genera (contaminants removed)\n", length(tr), length(K)))

ib <- which(sm$cohort_id %in% c("PORTUGAL_74247","SPAIN_74796") &
            (sm$dis==0 | grepl("periodontitis", sm$label, ignore.case=TRUE)))
oof <- rep(NA_real_, length(ib))
for (f in c("PORTUGAL_74247","SPAIN_74796")) {
  te <- which(sm$cohort_id[ib]==f); t2 <- which(sm$cohort_id[ib]!=f)
  fit <- glm(sm$dis[ib][t2] ~ P[ib[t2],2] + P[ib[t2],3], family=binomial)
  oof[te] <- as.numeric(cbind(1, P[ib[te],2], P[ib[te],3]) %*% coef(fit))
}
cat(sprintf("cross-cohort AUC (PC2+PC3): %.3f\n", auc(oof, sm$dis[ib])))
final <- glm(sm$dis[ib] ~ P[ib,2] + P[ib,3], family=binomial)
cal   <- glm(sm$dis[ib] ~ oof, family=binomial)

p <- plogis(predict(cal, data.frame(oof=oof))); d <- sm$dis[ib]
cat(sprintf("\n%10s %8s %8s\n","threshold","sens","spec"))
for (t in c(0.3,0.5,0.7)) { pr <- p>=t
  cat(sprintf("%10.2f %8.3f %8.3f\n", t, sum(pr&d==1)/sum(d==1), sum(!pr&d==0)/sum(d==0))) }

rv  <- which(sm$cohort_id=="RELVAS21_EU")
sr  <- as.numeric(cbind(1, P[rv,2], P[rv,3]) %*% coef(final))
pr  <- plogis(predict(cal, data.frame(oof=sr)))
lab <- ifelse(sm$class[rv]=="train", "0_healthy", sub(".*Grade_","grade_",sm$label[rv]))
cat("\n=== EXTERNAL VALIDATION: Relvas21, never seen ===\n")
cat(sprintf("%-12s %4s %8s %7s %7s %7s\n","grade","n","mean p","%LOW","%INT","%ELEV"))
for (g in sort(unique(lab))) { i <- lab==g; q <- pr[i]
  cat(sprintf("%-12s %4d %8.3f %6.0f%% %6.0f%% %6.0f%%\n", g, sum(i), mean(q),
      100*mean(q<0.30), 100*mean(q>=0.30&q<0.70), 100*mean(q>=0.70))) }
o <- as.integer(factor(lab))
cat(sprintf("\nSpearman vs grade: %.3f (p %.3g) | AUC healthy vs disease: %.3f\n",
    cor(pr,o,method="spearman"), suppressWarnings(cor.test(pr,o,method="spearman")$p.value),
    auc(pr, as.integer(lab!="0_healthy"))))

agem <- cv.glmnet(X[tr,], sm$age[tr], alpha=0.5, nfolds=5, standardize=TRUE)
saveRDS(list(genera=K, contaminants=CONT, pca=pc, oh_model=final, oh_calib=cal,
             age_model=agem, axis_names=c(PC2="periodontal pathogen axis",
             PC3="commensal / aerobe axis"), ref=list(mean=colMeans(X[tr,]), sd=apply(X[tr,],2,sd)),
             loadings=pc$rotation[,2:3], built=Sys.Date()), "models/oralclock_v2.rds")
cat("\nsaved: models/oralclock_v2.rds\n")
