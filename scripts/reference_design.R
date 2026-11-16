# --------------------------------------------------------------------------
# reference_design.R
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
B <- readRDS("models/oralclock_final.rds")
L <- log(M[,B$genera,drop=FALSE]+0.5); X <- L - rowMeans(L)
auc <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }
sc <- as.numeric(predict(B$model, X, s="lambda.min"))

NITRATE <- intersect(c("Rothia","Neisseria","Haemophilus","Actinomyces",
                       "Veillonella","Granulicatella"), colnames(X))
ni <- rowMeans(X[, NITRATE, drop=FALSE])
cat("nitrate-reducer panel:", paste(NITRATE, collapse=", "), "\n")

rv  <- which(sm$cohort_id=="RELVAS21_EU")
lab <- ifelse(sm$class[rv]=="train","0_healthy", sub(".*Grade_","grade_",sm$label[rv]))
cat("\n--- does the panel fall with worsening grade? ---\n")
print(round(tapply(ni[rv], lab, median), 3))
cat(sprintf("Spearman panel vs grade: %.3f (negative = falls with disease)\n",
    cor(ni[rv], as.integer(factor(lab)), method="spearman")))

EUR <- c("SLL2_ES","TRIAD_ES","SPAIN_83278","SOUTHKORE_13236")
refs <- list(
  `local 17 (Relvas healthy)` = rv[lab=="0_healthy"],
  `matched: Iberian V3-V4`    = which(sm$class=="train" & sm$cohort_id %in% c("SLL2_ES","TRIAD_ES","SPAIN_83278")),
  `matched + Korea V3-V4`     = which(sm$class=="train" & sm$cohort_id %in% EUR),
  `global: all healthy`       = which(sm$class=="train"))

cat(sprintf("\n%-28s %5s %8s %8s %8s %8s %9s\n",
    "reference panel","n","%ELEV h","%ELEV g1","%ELEV g2","%ELEV g3","AUC h vs 2+3"))
for (nm in names(refs)) {
  r <- refs[[nm]]; q90 <- quantile(sc[r], 0.90); q50 <- quantile(sc[r], 0.50)
  b <- ifelse(sc[rv] >= q90, "ELEV", ifelse(sc[rv] >= q50, "INT", "LOW"))
  k <- lab %in% c("0_healthy","grade_2","grade_3")
  cat(sprintf("%-28s %5d %7.0f%% %7.0f%% %7.0f%% %7.0f%% %9.3f\n", nm, length(r),
      100*mean(b[lab=="0_healthy"]=="ELEV"), 100*mean(b[lab=="grade_1"]=="ELEV"),
      100*mean(b[lab=="grade_2"]=="ELEV"), 100*mean(b[lab=="grade_3"]=="ELEV"),
      auc(sc[rv][k], as.integer(lab[k]!="0_healthy"))))
}
cat("\n(AUC is reference-independent; only the banding changes. Target: ~10% of healthy\n",
    " flagged ELEV, with the highest possible proportion of grades 2-3 flagged.)\n")
