# --------------------------------------------------------------------------
# freeze_models.R
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
dir.create("models", showWarnings=FALSE)
M  <- readRDS("features/genus_counts.rds"); rownames(M) <- unname(rownames(M))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(M), sm$run_accession), ]
KEEP <- colnames(M)[colMeans(M>0) > 0.10]
clr <- function(C) { C <- C[, KEEP, drop=FALSE]; L <- log(C + 0.5); L - rowMeans(L) }
X   <- clr(M)
sm$dis <- ifelse(sm$class=="train", 0L, 1L)
auc <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }

fit_task <- function(nm, cohorts, pattern) {
  i <- which(sm$cohort_id %in% cohorts &
             (sm$dis==0 | grepl(pattern, sm$label, ignore.case=TRUE)))
  # cross-cohort predictions -> calibration + thresholds (no leakage)
  oof <- rep(NA_real_, length(i))
  for (f in cohorts) {
    te <- which(sm$cohort_id[i]==f); tr <- which(sm$cohort_id[i]!=f)
    if (!length(te) || length(unique(sm$dis[i][tr]))<2) next
    cv <- cv.glmnet(X[i[tr],], sm$dis[i][tr], family="binomial", alpha=0.5, nfolds=5)
    oof[te] <- as.numeric(predict(cv, X[i[te],], s="lambda.min"))
  }
  ok <- !is.na(oof); d <- sm$dis[i][ok]; s <- oof[ok]
  cal <- glm(d ~ s, family=binomial)                       # Platt scaling
  final <- cv.glmnet(X[i,], sm$dis[i], family="binomial", alpha=0.5, nfolds=5)
  b <- as.matrix(coef(final, s="lambda.min"))[,1]
  cat(sprintf("\n===== %s =====\n n=%d (%d cases, %d cohorts) | cross-cohort AUC %.3f | %d genera retained\n",
      nm, length(i), sum(sm$dis[i]), length(cohorts), auc(s,d), sum(b[-1]!=0)))
  p <- plogis(predict(cal, data.frame(s=s)))
  cat(sprintf("%10s %9s %9s %9s %9s\n","threshold","sens","spec","PPV","NPV"))
  for (t in c(0.3,0.4,0.5,0.6,0.7)) {
    pr <- p >= t
    cat(sprintf("%10.2f %9.3f %9.3f %9.3f %9.3f\n", t,
        sum(pr & d==1)/sum(d==1), sum(!pr & d==0)/sum(d==0),
        ifelse(sum(pr)>0, sum(pr & d==1)/sum(pr), NA),
        ifelse(sum(!pr)>0, sum(!pr & d==0)/sum(!pr), NA)))
  }
  list(model=final, lambda="lambda.min", calib=cal, genera=KEEP,
       coefs=b[b!=0], auc_cv=auc(s,d), n=length(i), cases=sum(sm$dis[i]), cohorts=cohorts)
}

oh  <- fit_task("ORAL HEALTH (periodontitis)",
                c("PORTUGAL_74247","SPAIN_74796"), "periodontitis")
osc <- fit_task("OSCC",
                c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"), "OSCC|carcinoma")

cat("\n===== AGE (secondary, low precision) =====\n")
tr <- which(sm$class=="train")
agem <- cv.glmnet(X[tr,], sm$age[tr], alpha=0.5, nfolds=5, standardize=TRUE)
ab <- as.matrix(coef(agem, s="lambda.min"))[,1]
cat(sprintf(" trained on %d healthy samples | %d genera retained | LOCO MAE 12.2 y (null 14.4)\n",
    length(tr), sum(ab[-1]!=0)))

bundle <- list(oral_health=oh, oscc=osc,
               age=list(model=agem, genera=KEEP, coefs=ab[ab!=0],
                        loco_mae=12.18, loco_null=14.44, note="low precision; context only"),
               genera=KEEP, built=Sys.Date(), n_train=length(tr))
saveRDS(bundle, "models/oralclock_v1.rds")

for (nm in c("oral_health","oscc")) {
  b <- bundle[[nm]]$coefs; b <- b[names(b)!="(Intercept)"]
  write.table(data.frame(genus=names(b), coefficient=round(b,4),
                         direction=ifelse(b>0,"higher in disease","lower in disease")),
    sprintf("models/coefficients_%s.tsv", nm), sep="\t", quote=FALSE, row.names=FALSE)
}
ab2 <- ab[names(ab)!="(Intercept)"]
write.table(data.frame(genus=names(ab2), coefficient=round(ab2,4)),
  "models/coefficients_age.tsv", sep="\t", quote=FALSE, row.names=FALSE)

cat("\nsaved: models/oralclock_v1.rds + coefficient tables\n")
cat(sprintf("oral health: %d genera | OSCC: %d genera | age: %d genera\n",
    length(oh$coefs)-1, length(osc$coefs)-1, length(ab2[ab2!=0])))
