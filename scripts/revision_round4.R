# ---------------------------------------------------------------------------
# revision_round4.R
# Analyses requested by the second pre-submission review.
# Francisco M. Codoner
# ---------------------------------------------------------------------------
suppressMessages({library(glmnet); library(vegan)})
src <- readLines("scripts/reported_values.R")
j <- grep("^folds *<-", src)[1]; k <- grep("^\\s*R\\$", src); i <- k[k > j][1] - 1
stopifnot(!is.na(j), !is.na(i), i > j)
eval(parse(text = paste(src[1:i], collapse = "\n")), envir = globalenv())
stopifnot(is.function(loco), is.function(fit))
mae <- function(p,a) mean(abs(p-a))
auc2 <- function(s,d){ r<-rank(s); n1<-sum(d==1); n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0) }

cat("##### A. PER-COHORT AGE-MODEL PERFORMANCE (new Supplementary Table) #####\n\n")
V <- readRDS("results/age_predictions.rds")$loco
tb <- do.call(rbind, lapply(sort(unique(V$cohort)), function(f){
  ii <- V$cohort==f; a <- V$actual[ii]; p <- V$pred[ii]
  data.frame(cohort=f, n=sum(ii), age_min=min(a), age_max=max(a),
             age_median=median(a), null_MAE=round(mean(abs(a-mean(V$actual[!ii]))),3),
             MAE=round(mae(p,a),3), r=round(cor(p,a),3))}))
tb$improvement_pct <- round(100*(tb$null_MAE-tb$MAE)/tb$null_MAE,1)
print(tb, row.names=FALSE)
write.csv(tb, "results/S7_per_cohort_age.csv", row.names=FALSE)
cat(sprintf("\nsample-weighted  MAE %.3f  r %.3f\n", mae(V$pred,V$actual), cor(V$pred,V$actual)))
cat(sprintf("cohort-weighted  MAE %.3f  (mean of per-cohort MAE)\n", mean(tb$MAE)))
cat(sprintf("median within-cohort r %.3f   mean %.3f\n", median(tb$r), mean(tb$r)))
cat(sprintf("age range of eligible cohorts: %d to %d years; excluded cohorts span <40 y\n",
            min(tb$age_min), max(tb$age_max)))

cat("\n##### B. WITHIN-LABORATORY vs CROSS-COHORT AT MATCHED TRAINING SIZE #####\n\n")
cat("isolates laboratory from training-set size: the cross-cohort model is\n")
cat("retrained on a random subset of OTHER cohorts of the same size.\n\n")
big <- names(nn)[nn>=200 & ag[names(nn)]>=40]
set.seed(11)
for (f in big) {
  ii <- which(co==f); ntr <- round(length(ii)*4/5)
  kk <- FF(length(ii)); pw <- numeric(length(ii))
  for (q in 1:5) { cv <- fit(XG[tr,][ii[kk!=q],,drop=FALSE], y[ii[kk!=q]])
    pw[kk==q] <- as.numeric(predict(cv, XG[tr,][ii[kk==q],,drop=FALSE], s="lambda.min")) }
  out <- which(co!=f); rs <- replicate(5, {
    s <- sample(out, ntr); cv <- fit(XG[tr,][s,,drop=FALSE], y[s])
    p <- as.numeric(predict(cv, XG[tr,][ii,,drop=FALSE], s="lambda.min"))
    c(mae(p,y[ii]), cor(p,y[ii])) })
  cat(sprintf("%-16s n=%3d  train=%3d | within-lab MAE %6.3f r %.3f | cross-cohort same size MAE %6.3f r %.3f\n",
      f, length(ii), ntr, mae(pw,y[ii]), cor(pw,y[ii]), mean(rs[1,]), mean(rs[2,])))
}

cat("\n##### C. DISPERSION ACROSS COHORTS (PERMANOVA assumption) #####\n\n")
RM <- readRDS("features/genus_relab.rds"); rownames(RM) <- unname(rownames(RM))
smr <- sm[match(rownames(RM), sm$run_accession), ]
pv <- colMeans(RM>0); RM <- RM[, pv>0.10, drop=FALSE]; RM <- RM/rowSums(RM)
cat(sprintf("samples %d, genera %d; computing Bray-Curtis...\n", nrow(RM), ncol(RM)))
D <- vegdist(RM, method="bray")
bd <- betadisper(D, factor(smr$cohort_id))
pt <- permutest(bd, permutations=199)
print(pt)
cat(sprintf("\nmean distance-to-centroid range across cohorts: %.3f to %.3f\n",
            min(tapply(bd$distances, bd$group, mean)), max(tapply(bd$distances, bd$group, mean))))

cat("\n##### D. DISEASE FOLDS: TRAINING COMPOSITION AND COUNTS #####\n\n")
sets <- list(OSCC=c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
             periodontitis=c("PORTUGAL_74247","SPAIN_74796"))
rowsD <- list()
for (g in names(sets)) for (f in sets[[g]]) {
  keep <- sm$cohort_id %in% sets[[g]] & (sm$dis==0 |
    grepl(if(g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  te <- which(keep & sm$cohort_id==f); t2 <- which(keep & sm$cohort_id!=f)
  if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
  rowsD[[paste(g,f)]] <- data.frame(condition=g, test_cohort=f,
    test_n=length(te), test_cases=sum(sm$dis[te]==1), test_controls=sum(sm$dis[te]==0),
    train_cohorts=length(unique(sm$cohort_id[t2])),
    train_cohort_names=paste(sort(unique(sm$cohort_id[t2])), collapse="; "),
    train_n=length(t2), train_cases=sum(sm$dis[t2]==1), train_controls=sum(sm$dis[t2]==0))
}
DD <- do.call(rbind, rowsD); print(DD[,1:10], row.names=FALSE)
write.csv(DD, "results/S8_disease_folds.csv", row.names=FALSE)

cat("\n##### E. DISEASE MODELS WITH FILTERS REFIT INSIDE TRAINING FOLDS #####\n\n")
for (g in names(sets)) for (f in sets[[g]]) {
  keep <- sm$cohort_id %in% sets[[g]] & (sm$dis==0 |
    grepl(if(g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  te <- which(keep & sm$cohort_id==f); t2 <- which(keep & sm$cohort_id!=f)
  if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
  cv <- cv.glmnet(XG[t2,], sm$dis[t2], family="binomial", alpha=.5, foldid=FF(length(t2)))
  a0 <- auc2(as.numeric(predict(cv, XG[te,], s="lambda.min")), sm$dis[te])
  Kin <- setdiff(colnames(G)[colMeans(G[t2,,drop=FALSE]>0)>0.10], CT)
  Xi <- clr(G[,Kin,drop=FALSE])
  cv2 <- cv.glmnet(Xi[t2,], sm$dis[t2], family="binomial", alpha=.5, foldid=FF(length(t2)))
  a1 <- auc2(as.numeric(predict(cv2, Xi[te,], s="lambda.min")), sm$dis[te])
  cat(sprintf("%-28s global filter AUC %.3f | fold-internal (%d genera) AUC %.3f | diff %+.3f\n",
              paste(g,f), a0, length(Kin), a1, a1-a0))
}

cat("\n##### F. MATCHING AND RAREFACTION DIAGNOSTICS #####\n\n")
RC <- readRDS("results/disease_roc.rds"); set.seed(1)
Mr <- G[sm$depth>=10000,,drop=FALSE]
for (nm in names(RC)) {
  s<-RC[[nm]]$score; d<-RC[[nm]]$label; agv<-RC[[nm]]$age
  ca<-which(d==1); ct<-which(d==0); mn<-c(); md<-c(); ma<-c()
  for (rep in 1:200){ pool<-ct; pc<-c()
    for (q in sample(ca)){ cand<-pool[abs(agv[pool]-agv[q])<=5]; if(!length(cand)) next
      jj<-if(length(cand)==1) cand else sample(cand,1); pc<-c(pc,q,jj); pool<-setdiff(pool,jj)}
    if (length(pc)>=10){ mn<-c(mn,length(pc))
      md<-c(md, mean(agv[pc[d[pc]==1]])-mean(agv[pc[d[pc]==0]])); ma<-c(ma,auc2(s[pc],d[pc])) } }
  coh <- sub("^(OSCC|periodontitis) ","",nm)
  idx <- which(sm$cohort_id==coh)
  cat(sprintf("%-28s unmatched age diff %5.1f y | matched n=%3.0f  age diff %+4.1f y  AUC %.3f [%.3f, %.3f]\n",
    nm, mean(agv[d==1])-mean(agv[d==0]), median(mn), median(md), median(ma),
    quantile(ma,.25), quantile(ma,.75)))
  cat(sprintf("%-28s rarefied: %d of %d samples retained at depth >= 10,000\n", "",
    sum(rownames(G)[idx] %in% rownames(Mr)), length(idx)))
}

cat("\n##### G. GRADIENT REFERENCE: PROVENANCE AND IN-SAMPLE FRACTION #####\n\n")
ib <- which(sm$cohort_id %in% c("PORTUGAL_74247","SPAIN_74796") &
            (sm$dis==0 | grepl("periodontitis", sm$label, ignore.case=TRUE)))
cat(sprintf("periodontitis classifier trained on %d samples from %d cohorts (%d cases, %d controls)\n",
   length(ib), length(unique(sm$cohort_id[ib])), sum(sm$dis[ib]==1), sum(sm$dis[ib]==0)))
insample <- intersect(ib, tr)
cat(sprintf("reference distribution: %d samples (all training-set samples), scored by that classifier\n", length(tr)))
cat(sprintf("of these, %d (%.1f%%) were in the classifier's own training set\n",
            length(insample), 100*length(insample)/length(tr)))
GD <- readRDS("results/gradient_detail.rds")
hs <- GD$healthy; keep2 <- setdiff(seq_along(tr), match(insample, tr))
cat(sprintf("90th percentile of reference: all %.3f | excluding in-sample %.3f\n",
            quantile(hs,.90), quantile(hs[keep2],.90)))
lv <- levels(GD$grade)
for (g in lv) cat(sprintf("  %-9s flagged %5.1f%% (all ref) | %5.1f%% (in-sample excluded)\n", g,
  100*mean(GD$score[GD$grade==g] >= quantile(hs,.90)),
  100*mean(GD$score[GD$grade==g] >= quantile(hs[keep2],.90))))
cat("\nthreshold and flag-rate uncertainty, 2000 bootstrap draws of the 17 local controls:\n")
loc <- GD$local_ctrl; set.seed(3)
bq <- replicate(2000, quantile(sample(loc, length(loc), replace=TRUE), .90))
cat(sprintf("  local threshold %.3f [%.3f, %.3f]\n", quantile(loc,.90), quantile(bq,.025), quantile(bq,.975)))
for (g in lv) { fr <- sapply(bq, function(th) 100*mean(GD$score[GD$grade==g] >= th))
  cat(sprintf("  %-9s flagged %5.1f%% [%5.1f, %5.1f]\n", g,
      100*mean(GD$score[GD$grade==g] >= quantile(loc,.90)), quantile(fr,.025), quantile(fr,.975))) }
cat(sprintf("\nclinical grade counts: %s\n", paste(sprintf("%s n=%d", lv, as.integer(table(GD$grade))), collapse="  ")))
cat("\ndone\n")
