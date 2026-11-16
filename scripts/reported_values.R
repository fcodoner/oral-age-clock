# --------------------------------------------------------------------------
# reported_values.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages({library(glmnet); library(vegan)})
dir.create("results", showWarnings=FALSE)
FF <- function(n) rep(1:5, length.out=n)          # fixed folds, no RNG
R  <- list()
CT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
        "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
        "Burkholderia","Cutibacterium","Escherichia")
G  <- readRDS("features/genus_counts.rds"); rownames(G) <- unname(rownames(G))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(G), sm$run_accession), ]
sm$dis <- ifelse(sm$class=="train", 0L, 1L); sm$depth <- rowSums(G)
clr <- function(M){L <- log(M+0.5); L-rowMeans(L)}
K  <- setdiff(colnames(G)[colMeans(G>0)>0.10], CT)
XG <- clr(G[,K,drop=FALSE])
auc <- function(s,d){r<-rank(s);n1<-sum(d==1);n0<-sum(d==0)
  if(n1<2||n0<2) return(NA_real_); (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0)}
fit <- function(X,y,fam="gaussian") cv.glmnet(X,y,family=fam,alpha=.5,foldid=FF(nrow(X)),
         standardize=TRUE)
tr <- which(sm$class=="train"); y <- sm$age[tr]; co <- sm$cohort_id[tr]
ag <- tapply(y,co,function(v) diff(range(v))); nn <- table(co)
folds <- names(ag)[ag>=40 & nn[names(ag)]>=40]

loco <- function(X, yy, cc, ff) {
  P<-c();A<-c();C<-c()
  for (f in ff){ te<-cc==f
    cv<-fit(X[!te,,drop=FALSE],yy[!te])
    P<-c(P,as.numeric(predict(cv,X[te,,drop=FALSE],s="lambda.min")));A<-c(A,yy[te]);C<-c(C,rep(f,sum(te)))}
  wr<-sapply(unique(C),function(f){i<-C==f;suppressWarnings(cor(P[i],A[i]))})
  nu<-sapply(ff,function(f) mean(abs(yy[cc==f]-mean(yy[cc!=f]))))
  .LAST_PRED <<- list(pred=P, actual=A, cohort=C)
  c(MAE=mean(abs(P-A)), nullMAE=weighted.mean(nu,table(C)[ff]), r=cor(P,A),
    R2=1-sum((P-A)^2)/sum((A-mean(A))^2), medWithinR=median(wr,na.rm=TRUE))
}
cat("[1/7] age, leave-one-cohort-out\n")
R$age_loco <- loco(XG[tr,], y, co, folds)
PRED_loco <- .LAST_PRED
cat("    LOCO folds used:", length(folds), "of", length(unique(co)), "cohorts\n")

cat("[2/7] age, random cross-validation\n")
k <- FF(length(tr)); Q <- numeric(length(tr))
for (i in 1:5){ cv<-fit(XG[tr,][k!=i,],y[k!=i]); Q[k==i]<-as.numeric(predict(cv,XG[tr,][k==i,],s="lambda.min")) }
R$age_randomcv <- c(MAE=mean(abs(Q-y)), r=cor(Q,y), R2=1-sum((Q-y)^2)/sum((y-mean(y))^2))
PRED_rand <- list(pred=Q, actual=y, cohort=co)
saveRDS(list(loco=PRED_loco, rand=PRED_rand, folds=folds), "results/age_predictions.rds")
cat("    saved results/age_predictions.rds\n")

cat("[3/7] age, feature-space variants\n")
P <- t(as.matrix(read.delim(gzfile("picrust/pathways/path_abun_unstrat.tsv.gz"),
      row.names=1, check.names=FALSE)))
ids <- intersect(rownames(P), rownames(G))
XP <- clr(P[ids, colMeans(P[ids,]>0)>0.10, drop=FALSE])
mt <- match(ids, rownames(G)); trp <- which(sm$class[mt]=="train")
yp <- sm$age[mt][trp]; cp <- sm$cohort_id[mt][trp]
fp <- intersect(folds, unique(cp))
R$age_pathway <- loco(XP[trp,], yp, cp, fp)
R$age_both    <- loco(cbind(XG[mt,][trp,], XP[trp,]), yp, cp, fp)
zc <- function(X,cc){o<-X; for(c in unique(cc)){i<-cc==c
  m<-colMeans(X[i,,drop=FALSE]); s<-apply(X[i,,drop=FALSE],2,sd); s[s==0]<-1
  o[i,]<-sweep(sweep(X[i,,drop=FALSE],2,m,"-"),2,s,"/")}; o}
R$age_zscore  <- loco(zc(XG[tr,],co), y, co, folds)
R$age_noAGP   <- { ff<-setdiff(folds,"AGP_ORAL"); i<-co!="AGP_ORAL"
                   loco(XG[tr,][i,], y[i], co[i], intersect(ff,unique(co[i]))) }

cat("[4/7] age, ASV vs genus (V3-V4 paired)\n")
EX <- "CHAUDHARI20_IN__paired__V3-V4_seqtab.rds"
ffs <- list.files("final_split/train", pattern="paired__V3-V4_seqtab.rds$", full.names=TRUE)
ffs <- ffs[basename(ffs)!=EX]
tx <- readRDS("taxonomy/taxa_oriented.rds")$tax
bb <- !is.na(tx[,"Phylum"]) & !is.na(tx[,"Kingdom"]) & tx[,"Kingdom"]=="Bacteria"
names(bb) <- rownames(tx)
pr <- list(); nt <- 0L
for (p in ffs){ st<-readRDS(p); st<-st[,bb[colnames(st)],drop=FALSE]
  st<-st[rowSums(st)>=3000,,drop=FALSE]; nt<-nt+nrow(st); pr[[p]]<-colSums(st>0) }
v <- unlist(pr,use.names=TRUE); names(v) <- sub("^.*\\.(?=[ACGT])","",names(v),perl=TRUE)
pv <- tapply(v,names(v),sum); sel <- names(pv)[pv >= 0.05*nt]
A <- NULL; md <- list()
for (p in ffs){ st<-readRDS(p); st<-st[,bb[colnames(st)],drop=FALSE]
  st<-st[rowSums(st)>=3000,,drop=FALSE]
  m<-matrix(0L,nrow(st),length(sel),dimnames=list(rownames(st),sel))
  cm<-intersect(colnames(st),sel); m[,cm]<-st[,cm,drop=FALSE]; A<-rbind(A,m)
  md[[p]]<-data.frame(id=rownames(st),coh=strsplit(sub("_seqtab.rds$","",basename(p)),"__")[[1]][1]) }
mdd <- do.call(rbind, md); ya <- sm$age[match(mdd$id, sm$run_accession)]
aga <- tapply(ya,mdd$coh,function(z) diff(range(z))); na_ <- table(mdd$coh)
fa <- names(aga)[aga>=40 & na_[names(aga)]>=40]
R$asv_level   <- loco(clr(A), ya, mdd$coh, fa)
R$genus_sameF <- loco(XG[match(mdd$id, rownames(G)),], ya, mdd$coh, fa)

cat("[5/7] age, within laboratory\n")
big <- names(nn)[nn>=200 & ag[names(nn)]>=40]
R$age_withinlab <- t(sapply(big, function(f){
  i<-which(co==f); kk<-FF(length(i)); p<-numeric(length(i))
  for (j in 1:5){cv<-fit(XG[tr,][i[kk!=j],],y[i[kk!=j]])
    p[kk==j]<-as.numeric(predict(cv,XG[tr,][i[kk==j],],s="lambda.min"))}
  c(n=length(i), nullMAE=mean(abs(y[i]-mean(y[i]))), MAE=mean(abs(p-y[i])),
    r=cor(p,y[i]), R2=1-sum((p-y[i])^2)/sum((y[i]-mean(y[i]))^2))}))

cat("[6/7] disease, cross-cohort with age and depth controls\n")
set.seed(1)
Mr <- G[sm$depth>=10000,,drop=FALSE]; Mr <- rrarefy(Mr,10000); Xr <- clr(Mr[,K,drop=FALSE])
sets <- list(OSCC=c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
             periodontitis=c("PORTUGAL_74247","SPAIN_74796"))
out <- list()
ROC <- list()
for (g in names(sets)) for (f in sets[[g]]) {
  keep <- sm$cohort_id %in% sets[[g]] & (sm$dis==0 |
    grepl(if(g=="OSCC") "OSCC|carcinoma" else g, sm$label, ignore.case=TRUE))
  te<-which(keep&sm$cohort_id==f); t2<-which(keep&sm$cohort_id!=f)
  if (sum(sm$dis[te]==1)<5||sum(sm$dis[te]==0)<5) next
  d<-sm$dis[te]
  cv<-cv.glmnet(XG[t2,],sm$dis[t2],family="binomial",alpha=.5,foldid=FF(length(t2)))
  s<-as.numeric(predict(cv,XG[te,],s="lambda.min")); sz<-(s-mean(s))/sd(s)
  cf<-summary(glm(d~sm$age[te]+sz,family=binomial))$coefficients
  ad<-auc(sm$depth[te],d); ad<-max(ad,1-ad)
  set.seed(1); ca<-which(d==1); ct_<-which(d==0); acc<-c()
  for (rep in 1:200){ pool<-ct_; pc<-c(); for (i in sample(ca)){
      cand<-pool[abs(sm$age[te][pool]-sm$age[te][i])<=5]; if(!length(cand)) next
      j<-if(length(cand)==1) cand else sample(cand,1); pc<-c(pc,i,j); pool<-setdiff(pool,j)}
    if (length(pc)>=10) acc<-c(acc,auc(s[pc],d[pc])) }
  ter<-te[rownames(G)[te] %in% rownames(Mr)]; trr<-t2[rownames(G)[t2] %in% rownames(Mr)]
  ar<-NA_real_
  if (length(ter)>=10 && sum(sm$dis[ter]==1)>=5 && sum(sm$dis[ter]==0)>=5){
    i1<-match(rownames(G)[trr],rownames(Mr)); i2<-match(rownames(G)[ter],rownames(Mr))
    c2<-cv.glmnet(Xr[i1,],sm$dis[trr],family="binomial",alpha=.5,foldid=FF(length(i1)))
    ar<-auc(as.numeric(predict(c2,Xr[i2,],s="lambda.min")),sm$dis[ter]) }
  ROC[[paste(g,f)]] <- list(score=s, label=d, age=sm$age[te])
  out[[paste(g,f)]] <- c(n=length(te), cases=sum(d), AUC=auc(s,d),
    AUC_age=auc(sm$age[te],d), OR_SD=exp(cf["sz","Estimate"]), p=cf["sz","Pr(>|z|)"],
    AUC_matched=mean(acc,na.rm=TRUE), AUC_depth=ad, AUC_rarefied=ar)
}
R$disease <- do.call(rbind,out)
saveRDS(ROC, "results/disease_roc.rds")

cat("[7/7] oral-health gradient\n")
ib<-which(sm$cohort_id %in% c("PORTUGAL_74247","SPAIN_74796") &
   (sm$dis==0|grepl("periodontitis",sm$label,ignore.case=TRUE)))
cv<-cv.glmnet(XG[ib,],sm$dis[ib],family="binomial",alpha=.5,foldid=FF(length(ib)))
rv<-which(sm$cohort_id=="RELVAS21_EU")
s<-as.numeric(predict(cv,XG[rv,],s="lambda.min"))
lab<-factor(ifelse(sm$class[rv]=="train","grade 0",sub(".*Grade_","grade ",sm$label[rv])),
  levels=c("grade 0","grade 1","grade 2","grade 3"))
ct_<-suppressWarnings(cor.test(s,as.integer(lab),method="spearman"))
d23<-lab %in% c("grade 2","grade 3"); kk<-lab %in% c("grade 0","grade 2","grade 3")
hs<-as.numeric(predict(cv,XG[tr,],s="lambda.min")); loc<-s[lab=="grade 0"]
R$gradient <- c(rho=ct_$estimate, p=ct_$p.value, AUC_any=auc(s,as.integer(lab!="grade 0")),
  AUC_g3=auc(s[lab %in% c("grade 0","grade 3")],as.integer(lab[lab %in% c("grade 0","grade 3")]=="grade 3")),
  AUC_g23=auc(s[kk],as.integer(d23[kk])), nfeat=sum(as.matrix(coef(cv,s="lambda.min"))[-1,1]!=0))
R$anchoring <- rbind(published=sapply(levels(lab),function(g) 100*mean(s[lab==g]>=quantile(hs,.90))),
                     local    =sapply(levels(lab),function(g) 100*mean(s[lab==g]>=quantile(loc,.90))))
saveRDS(list(score=s, grade=lab, healthy=hs, local_ctrl=loc,
  thr_published=quantile(hs,.90), thr_local=quantile(loc,.90)),
  "results/gradient_detail.rds")
saveRDS(R,"results/reported_values.rds")

cat("\n================ REPORTED VALUES ================\n")
pr2 <- function(nm,v) cat(sprintf("%-22s %s\n", nm,
  paste(sprintf("%s=%.3f", names(v), v), collapse="  ")))
for (nm in c("age_loco","age_randomcv","age_pathway","age_both","age_zscore","age_noAGP",
             "asv_level","genus_sameF","gradient")) pr2(nm, R[[nm]])
cat("\n-- age within laboratory --\n"); print(round(R$age_withinlab,3))
cat("\n-- disease --\n"); print(round(R$disease,3))
cat("\n-- anchoring (% flagged elevated) --\n"); print(round(R$anchoring,1))
cat("\nsaved: results/reported_values.rds\n")
