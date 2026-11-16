# --------------------------------------------------------------------------
# figures.R
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
dir.create("figures", showWarnings=FALSE)
CT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
        "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
        "Burkholderia","Cutibacterium","Escherichia")
G  <- readRDS("features/genus_counts.rds"); rownames(G) <- unname(rownames(G))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(G), sm$run_accession), ]
sm$dis <- ifelse(sm$class=="train", 0L, 1L)
K  <- setdiff(colnames(G)[colMeans(G>0)>0.10], CT)
L  <- log(G[,K,drop=FALSE]+0.5); X <- L - rowMeans(L)
tr <- which(sm$class=="train"); y <- sm$age[tr]; co <- sm$cohort_id[tr]
ag <- tapply(y, co, function(v) diff(range(v))); nn <- table(co)
folds <- names(ag)[ag>=40 & nn[names(ag)]>=40]
PAL <- c("#2E74B5","#C62828","#2E7D32","#E8871A","#6A4C93","#777777")
RV <- readRDS("results/reported_values.rds")   # single source of truth
op <- function(f,w,h){pdf(file.path("figures",f),width=w,height=h,pointsize=9)
  par(family="sans")}

## ---- Figure 1 : database ----
op("Fig1.pdf",7.2,7.0); layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
par(mar=c(4,8,2,1))
st <- c("Candidate runs"=5685,"Status assigned"=4491,"Denoised"=4437,"Analysed"=4426)
barplot(rev(st),horiz=TRUE,las=1,col=PAL[1],border=NA,xlab="samples",main="a  Assembly")
par(mar=c(4,4,2,1))
hist(sm$age[tr],breaks=30,col=PAL[1],border="white",main="b  Age, healthy set",
     xlab="age (years)"); abline(v=median(sm$age[tr]),lty=2)
par(mar=c(8,4,2,1))
R <- G/rowSums(G); m <- sort(colMeans(R[tr,]),decreasing=TRUE)[1:15]
barplot(100*m,las=2,col=PAL[1],border=NA,ylab="mean relative abundance (%)",
        main="c  Mean composition, healthy set",cex.names=.8)
dev.off()

## ---- Figure 2 : batch and evaluation ----
op("Fig2.pdf",7.2,3.6); par(mfrow=c(1,2),mar=c(4,4,2.5,1))
v <- rbind(`batch first`=c(0.14,0.05,6.15,1.70,27.69),
           `age first`  =c(0.84,2.30,5.71,2.15,24.74))
barplot(v,beside=TRUE,names.arg=c("age","sex","region","route","cohort"),
        col=PAL[c(1,2)],border=NA,ylab="variance explained (%)",las=1,
        main="a  Variance partition"); legend("topleft",rownames(v),fill=PAL[c(1,2)],
        border=NA,bty="n",cex=.85)
P <- c(); A <- c()
for (f in folds) { te <- co==f
  cv <- cv.glmnet(X[tr,][!te,],y[!te],alpha=.5,nfolds=5,standardize=TRUE)
  P <- c(P,as.numeric(predict(cv,X[tr,][te,],s="lambda.min"))); A <- c(A,y[te]) }
k <- sample(rep(1:5,length.out=length(tr))); Q <- numeric(length(tr))
for (i in 1:5){ cv <- cv.glmnet(X[tr,][k!=i,],y[k!=i],alpha=.5,nfolds=5,standardize=TRUE)
  Q[k==i] <- as.numeric(predict(cv,X[tr,][k==i,],s="lambda.min")) }
plot(y,Q,pch=16,col=adjustcolor(PAL[3],.25),cex=.5,xlab="chronological age (years)",
     ylab="predicted age (years)",las=1,main="b  Evaluation scheme")
points(A,P,pch=16,col=adjustcolor(PAL[2],.25),cex=.5); abline(0,1,lty=2)
legend("topleft",c(sprintf("random CV  MAE %.1f",mean(abs(Q-y))),
       sprintf("LOCO  MAE %.1f",mean(abs(P-A)))),col=PAL[c(3,2)],pch=16,bty="n",cex=.85)
dev.off()
saveRDS(list(P=P,A=A,Q=Q,y=y),"figures/fig2_data.rds")

## ---- Figure 3 : ASV vs genus ----
op("Fig3.pdf",4.6,4.2); par(mar=c(4,4,2.5,1))
gm <- RV$genus_sameF["medWithinR"]; am <- RV$asv_level["medWithinR"]
plot(range(0,.55),c(.8,2.2),type="n",axes=FALSE,xlab="within-cohort correlation with age",
     ylab="",main="Taxonomic resolution")
axis(1); axis(2,at=c(1,2),labels=c("ASV","genus"),las=1)
segments(am,1,gm,2,col=adjustcolor("grey40",.6),lwd=2)
points(am,1,pch=16,col=PAL[2],cex=1.6); points(gm,2,pch=16,col=PAL[1],cex=1.6)
text(am,.88,sprintf("%.3f",am),cex=.9); text(gm,2.12,sprintf("%.3f",gm),cex=.9)
dev.off()

## ---- Figure 4 : meta-analysis ----
op("Fig4.pdf",7.2,4.6); par(mfrow=c(1,2),mar=c(4,11,2.5,1))
gen <- c(Alloprevotella=-.158,Haemophilus=-.157,Sneathia=-.134,Lautropia=-.129,
  `Absconditabacteria SR1`=-.109,Bergeyella=-.090,`Saccharibacteria TM7`=-.071,
  Lactobacillus=.080,Desulfovibrio=.085,`Bacteroidaceae [G-1]`=.086,
  `Veillonellaceae [G-1]`=.088,Anaeroglobus=.098,Bifidobacterium=.119)
pw <- c(`heme b biosynthesis (uroporph.)`=-.124,`heme b biosynthesis II`=-.123,
  `heme b biosynthesis V`=-.120,`(Kdo)2-lipid A biosynthesis`=-.138,
  `glycogen degradation I`=-.128,`D-erythronate degradation II`=-.115,
  `NAD salvage III`=-.111,`taurine degradation`=-.101,`tRNA processing`=-.098,
  `2-carboxy-1,4-naphthoquinol`=-.096,`N-acetylneuraminate degr.`=-.088,
  `2,3-butanediol biosynthesis`=.083,`UMP biosynthesis II`=.117,
  `UMP biosynthesis I`=.125,`UMP biosynthesis III`=.128,
  `(R,R)-butanediol biosynthesis`=.143)
ff <- function(v,ttl){ o <- order(v); v <- v[o]
  plot(range(v)*1.25,c(1,length(v)),type="n",yaxt="n",ylab="",
       xlab="pooled Spearman rho with age",main=ttl)
  abline(v=0,col="grey70"); segments(0,seq_along(v),v,seq_along(v),col="grey60")
  points(v,seq_along(v),pch=16,col=ifelse(v<0,PAL[1],PAL[2]),cex=1.1)
  axis(2,at=seq_along(v),labels=names(v),las=1,cex.axis=.7,tick=FALSE) }
ff(gen,"a  Genera (q < 0.05)"); ff(pw,"b  Inferred pathways (q < 0.05)")
dev.off()

## ---- Figure 5 : disease ----
op("Fig5.pdf",7.2,3.6); par(mfrow=c(1,3),mar=c(4,4,2.5,1))
auc <- function(s,d){r<-rank(s);n1<-sum(d==1);n0<-sum(d==0)
  (sum(r[d==1])-n1*(n1+1)/2)/(n1*n0)}
sets <- list(OSCC=c("KUMPITSCH20_EU","LEE17_CN","TAIWAN_39064","WOLF17_EU"),
             periodontitis=c("PORTUGAL_74247","SPAIN_74796"))
plot(0:1,0:1,type="n",xlab="1 - specificity",ylab="sensitivity",main="a  Cross-cohort ROC")
abline(0,1,lty=3,col="grey"); i <- 0
for (gname in names(sets)) for (f in sets[[gname]]) {
  keep <- sm$cohort_id %in% sets[[gname]] & (sm$dis==0 |
    grepl(if(gname=="OSCC") "OSCC|carcinoma" else gname, sm$label, ignore.case=TRUE))
  te <- which(keep & sm$cohort_id==f); t2 <- which(keep & sm$cohort_id!=f)
  if (sum(sm$dis[te]==1)<5 || sum(sm$dis[te]==0)<5) next
  cv <- cv.glmnet(X[t2,],sm$dis[t2],family="binomial",alpha=.5,nfolds=5)
  s <- as.numeric(predict(cv,X[te,],s="lambda.min")); d <- sm$dis[te]
  o <- order(-s); i <- i+1
  lines(cumsum(d[o]==0)/sum(d==0), cumsum(d[o]==1)/sum(d==1),
        col=PAL[(i-1)%%6+1], lwd=1.8) }
legend("bottomright",c(sets$OSCC,sets$periodontitis),col=PAL[1:6],lwd=1.8,bty="n",cex=.6)
cn <- sub("^\\S+ ","",rownames(RV$disease))
D <- RV$disease; mb <- D[,"AUC"]; agx <- D[,"AUC_age"]
barplot(rbind(agx,mb),beside=TRUE,names.arg=rep("",6),col=PAL[c(6,1)],border=NA,
        ylim=c(0,1),ylab="AUC",main="b  Age alone vs microbiome"); abline(h=.5,lty=3)
text(seq(2,by=3,length.out=6),-.04,cn,srt=45,adj=1,xpd=TRUE,cex=.6)
legend("topleft",c("age only","microbiome"),fill=PAL[c(6,1)],border=NA,bty="n",cex=.75)
rar <- D[,"AUC_rarefied"]
plot(mb,rar,pch=16,col=PAL[1],cex=1.3,xlim=c(.6,.95),ylim=c(.6,.95),
     xlab="AUC, full depth",ylab="AUC, rarefied to 10,000",main="c  Depth control")
abline(0,1,lty=2); dev.off()

## ---- Figure 6 : oral-health gradient ----
op("Fig6.pdf",7.2,3.6); par(mfrow=c(1,2),mar=c(4,4,2.5,1))
ib <- which(sm$cohort_id %in% c("PORTUGAL_74247","SPAIN_74796") &
      (sm$dis==0 | grepl("periodontitis",sm$label,ignore.case=TRUE)))
cv <- cv.glmnet(X[ib,],sm$dis[ib],family="binomial",alpha=.5,nfolds=5)
rv <- which(sm$cohort_id=="RELVAS21_EU")
s  <- as.numeric(predict(cv,X[rv,],s="lambda.min"))
lab <- ifelse(sm$class[rv]=="train","grade 0",sub(".*Grade_","grade ",sm$label[rv]))
lab <- factor(lab,levels=c("grade 0","grade 1","grade 2","grade 3"))
boxplot(s~lab,col=adjustcolor(PAL[1],.4),border=PAL[1],ylab="oral-health score",
        xlab="",las=1,main="a  Unseen cohort"); 
mtext(sprintf("rho = %.3f", cor(s,as.integer(lab),method="spearman")),3,-1.2,cex=.75)
text(1:4,par("usr")[3]+.3,sprintf("%.0f y",tapply(sm$age[rv],lab,mean)),cex=.75,col="grey30")
hs <- as.numeric(predict(cv,X[tr,],s="lambda.min"))
loc <- s[lab=="grade 0"]
fl <- function(ref) sapply(levels(lab),function(g)
  100*mean(s[lab==g] >= quantile(ref,.90)))
barplot(rbind(fl(hs),fl(loc)),beside=TRUE,col=PAL[c(6,2)],border=NA,
        ylab="% flagged elevated",las=1,main="b  Reference anchoring",
        names.arg=levels(lab),cex.names=.8)
legend("topleft",c("published reference","same-laboratory controls"),
       fill=PAL[c(6,2)],border=NA,bty="n",cex=.75); dev.off()

cat("written to figures/:\n"); print(list.files("figures"))
