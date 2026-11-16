# --------------------------------------------------------------------------
# figure3_strata.R
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
FF <- function(n) rep(1:5, length.out=n)
CT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
        "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
        "Burkholderia","Cutibacterium","Escherichia")
G  <- readRDS("features/genus_counts.rds"); rownames(G) <- unname(rownames(G))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(G), sm$run_accession), ]
clr <- function(M){L <- log(M+0.5); L-rowMeans(L)}
K  <- setdiff(colnames(G)[colMeans(G>0)>0.10], CT)
tx <- readRDS("taxonomy/taxa_oriented.rds")$tax
bb <- !is.na(tx[,"Phylum"]) & !is.na(tx[,"Kingdom"]) & tx[,"Kingdom"]=="Bacteria"
names(bb) <- rownames(tx)
# four groups labelled paired are forward-only; treat them as single-end
FWD <- c("PORTUGAL_74247__paired__V1-V2","SPAIN_74796__paired__V1-V2",
         "CHAUDHARI20_IN__paired__V3-V4","ANTARCTIC_19274__paired__V4")

fs <- list.files("final_split/train", pattern="rds$", full.names=TRUE)
key <- sapply(fs, function(p){ b <- sub("_seqtab.rds$","",basename(p))
  pt <- strsplit(b,"__")[[1]]; rt <- if (b %in% FWD) "single" else pt[2]
  paste(pt[3], rt) })
res <- list()
for (st in unique(key)) {
  gs <- fs[key==st]; if (length(gs) < 2) next
  pr <- list(); nt <- 0L
  for (p in gs){ x<-readRDS(p); x<-x[,bb[colnames(x)],drop=FALSE]
    x<-x[rowSums(x)>=3000,,drop=FALSE]; if(!nrow(x)) next
    nt<-nt+nrow(x); pr[[p]]<-colSums(x>0) }
  if (nt < 150) next
  v <- unlist(pr,use.names=TRUE); names(v) <- sub("^.*\\.(?=[ACGT])","",names(v),perl=TRUE)
  pv <- tapply(v,names(v),sum); sel <- names(pv)[pv >= 0.05*nt]
  if (length(sel) < 50) next
  A <- NULL; md <- list()
  for (p in gs){ x<-readRDS(p); x<-x[,bb[colnames(x)],drop=FALSE]
    x<-x[rowSums(x)>=3000,,drop=FALSE]; if(!nrow(x)) next
    m<-matrix(0L,nrow(x),length(sel),dimnames=list(rownames(x),sel))
    cm<-intersect(colnames(x),sel); m[,cm]<-x[,cm,drop=FALSE]; A<-rbind(A,m)
    md[[p]]<-data.frame(id=rownames(x), coh=strsplit(sub("_seqtab.rds$","",basename(p)),"__")[[1]][1]) }
  mdd <- do.call(rbind, md)
  y <- sm$age[match(mdd$id, sm$run_accession)]; co <- mdd$coh
  ok <- !is.na(y); A<-A[ok,,drop=FALSE]; y<-y[ok]; co<-co[ok]; mdd<-mdd[ok,]
  ag <- tapply(y,co,function(z) diff(range(z))); nn <- table(co)
  fo <- names(ag)[ag>=40 & nn[names(ag)]>=40]
  if (!length(fo) || length(unique(co)) < 2) next
  XG <- clr(G[match(mdd$id, rownames(G)), K, drop=FALSE]); XA <- clr(A)
  pf <- function(X) sapply(fo, function(f){ te<-co==f
    if (sum(!te) < 40) return(NA_real_)
    cv<-cv.glmnet(X[!te,,drop=FALSE],y[!te],alpha=.5,foldid=FF(sum(!te)),standardize=TRUE)
    suppressWarnings(cor(as.numeric(predict(cv,X[te,,drop=FALSE],s="lambda.min")), y[te])) })
  rg<-pf(XG); ra<-pf(XA); keep<-!is.na(rg)&!is.na(ra)
  if (!any(keep)) next
  res[[st]] <- data.frame(stratum=st, cohort=fo[keep], genus=rg[keep], asv=ra[keep],
                          n=as.integer(nn[fo[keep]]), nASV=ncol(XA))
  cat(sprintf("%-18s %4d samples | %5d ASVs | %d folds\n", st, nrow(A), ncol(XA), sum(keep)))
}
D <- do.call(rbind, res); rownames(D) <- NULL
D$diff <- D$genus - D$asv
saveRDS(D, "results/fig3_strata.rds")

cat("\n"); print(format(D[,c("stratum","cohort","n","genus","asv","diff")], digits=3), row.names=FALSE)
wt <- suppressWarnings(wilcox.test(D$genus, D$asv, paired=TRUE))
bs <- replicate(5000, median(sample(D$diff, nrow(D), replace=TRUE)))
cat(sprintf("\npairs %d | genus higher in %d (%.0f%%) | median difference %.3f [95%% CI %.3f, %.3f]\n",
    nrow(D), sum(D$diff>0), 100*mean(D$diff>0), median(D$diff),
    quantile(bs,.025), quantile(bs,.975)))
cat(sprintf("paired Wilcoxon signed-rank P = %.4f\n", wt$p.value))
cat(sprintf("median within-cohort r: genus %.3f | ASV %.3f\n", median(D$genus), median(D$asv)))

PAL <- c("#2E74B5","#C62828")
pdf("figures/Fig3.pdf", width=5.4, height=4.8, pointsize=9); par(family="sans")
par(mar=c(4.2,4.2,2.8,8))
plot(c(0.75,2.25), range(c(D$genus,D$asv))*c(0.8,1.15), type="n", xaxt="n",
     xlab="", ylab="within-cohort correlation with age", main="")
title(main="Feature representation", adj=0, line=1.2)
axis(1, at=1:2, labels=c("sequence variant","genus"))
segments(1, D$asv, 2, D$genus, lwd=1.7,
  col=ifelse(D$diff>0, adjustcolor(PAL[1],.6), adjustcolor(PAL[2],.6)))
points(rep(1,nrow(D)), D$asv, pch=16, col=PAL[2], cex=1.05)
points(rep(2,nrow(D)), D$genus, pch=16, col=PAL[1], cex=1.05)
segments(.88,median(D$asv),1.12,median(D$asv),lwd=3,col=PAL[2])
segments(1.88,median(D$genus),2.12,median(D$genus),lwd=3,col=PAL[1])
text(2.3, D$genus, D$cohort, adj=0, cex=.55, xpd=NA, col="grey30")
mtext(sprintf("median %.2f", median(D$asv)), 1, at=1, line=2.2, cex=.8, col=PAL[2])
mtext(sprintf("median %.2f", median(D$genus)), 1, at=2, line=2.2, cex=.8, col=PAL[1])
legend("topleft", bty="n", cex=.76,
  legend=c(sprintf("genus higher in %d of %d cohorts", sum(D$diff>0), nrow(D)),
           sprintf("median difference %.2f (P = %.3f)", median(D$diff), wt$p.value)))
dev.off()
cat("\nfigures/Fig3.pdf written\n")
