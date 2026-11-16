# --------------------------------------------------------------------------
# figure3.R
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
tx <- readRDS("taxonomy/taxa_oriented.rds")$tax
bb <- !is.na(tx[,"Phylum"]) & !is.na(tx[,"Kingdom"]) & tx[,"Kingdom"]=="Bacteria"
names(bb) <- rownames(tx)

EX  <- "CHAUDHARI20_IN__paired__V3-V4_seqtab.rds"
ffs <- list.files("final_split/train", pattern="paired__V3-V4_seqtab.rds$", full.names=TRUE)
ffs <- ffs[basename(ffs) != EX]
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
  md[[p]]<-data.frame(id=rownames(st), coh=strsplit(sub("_seqtab.rds$","",basename(p)),"__")[[1]][1]) }
mdd <- do.call(rbind, md)
y  <- sm$age[match(mdd$id, sm$run_accession)]; co <- mdd$coh
K  <- setdiff(colnames(G)[colMeans(G>0)>0.10], CT)
XG <- clr(G[match(mdd$id, rownames(G)), K, drop=FALSE])
XA <- clr(A)
ag <- tapply(y,co,function(z) diff(range(z))); nn <- table(co)
folds <- names(ag)[ag>=40 & nn[names(ag)]>=40]
cat(sprintf("samples %d | genera %d | ASVs %d | folds %d\n",
    nrow(A), ncol(XG), ncol(XA), length(folds)))

perfold <- function(X){ sapply(folds, function(f){ te <- co==f
  cv <- cv.glmnet(X[!te,,drop=FALSE], y[!te], alpha=.5, foldid=FF(sum(!te)), standardize=TRUE)
  p  <- as.numeric(predict(cv, X[te,,drop=FALSE], s="lambda.min"))
  suppressWarnings(cor(p, y[te])) }) }
rg <- perfold(XG); ra <- perfold(XA)
wt <- suppressWarnings(wilcox.test(rg, ra, paired=TRUE))
saveRDS(list(genus=rg, asv=ra, n=nrow(A), wilcox=wt), "results/fig3_perfold.rds")

cat("\n%-20s %8s %8s %9s\n")
cat(sprintf("%-20s %8s %8s %9s\n","cohort","genus r","ASV r","diff"))
for (f in folds) cat(sprintf("%-20s %8.3f %8.3f %9.3f\n", f, rg[f], ra[f], rg[f]-ra[f]))
cat(sprintf("\nmedian: genus %.3f | ASV %.3f | genus higher in %d of %d folds\n",
    median(rg), median(ra), sum(rg>ra), length(folds)))
cat(sprintf("paired Wilcoxon signed-rank: V = %s, P = %.4f\n", wt$statistic, wt$p.value))

PAL <- c("#2E74B5","#C62828")
pdf("figures/Fig3.pdf", width=4.8, height=4.6, pointsize=9); par(family="sans")
par(mar=c(4,4.2,2.6,7.5))
plot(c(0.75,2.25), range(c(rg,ra))*c(0.85,1.15), type="n", xaxt="n",
     xlab="", ylab="within-cohort correlation with age", main="")
title(main="Feature representation", adj=0, line=1.1)
axis(1, at=1:2, labels=c("ASV","genus"))
segments(1, ra, 2, rg, col=ifelse(rg>ra, adjustcolor(PAL[1],.65), adjustcolor(PAL[2],.65)), lwd=1.8)
points(rep(1,length(ra)), ra, pch=16, col=PAL[2], cex=1.1)
points(rep(2,length(rg)), rg, pch=16, col=PAL[1], cex=1.1)
segments(0.88,median(ra),1.12,median(ra), lwd=3, col=PAL[2])
segments(1.88,median(rg),2.12,median(rg), lwd=3, col=PAL[1])
text(2.3, rg, names(rg), adj=0, cex=0.62, xpd=NA, col="grey30")
mtext(sprintf("median %.2f", median(ra)), side=1, at=1, line=2.1, cex=.8, col=PAL[2])
mtext(sprintf("median %.2f", median(rg)), side=1, at=2, line=2.1, cex=.8, col=PAL[1])
legend("topleft", bty="n", cex=.78,
  legend=sprintf("genus higher in %d of %d cohorts (P = %.3f)", sum(rg>ra), length(folds), wt$p.value))
dev.off()

cat("\n---- paste into Results, section 'Taxonomic aggregation level determines transferability' ----\n")
cat(sprintf(
"Within a single region and processing route, where sequence variants are directly comparable across cohorts, we compared variant-level and genus-level features on identical samples and identical folds (%d samples, %d cohorts). Genus-level features gave a higher within-cohort correlation with age in %d of %d cohorts (median %.2f against %.2f; paired Wilcoxon signed-rank P = %.3f), with the two representations differing by up to %.2f in a single cohort (Fig. 3).\n",
nrow(A), length(folds), sum(rg>ra), length(folds), median(rg), median(ra), wt$p.value, max(rg-ra)))
