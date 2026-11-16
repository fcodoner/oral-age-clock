# --------------------------------------------------------------------------
# figure6.R
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
CT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
        "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
        "Burkholderia","Cutibacterium","Escherichia")
G  <- readRDS("features/genus_counts.rds"); rownames(G) <- unname(rownames(G))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(G), sm$run_accession), ]
sm$dis <- ifelse(sm$class=="train", 0L, 1L)
K <- setdiff(colnames(G)[colMeans(G>0)>0.10], CT)
L <- log(G[,K,drop=FALSE]+0.5); X <- L - rowMeans(L)
tr <- which(sm$class=="train")
PAL <- c("#2E74B5","#C62828","#2E7D32","#E8871A","#6A4C93","#777777")

ib <- which(sm$cohort_id %in% c("PORTUGAL_74247","SPAIN_74796") &
      (sm$dis==0 | grepl("periodontitis", sm$label, ignore.case=TRUE)))
# Fixed fold assignment: deterministic, reproducible on any machine.
fid <- rep(1:5, length.out=length(ib))
cv <- cv.glmnet(X[ib,], sm$dis[ib], family="binomial", alpha=.5, foldid=fid)
rv <- which(sm$cohort_id=="RELVAS21_EU")
s  <- as.numeric(predict(cv, X[rv,], s="lambda.min"))
lab <- ifelse(sm$class[rv]=="train","grade 0", sub(".*Grade_","grade ", sm$label[rv]))
lab <- factor(lab, levels=c("grade 0","grade 1","grade 2","grade 3"))
age <- tapply(sm$age[rv], lab, mean)
ct  <- suppressWarnings(cor.test(s, as.integer(lab), method="spearman"))

pdf("figures/Fig6.pdf", width=7.2, height=3.8, pointsize=9); par(family="sans")
par(mfrow=c(1,2), mar=c(5.2,4,2.5,1))
boxplot(s ~ lab, col=adjustcolor(PAL[1],.35), border=PAL[1], outcol=PAL[1],
        ylab="oral-health score", xlab="", las=1, xaxt="n", main="a  Unseen cohort")
axis(1, at=1:4, labels=FALSE)
mtext(levels(lab), side=1, at=1:4, line=0.9, cex=.9)
mtext(sprintf("(%.0f y)", age), side=1, at=1:4, line=1.9, cex=.78, col="grey35")
mtext("clinician-assigned grade (mean age)", side=1, line=3.4, cex=.8)
legend("bottomright", bty="n", cex=.8, inset=c(0,0.02),
       legend=c(sprintf("rho = %.2f", ct$estimate),
                sprintf("P = %.1e", ct$p.value)))

hs  <- as.numeric(predict(cv, X[tr,], s="lambda.min"))
loc <- s[lab=="grade 0"]
fl  <- function(ref) sapply(levels(lab), function(g) 100*mean(s[lab==g] >= quantile(ref,.90)))
par(mar=c(5.2,4,2.5,1))
bp <- barplot(rbind(fl(hs), fl(loc)), beside=TRUE, col=PAL[c(6,2)], border=NA,
        ylab="% flagged elevated", las=1, ylim=c(0,95),
        main="b  Reference anchoring", names.arg=levels(lab), cex.names=.85)
v <- rbind(fl(hs), fl(loc))
text(bp, v+3.5, sprintf("%.0f", v), cex=.7, col=PAL[c(6,2)])
legend("topleft", c("reference from published cohorts","controls from the same laboratory"),
       fill=PAL[c(6,2)], border=NA, bty="n", cex=.75)
dev.off()
cat("figures/Fig6.pdf rewritten\n")
