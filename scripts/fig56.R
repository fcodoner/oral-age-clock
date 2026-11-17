# --------------------------------------------------------------------------
# fig56.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
rv <- readRDS("results/reported_values.rds")
D  <- rv$disease
roc <- function(s,d){ o<-order(-s); d<-d[o]
  list(fpr=c(0,cumsum(1-d)/sum(1-d)), tpr=c(0,cumsum(d)/sum(d))) }
short <- function(x) sub("^(OSCC|periodontitis) ","",x)

## ---------------- Figure 5 ----------------
RC <- readRDS("results/disease_roc.rds")
cols <- c("#2c6fa8","#b1483f","#3f8f5a","#d08a1e","#6f5ba8","#7a7a7a")
pdf("figures/Fig5.pdf", width=12.6, height=4.5)
layout(matrix(1:3,1,3), widths=c(1.05,1.15,0.95))

par(mar=c(4.2,4.4,3.0,0.8), xpd=NA)
plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="1 - specificity", ylab="sensitivity",
     cex.lab=0.98, cex.axis=0.88)
par(xpd=FALSE); abline(0,1,lty=3,col="grey65"); par(xpd=NA)
for (i in seq_along(RC)) { r <- roc(RC[[i]]$score, RC[[i]]$label)
  lines(r$fpr, r$tpr, col=cols[i], lwd=2.1) }
mtext("a  Cross-cohort ROC", side=3, line=1.2, adj=0, at=-0.17, font=2, cex=1.02)
legend("bottomright", bty="n", cex=0.70, lwd=2.1, col=cols[seq_along(RC)],
  legend=sprintf("%-16s %.2f", short(names(RC)), D[names(RC),"AUC"]))

par(mar=c(8.6,4.4,3.0,0.8))
bp <- barplot(rbind(D[,"AUC_age"], D[,"AUC"]), beside=TRUE, border=NA,
  col=c("#9a9a9a","#2c6fa8"), ylim=c(0,1), ylab="AUC", las=2,
  names.arg=short(rownames(D)), cex.names=0.80, cex.axis=0.88)
par(xpd=FALSE); abline(h=0.5, lty=3, col="grey35"); par(xpd=NA)
mtext("b  Age alone vs microbiome", side=3, line=1.2, adj=0, at=bp[1,1]-1.2, font=2, cex=1.02)
legend("topleft", c("age only","microbiome"), fill=c("#9a9a9a","#2c6fa8"),
  border=NA, bty="n", cex=0.80)

par(mar=c(4.2,4.4,3.0,1.2))
plot(D[,"AUC"], D[,"AUC_rarefied"], pch=19, cex=1.25, col="#2c6fa8",
  xlim=c(0.6,0.95), ylim=c(0.6,0.95), xlab="AUC, full depth",
  ylab="AUC, rarefied to 10,000", cex.lab=0.98, cex.axis=0.88)
par(xpd=FALSE); abline(0,1,lty=2,col="grey25"); par(xpd=NA)
mtext("c  Depth control", side=3, line=1.2, adj=0, at=0.555, font=2, cex=1.02)
dev.off()
cat("figures/Fig5.pdf written\n")

## ---------------- Figure 6 ----------------
GD <- readRDS("results/gradient_detail.rds")
g  <- rv$gradient; AN <- rv$anchoring
SM <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
av <- SM$age[SM$cohort_id=="RELVAS21_EU"]
ma <- if (length(av)==length(GD$score)) round(tapply(av, GD$grade, mean, na.rm=TRUE)) else NULL

pdf("figures/Fig6.pdf", width=11.0, height=5.0)
layout(matrix(1:2,1,2), widths=c(1,1.08))

par(mar=c(5.4,4.4,3.0,1.0), xpd=NA)
lv <- levels(GD$grade)
boxplot(GD$score ~ GD$grade, col="#cfe0ef", border="#2c6fa8", lwd=1.5,
  xaxt="n", xlab="", ylab="oral-health score", cex.axis=0.9, cex.lab=0.98)
axis(1, at=seq_along(lv), labels=FALSE)
mtext(lv, side=1, line=0.9, at=seq_along(lv), cex=0.88)
if (!is.null(ma)) mtext(sprintf("(%d y)", ma), side=1, line=1.9, at=seq_along(lv), cex=0.80, col="grey35")
mtext("clinician-assigned grade (mean age)", side=1, line=3.4, cex=0.92)
mtext("a  Unseen cohort", side=3, line=1.2, adj=0, at=0.25, font=2, cex=1.04)
u <- par("usr")
text(u[2], u[3]+diff(u[3:4])*0.06,
  sprintf("rho = %.2f\nP = %.1e", g[["rho.rho"]], g[["p"]]),
  adj=c(1.1,0), cex=0.82, col="grey20")

par(mar=c(4.4,4.4,3.0,1.0))
## Bootstrap the 17 local controls to put limits on each local flag rate.
set.seed(3)
bq <- replicate(2000, quantile(sample(GD$local_ctrl, length(GD$local_ctrl), replace=TRUE), 0.90))
ci <- sapply(seq_along(lv), function(j) {
        fr <- sapply(bq, function(th) 100*mean(GD$score[GD$grade == lv[j]] >= th))
        quantile(fr, c(0.025, 0.975)) })
yhi <- max(c(AN, ci))
bp <- barplot(AN, beside=TRUE, border=NA, col=c("#9a9a9a","#b1483f"),
  ylab="% flagged elevated", ylim=c(0, yhi*1.18), cex.axis=0.9, cex.names=0.92)
arrows(bp[2,], ci[1,], bp[2,], ci[2,], angle=90, code=3, length=0.035, lwd=1.3, col="#8d3830")
for (j in 1:ncol(AN)) for (i in 1:2) {
  yv <- if (i == 2) max(AN[i,j], ci[2,j]) else AN[i,j]
  text(bp[i,j], yv + yhi*0.04, sprintf("%.0f", AN[i,j]),
       cex=0.76, col=c("grey35","#b1483f")[i]) }
mtext("b  Reference anchoring", side=3, line=1.2, adj=0, at=bp[1,1]-1.0, font=2, cex=1.04)
legend("topleft", c("reference from training samples","controls from the same study"),
  fill=c("#9a9a9a","#b1483f"), border=NA, bty="n", cex=0.80)
## The local threshold is the 90th percentile of the same controls it is then applied to,
## so the grade 0 bar is circular by construction. Overlay the leave-one-out estimate.
g0  <- which(GD$grade == lv[1])
loo <- 100*mean(sapply(g0, function(k) GD$score[k] >= quantile(GD$score[setdiff(g0, k)], 0.90)))
segments(bp[2,1], AN[2,1], bp[2,1], loo, col="#b1483f", lwd=1.5, lty=3)
points(bp[2,1], loo, pch=18, cex=1.6, col="#b1483f")
text(bp[2,1], loo + yhi*0.045, sprintf("%.0f", loo), cex=0.74, col="#b1483f", font=3)
legend("topright", legend="grade 0, leave-one-out threshold", pch=18, col="#b1483f",
  bty="n", cex=0.74, pt.cex=1.4)
cat(sprintf("grade 0 flagged: %.1f%% in-sample, %.1f%% leave-one-out\n", AN[2,1], loo))
dev.off()
cat("figures/Fig6.pdf written\n")
cat(sprintf("\nrho=%.3f p=%.2e | published %s | local %s\n", g[["rho.rho"]], g[["p"]],
  paste(round(AN[1,],1),collapse="/"), paste(round(AN[2,],1),collapse="/")))
