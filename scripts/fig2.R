# --------------------------------------------------------------------------
# fig2.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
M  <- readRDS("results/matched_comparison.rds")
rv <- readRDS("results/reported_values.rds")
mae <- function(p,a) mean(abs(p-a))
L <- mae(M$loco,M$y); R <- mae(M$rand,M$y)
rl <- cor(M$loco,M$y); rr <- cor(M$rand,M$y)
gp <- M$boot$gap
cat(sprintf("LOCO %.3f r %.3f | rand %.3f r %.3f | gap %.2f [%.2f, %.2f] | n=%d, %d cohorts\n",
    L, rl, R, rr, mean(gp), quantile(gp,.025), quantile(gp,.975), length(M$y), length(unique(M$cohort))))

## panel a values are read from the PERMANOVA result, never typed in by hand.
## Run scripts/permanova.R first; it writes results/permanova.rds.
PV <- readRDS("results/permanova.rds")
ti <- c("age","sex2","region","route","cohort_id")
tm <- c("age","sex","region","route","cohort")
pc <- function(A) 100 * A$R2[match(ti, rownames(A))]
bf <- pc(PV$study_first)   # age fitted last, after the study terms
af <- pc(PV$age_first)     # age fitted first
stot <- sum(bf[3:5])       # region + route + cohort, study terms fitted first
cat(sprintf("panel a: study terms first %s | age first %s | study total %.2f%%\n",
    paste(sprintf("%.2f", bf), collapse="/"),
    paste(sprintf("%.2f", af), collapse="/"), stot))

pdf("figures/Fig2.pdf", width=10.4, height=5.0)
layout(matrix(1:2,1,2), widths=c(1,1.05))
par(mar=c(4.0,4.2,3.0,1.0), xpd=NA)
bp <- barplot(rbind(bf,af), beside=TRUE, names.arg=tm, border=NA,
  col=c("#2c6fa8","#b1483f"), ylab="variance explained (%)",
  ylim=c(0, max(bf,af)*1.18), cex.names=0.95, cex.axis=0.9)
mtext("a  Variance partition", side=3, line=1.2, adj=0, at=bp[1,1]-0.9, font=2, cex=1.05)
legend("topleft", c("study terms first","age first"), fill=c("#2c6fa8","#b1483f"),
  border=NA, bty="n", cex=0.88, inset=c(0.02,0.02))
text(colMeans(bp)[5], max(bf)*1.07, sprintf("%.1f%% study terms", stot),
  cex=0.78, col="grey30")

par(mar=c(4.0,4.2,3.0,1.0))
rg <- range(c(M$y, M$loco, M$rand))
plot(M$y, M$rand, pch=19, cex=0.34, col=adjustcolor("#3f8f5a",0.40),
  xlim=rg, ylim=rg, xlab="chronological age (years)", ylab="predicted age (years)",
  cex.axis=0.9, cex.lab=0.98)
points(M$y, M$loco, pch=19, cex=0.34, col=adjustcolor("#b1483f",0.40))
par(xpd=FALSE); abline(0,1,lty=2,col="grey25"); par(xpd=NA)
lines(lowess(M$y, M$rand, f=0.6), col="#2f6b45", lwd=2.4)
lines(lowess(M$y, M$loco, f=0.6), col="#8d3830", lwd=2.4)
mtext("b  Evaluation scheme", side=3, line=1.2, adj=0, at=rg[1]-diff(rg)*0.13, font=2, cex=1.05)
leg  <- c(sprintf("random CV   MAE %.2f   r %.2f", R, rr),
          sprintf("LOCO             MAE %.2f   r %.2f", L, rl))
lcol <- c("#2f6b45","#8d3830")
if (file.exists("results/revision_round1.rds")) {
  RR   <- readRDS("results/revision_round1.rds")$matched
  leg  <- c(leg, sprintf("18-fold CV  MAE %.2f   r %.2f", mae(RR$r18, RR$y), cor(RR$r18, RR$y)))
  lcol <- c(lcol, "#5aa87a")
}
legend("topleft", bty="n", cex=0.80, text.col=lcol, legend=leg)
legend("bottomright", bty="n", cex=0.74, text.col="grey25",
  legend=c(sprintf("difference %.2f y  [%.2f, %.2f]", mean(gp), quantile(gp,.025), quantile(gp,.975)),
           sprintf("%d samples, %d cohorts", length(M$y), length(unique(M$cohort)))))

dev.off()
cat("figures/Fig2.pdf written\n")
