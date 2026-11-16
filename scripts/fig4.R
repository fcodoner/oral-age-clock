# --------------------------------------------------------------------------
# fig4.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
M <- readRDS("results/meta_detail.rds")
nms <- names(M)
TT <- c(genus="a  Genera", pathway="b  Inferred pathways")
pdf("figures/Fig4.pdf", width=13.6, height=5.8)
layout(matrix(seq_along(nms),1,length(nms)), widths=c(1,1.42)[seq_along(nms)])
for (k in seq_along(nms)) {
  nm <- nms[k]; d <- M[[nm]]; d <- d[d$q < 0.05, ]; d <- d[order(d$rho), ]
  lab <- substr(d$feature, 1, 52)
  lm <- min(max(nchar(lab))*0.40 + 1.0, 23)
  par(mar=c(3.8, lm, 2.4, 1.8), xpd=NA)
  xr <- range(c(d$lo, d$hi, 0)); pad <- diff(xr)*0.06
  plot(NA, xlim=c(xr[1]-pad, xr[2]+pad), ylim=c(0.4, nrow(d)+0.6),
       axes=FALSE, xlab="", ylab="")
  par(xpd=FALSE); abline(v=0, col="grey60", lwd=1); par(xpd=NA)
  cols <- ifelse(d$rho > 0, "#b1483f", "#2c6fa8")
  segments(d$lo, seq_len(nrow(d)), d$hi, seq_len(nrow(d)), col=cols, lwd=1.9)
  points(d$rho, seq_len(nrow(d)), pch=19, cex=1.05, col=cols)
  axis(1, cex.axis=0.84, mgp=c(2,0.55,0))
  mtext(lab, side=2, at=seq_len(nrow(d)), las=1, line=0.4, cex=0.74, adj=1)
  mtext("pooled Spearman rho with age", side=1, line=2.2, cex=0.84)
  box(bty="l", col="grey70")
  ttl <- if (nm %in% names(TT)) TT[[nm]] else paste0(letters[k], "  ", nm)
  mtext(ttl, side=3, line=0.7, adj=0, at=grconvertX(0.012,"nfc","user"), font=2, cex=1.02)
  cat(sprintf("%s: %d at q<0.05 | rho %.3f to %.3f | widest CI %.3f\n",
      nm, nrow(d), min(d$rho), max(d$rho), max(d$hi-d$lo)))
}
dev.off(); cat("figures/Fig4.pdf written\n")
