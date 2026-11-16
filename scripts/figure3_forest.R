# --------------------------------------------------------------------------
# figure3_forest.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
D <- readRDS("results/fig3_strata.rds")
cat("columns:", paste(names(D), collapse=", "), "\n"); print(utils::head(D,3))

nm <- names(D)
gc <- nm[nm %in% c("genus","genus_r","r_genus")][1]
ac <- nm[nm %in% c("asv","asv_r","r_asv")][1]
cc <- nm[grepl("cohort|group|fold", nm, ignore.case=TRUE)][1]
sc <- nm[grepl("strat|region|route", nm, ignore.case=TRUE)][1]
stopifnot(!is.na(gc), !is.na(ac), !is.na(cc))
D$.g <- D[[gc]]; D$.a <- D[[ac]]; D$.d <- D$.g - D$.a
D$.c <- as.character(D[[cc]])
D$.s <- if (is.na(sc)) "all cohorts" else as.character(D[[sc]])

w  <- wilcox.test(D$.g, D$.a, paired=TRUE, conf.int=TRUE)
md <- median(D$.d); ci <- as.numeric(w$conf.int)
nhi <- sum(D$.d > 0)

ord <- names(sort(tapply(D$.s, D$.s, length), decreasing=TRUE))
D$.s <- factor(D$.s, levels=ord)
D <- D[order(D$.s, -D$.d), ]; D$.s <- as.character(D$.s)
rows <- list(); y <- 0
for (s in unique(D$.s)) {
  y <- y - 1; rows[[length(rows)+1]] <- list(y=y, type="hdr", lab=s)
  for (i in which(D$.s == s)) { y <- y - 1
    rows[[length(rows)+1]] <- list(y=y, type="pt", lab=D$.c[i],
      d=D$.d[i], g=D$.g[i], a=D$.a[i]) }
}
y <- y - 1.4; rows[[length(rows)+1]] <- list(y=y, type="pool")
ylo <- y - 0.5

pdf("figures/Fig3.pdf", width=7.4, height=0.30*length(rows)+1.85)
par(mar=c(3.4,0.4,2.2,0.4), xpd=NA)
xr <- range(c(D$.d, ci, 0)); pad <- diff(xr)*0.06
xlo <- xr[1]-pad; xhi <- xr[2]+pad
Lx <- xlo - diff(xr)*0.42;  Rx <- xhi + diff(xr)*0.10
plot(NA, xlim=c(Lx, Rx + diff(xr)*0.34), ylim=c(ylo, 0.55),
     axes=FALSE, xlab="", ylab="")
segments(0, ylo+0.25, 0, 0.18, col="grey55", lwd=1)

for (r in rows) {
  if (r$type=="hdr") {
    text(Lx, r$y, r$lab, adj=0, font=3, cex=0.74, col="grey25")
  } else if (r$type=="pt") {
    text(Lx+diff(xr)*0.04, r$y, r$lab, adj=0, cex=0.72)
    text(Rx, r$y, sprintf("%+.2f / %+.2f", r$g, r$a), adj=0, cex=0.68, col="grey30")
    points(r$d, r$y, pch=19, cex=0.95,
           col=if (r$d>0) "#2c6fa8" else "#b1483f")
  } else {
    segments(ci[1], r$y, ci[2], r$y, lwd=2.4, col="grey20")
    points(md, r$y, pch=18, cex=1.9, col="grey10")
    text(Lx, r$y, "pooled (16 cohorts)", adj=0, font=2, cex=0.76)
    text(Rx, r$y, sprintf("%+.2f [%+.2f, %+.2f]", md, ci[1], ci[2]),
         adj=0, cex=0.70, col="grey20")
  }
}
axis(1, at=pretty(c(xlo,xhi), 6), cex.axis=0.78, mgp=c(2,0.6,0))
mtext("difference in within-cohort correlation with age  (genus - sequence variant)",
      side=1, line=1.8, cex=0.86)
mtext("Taxonomic resolution and cross-study transfer", side=3, line=1.1, adj=0,
      at=Lx, font=2, cex=1.0)
mtext(sprintf("genus higher in %d of %d cohorts; paired Wilcoxon P = %.2f", nhi, nrow(D), w$p.value),
      side=3, line=0.1, adj=0, at=Lx, cex=0.78, col="grey30")
text(Rx, 0.28, "genus r / ASV r", adj=0, cex=0.68, font=2, col="grey30")
dev.off()
cat(sprintf("\nFig3.pdf written | %d of %d | median %+.3f [%+.3f, %+.3f] | P = %.4f\n",
            nhi, nrow(D), md, ci[1], ci[2], w$p.value))
