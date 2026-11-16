# --------------------------------------------------------------------------
# figure4.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
PAL <- c("#2E74B5","#C62828")
gen <- c(Alloprevotella=-.158, Haemophilus=-.157, Sneathia=-.134, Lautropia=-.129,
  `Absconditabacteria SR1`=-.109, Bergeyella=-.090, `Saccharibacteria TM7`=-.071,
  Lactobacillus=.080, Desulfovibrio=.085, `Bacteroidaceae [G-1]`=.086,
  `Veillonellaceae [G-1]`=.088, Anaeroglobus=.098, Bifidobacterium=.119)
pw <- c(`(Kdo)2-lipid A biosynthesis`=-.138, `glycogen degradation I`=-.128,
  `heme b biosynthesis (uroporphyrinogen)`=-.124, `heme b biosynthesis II`=-.123,
  `heme b biosynthesis V`=-.120, `D-erythronate degradation II`=-.115,
  `NAD salvage III`=-.111, `taurine degradation`=-.101, `tRNA processing`=-.098,
  `2-carboxy-1,4-naphthoquinol biosynthesis`=-.096,
  `N-acetylneuraminate degradation`=-.088, `2,3-butanediol biosynthesis`=.083,
  `UMP biosynthesis II`=.117, `UMP biosynthesis I`=.125, `UMP biosynthesis III`=.128,
  `(R,R)-butanediol biosynthesis`=.143)

pdf("figures/Fig4.pdf", width=7.2, height=4.8, pointsize=9); par(family="sans")
layout(matrix(1:2, 1, 2), widths=c(1, 1.22))
ff <- function(v, ttl, lmar, cx) {
  par(mar=c(4.2, lmar, 2.6, 1.6))
  o <- order(v); v <- v[o]
  plot(range(v)*1.3, c(0.6, length(v)+0.4), type="n", yaxt="n", ylab="", main="",
       xlab="pooled Spearman rho with age")
  title(main=ttl, adj=0, line=1.1, cex.main=1.05)
  abline(v=0, col="grey75")
  segments(0, seq_along(v), v, seq_along(v), col="grey65")
  points(v, seq_along(v), pch=16, col=ifelse(v < 0, PAL[1], PAL[2]), cex=1.15)
  axis(2, at=seq_along(v), labels=names(v), las=1, cex.axis=cx, tick=FALSE, line=-0.4)
}
ff(gen, "a  Genera",           8.4, 0.80)
ff(pw,  "b  Inferred pathways", 11.6, 0.68)
par(xpd=NA)
legend("bottomright", inset=c(0, -0.22), horiz=TRUE, bty="n", cex=0.8,
       legend=c("declines with age", "increases with age"), pch=16, col=PAL[1:2])
dev.off()
cat("figures/Fig4.pdf rewritten\n")
