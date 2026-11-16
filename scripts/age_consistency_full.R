# --------------------------------------------------------------------------
# age_consistency_full.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
set.seed(42)
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
G  <- readRDS("features/genus_counts.rds"); rownames(G) <- unname(rownames(G))
pf <- if (file.exists("picrust/pathways/path_abun_descrip.tsv.gz"))
        "picrust/pathways/path_abun_descrip.tsv.gz" else "picrust/pathways/path_abun_unstrat.tsv.gz"
P  <- read.delim(gzfile(pf), row.names=1, check.names=FALSE, quote="", comment.char="")
desc <- if (grepl("descrip", pf)) setNames(P[[1]], rownames(P)) else setNames(rownames(P), rownames(P))
if (grepl("descrip", pf)) P <- P[,-1]
P <- t(as.matrix(P))
ids <- intersect(rownames(P), rownames(G))
P <- P[ids,]; G <- G[ids,]; sm <- sm[match(ids, sm$run_accession), ]
k <- sm$class=="train"; P <- P[k,]; G <- G[k,]; sm <- sm[k,]
CONT <- c("Variovorax","Acinetobacter","Stenotrophomonas","Pseudomonas","Ralstonia",
          "Sphingomonas","Bradyrhizobium","Methylobacterium","Delftia","Comamonas",
          "Burkholderia","Cutibacterium","Escherichia")
rel <- function(M) { M <- as.matrix(M); storage.mode(M) <- "double"
  M[!is.finite(M)] <- 0; r <- rowSums(M); r[r==0] <- 1; M/r }
FS <- list(genus   = rel(G)[, setdiff(colnames(G)[colMeans(G>0)>0.10], CONT), drop=FALSE],
           pathway = rel(P)[, colMeans(P>0) > 0.10, drop=FALSE])
y <- sm$age; co <- sm$cohort_id
ag <- tapply(y, co, function(v) diff(range(v))); nn <- table(co)
use <- names(ag)[ag >= 40 & nn[names(ag)] >= 40]
cat(sprintf("healthy samples %d | informative cohorts %d\n", nrow(G), length(use)))

# DerSimonian-Laird random-effects meta-analysis on Fisher-z
meta <- function(rho, n) {
  ok <- !is.na(rho) & n > 4; rho <- rho[ok]; n <- n[ok]
  if (length(rho) < 5) return(c(NA,NA,NA,NA,NA))
  z <- atanh(pmax(pmin(rho, .999), -.999)); v <- 1/(n-3); w <- 1/v
  zf <- sum(w*z)/sum(w); Q <- sum(w*(z-zf)^2)
  tau2 <- max(0, (Q-(length(z)-1))/(sum(w)-sum(w^2)/sum(w)))
  ws <- 1/(v+tau2); zp <- sum(ws*z)/sum(ws); se <- sqrt(1/sum(ws))
  c(tanh(zp), 2*pnorm(-abs(zp/se)), Q/max(Q,1e-9), tanh(zp-1.96*se), tanh(zp+1.96*se))
}
META <- list()
for (nm in names(FS)) {
  X <- FS[[nm]]
  rho <- sapply(use, function(c) { i <- co==c
    apply(X[i,,drop=FALSE], 2, function(v) {
      ok <- is.finite(v)
      if (sum(ok) < 10 || sd(v[ok]) == 0) NA_real_
      else suppressWarnings(cor(v[ok], y[i][ok], method="spearman")) }) })
  ns  <- as.integer(nn[use])
  res <- t(apply(rho, 1, function(r) meta(r, ns)))
  colnames(res) <- c("rho","p","het","lo","hi")
  pos <- rowSums(rho > 0, na.rm=TRUE); kk <- rowSums(!is.na(rho))
  d <- data.frame(feature=rownames(rho), rho=round(res[,"rho"],3),
                  lo=round(res[,"lo"],3), hi=round(res[,"hi"],3),
                  agree=paste0(pmax(pos,kk-pos),"/",kk),
                  q=signif(p.adjust(res[,"p"], "BH"), 3), stringsAsFactors=FALSE)
  ## Store the full MetaCyc descriptor and keep a truncated copy for figure
  ## labels only; truncating the stored name made distinct pathways collide.
  if (nm=="pathway") { ff <- unname(desc[d$feature])
                       d$feature_full <- ff
                       d$feature      <- substr(ff, 1, 58) }
  d <- d[!is.na(d$q), ]; d <- d[order(d$q, -abs(d$rho)), ]
  META[[nm]] <- d
  cat(sprintf("\n===== %s: %d features tested, %d significant at q<0.05 =====\n",
      nm, nrow(d), sum(d$q < 0.05)))
  print(head(d[d$q < 0.05, ], 18), row.names=FALSE)
  if (sum(d$q<0.05) > 18) cat(sprintf("... and %d more\n", sum(d$q<0.05)-18))
}
saveRDS(META, "results/meta_detail.rds")
cat("\nsaved results/meta_detail.rds\n")
