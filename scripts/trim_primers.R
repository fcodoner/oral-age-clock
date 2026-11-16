#!/usr/bin/env Rscript
# --------------------------------------------------------------------------
# trim_primers.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# trim_primers.R v3 - normalise ASVs to common 16S coordinates.
# v3: sequences whose primer is unrecognised take the group's modal cut, so a
# cohort is never left internally ragged. Groups below MINRATE are left alone.
args <- commandArgs(trailingOnly = TRUE)
dry   <- "--dry-run" %in% args
indir <- Sys.getenv("TRIM_IN","processed"); outdir <- Sys.getenv("TRIM_OUT","processed/trimmed")
p <- read.delim(Sys.getenv("TRIM_PARAMS","metadata/dada2/dada2_params.tsv"), stringsAsFactors = FALSE)
if (!dry) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

FWD <- list("V1-V2"="AGAGTTTGATCMTGGCTCAG", "V1-V3"="AGAGTTTGATCMTGGCTCAG",
            "V3-V4"="CCTACGGGNGGCWGCAG",    "V4"="GTGYCAGCMGCCGCGGTAA")
REV <- list("V1-V2"=c("TGCTGCCTCCCGTAGGAGT"),
            "V1-V3"=c("ATTACCGCGGCTGCTGG"),
            "V3-V4"=c("GACTACHVGGGTATCTAATCC","GGACTACNVGGGTWTCTAAT"),
            "V4"   =c("GGACTACNVGGGTWTCTAAT"))
MAXOFF <- 25L; MINFRAG <- 10L; MINRATE <- 0.30

IU <- c(A="A",C="C",G="G",T="T",R="[AG]",Y="[CT]",S="[GC]",W="[AT]",K="[GT]",
        M="[AC]",B="[CGT]",D="[AGT]",H="[ACT]",V="[ACG]",N="[ACGT]")
re <- function(s) paste0(IU[strsplit(toupper(s),"")[[1]]], collapse="")
rc <- function(s) paste(rev(strsplit(chartr("ACGTRYSWKMBDHVN","TGCAYRSWMKVHDBN",
                       toupper(s)),"")[[1]]), collapse="")
wq <- function(L,w,q){o<-order(L);L<-L[o];w<-w[o];L[which(cumsum(w)/sum(w)>=q)[1]]}

cut5 <- function(s, prim) {
  L <- nchar(prim); cut <- integer(length(s)); todo <- rep(TRUE, length(s))
  for (o in 0:MAXOFF) { if (!any(todo)) break
    h <- todo & grepl(paste0("^.{",o,"}",re(prim)), s); cut[h] <- o+L; todo[h] <- FALSE }
  for (k in 1:(L-MINFRAG)) { if (!any(todo)) break
    h <- todo & grepl(paste0("^", re(substring(prim,k+1,L))), s)
    cut[h] <- L-k; todo[h] <- FALSE }
  cut
}
cut3 <- function(s, rp) {
  L <- nchar(rp); cut <- integer(length(s)); todo <- rep(TRUE, length(s))
  for (o in 0:MAXOFF) { if (!any(todo)) break
    h <- todo & grepl(paste0(re(rp),".{",o,"}$"), s); cut[h] <- o+L; todo[h] <- FALSE }
  for (k in 1:(L-MINFRAG)) { if (!any(todo)) break
    h <- todo & grepl(paste0(re(substring(rp,1,L-k)),"$"), s)
    cut[h] <- L-k; todo[h] <- FALSE }
  cut
}
fill <- function(cut, w) {                 # modal cut for unrecognised sequences
  h <- cut > 0
  if (!any(h) || sum(w[h])/sum(w) < MINRATE) return(list(cut=integer(length(cut)), mode=0L))
  tb <- tapply(w[h], cut[h], sum)
  md <- as.integer(names(tb)[which.max(tb)])
  cut[!h] <- md
  list(cut=cut, mode=md)
}

cat(sprintf("%-34s %-6s %7s %8s %8s %6s %6s %8s %8s %5s\n",
    "group","region","n_ASV","fwd.hit","rev.hit","f.mode","r.mode","len.in","len.out","drop"))
for (i in seq_len(nrow(p))) {
  g <- p$group[i]; f <- file.path(indir, paste0(g,"_seqtab.rds"))
  if (!file.exists(f)) next
  reg <- p$region[i]
  if (is.null(FWD[[reg]])) { cat(sprintf("%-34s %-6s  no primer set\n",g,reg)); next }
  s <- readRDS(f); q <- colnames(s); w <- colSums(s); tot <- sum(w)
  r5 <- cut5(q, FWD[[reg]]); h5 <- 100*sum(w[r5>0])/tot
  r3 <- NULL
  for (rp in REV[[reg]]) { cc <- cut3(q, rc(rp))
    if (is.null(r3) || sum(w[cc>0]) > sum(w[r3>0])) r3 <- cc }
  h3 <- 100*sum(w[r3>0])/tot
  f5 <- fill(r5, w); f3 <- fill(r3, w)
  st <- f5$cut + 1L; en <- nchar(q) - f3$cut
  keep <- (en - st + 1L) >= 50L
  new <- substr(q, st, en)
  cat(sprintf("%-34s %-6s %7d %7.1f%% %7.1f%% %6d %6d %8d %8d %5d\n", g, reg, ncol(s),
      h5, h3, f5$mode, f3$mode, wq(nchar(q),w,.5), wq(nchar(new[keep]),w[keep],.5), sum(!keep)))
  if (!dry) {
    s <- s[, keep, drop=FALSE]; new <- new[keep]
    if (anyDuplicated(new)) {
      ix <- split(seq_along(new), new)
      m <- vapply(ix, function(j) rowSums(s[, j, drop=FALSE]), numeric(nrow(s)))
      if (is.null(dim(m))) m <- matrix(m, nrow=nrow(s))
      dimnames(m) <- list(rownames(s), names(ix)); s <- m
    } else colnames(s) <- new
    saveRDS(s, file.path(outdir, paste0(g,"_seqtab.rds")))
  }
}
cat(if (dry) "\nDRY RUN - nothing written.\n" else sprintf("\nWritten to %s/\n", outdir))
