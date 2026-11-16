# --------------------------------------------------------------------------
# loco_asv.R
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
MIN <- 3000; PREV_ASV <- 0.05; PREV_GEN <- 0.10

tx <- readRDS("taxonomy/taxa_oriented.rds")$tax
bact <- !is.na(tx[,"Phylum"]) & !is.na(tx[,"Kingdom"]) & tx[,"Kingdom"]=="Bacteria"
lab <- tx[,"Genus"]
for (r in c("Family","Order","Class","Phylum")) {
  i <- is.na(lab) & !is.na(tx[,r]); lab[i] <- paste0(tx[i,r], " (unclassified)")
}
names(bact) <- names(lab) <- rownames(tx)

# V3-V4 true-paired groups only; Chaudhari is forward-only despite its filename
EXCL <- "CHAUDHARI20_IN__paired__V3-V4_seqtab.rds"
ff <- list.files("final_split/train", pattern="paired__V3-V4_seqtab.rds$", full.names=TRUE)
ff <- ff[basename(ff) != EXCL]
cat("groups:", length(ff), "\n")

# pass 1 -- ASV prevalence, no full matrix in memory
pres <- list(); ntot <- 0L
for (p in ff) {
  st <- readRDS(p); s <- colnames(st); st <- st[, bact[s], drop=FALSE]
  st <- st[rowSums(st) >= MIN, , drop=FALSE]
  ntot <- ntot + nrow(st); pres[[p]] <- colSums(st > 0)
}
v  <- unlist(pres, use.names=TRUE); names(v) <- sub("^.*\\.", "", names(v))
pv <- tapply(v, names(v), sum)
sel <- names(pv)[pv >= PREV_ASV * ntot]
cat(sprintf("samples %d | ASVs kept at prevalence>=%.0f%%: %d\n", ntot, 100*PREV_ASV, length(sel)))

# pass 2 -- build ASV and genus matrices on the same samples
A <- matrix(0L, 0, length(sel), dimnames=list(NULL, sel))
glist <- list(); meta <- list()
for (p in ff) {
  g  <- sub("_seqtab.rds$","",basename(p))
  st <- readRDS(p); s <- colnames(st); st <- st[, bact[s], drop=FALSE]
  st <- st[rowSums(st) >= MIN, , drop=FALSE]
  a  <- matrix(0L, nrow(st), length(sel), dimnames=list(rownames(st), sel))
  cm <- intersect(colnames(st), sel); a[, cm] <- st[, cm, drop=FALSE]
  A <- rbind(A, a)
  glist[[g]] <- t(rowsum(t(st), lab[colnames(st)]))
  meta[[g]] <- data.frame(run_accession=rownames(st), cohort_id=strsplit(g,"__")[[1]][1],
                          stringsAsFactors=FALSE)
}
gen <- sort(unique(unlist(lapply(glist, colnames))))
G <- matrix(0L, nrow(A), length(gen), dimnames=list(rownames(A), gen))
r0 <- 0L
for (cs in glist) { G[(r0+1):(r0+nrow(cs)), colnames(cs)] <- cs; r0 <- r0 + nrow(cs) }
G <- G[, colMeans(G>0) > PREV_GEN, drop=FALSE]

sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
md <- do.call(rbind, meta)
md$age <- sm$age[match(md$run_accession, sm$run_accession)]
stopifnot(!any(is.na(md$age)), nrow(md)==nrow(A))
y <- md$age; co <- md$cohort_id
cat(sprintf("matrices: ASV %d x %d | genus %d x %d\n", nrow(A), ncol(A), nrow(G), ncol(G)))

ag <- tapply(y, co, function(v) diff(range(v))); nn <- table(co)
folds <- names(ag)[ag >= 40 & nn[names(ag)] >= 40]
cat("folds:", paste(folds, collapse=", "), "\n\n")

clr <- function(M) { L <- log(M + 0.5); L - rowMeans(L) }
run <- function(X, nm) {
  P <- c(); Y <- c(); C <- c()
  for (f in folds) {
    te <- co==f; tr <- !te
    cv <- cv.glmnet(X[tr,], y[tr], alpha=0.5, nfolds=5, standardize=TRUE)
    P <- c(P, as.numeric(predict(cv, X[te,], s="lambda.min")))
    Y <- c(Y, y[te]); C <- c(C, rep(f, sum(te)))
  }
  wr <- sapply(unique(C), function(f){i<-C==f; suppressWarnings(cor(P[i],Y[i]))})
  nu <- sapply(folds, function(f) mean(abs(y[co==f] - mean(y[co!=f]))))
  cat(sprintf("%-8s MAE %6.2f | nullMAE %6.2f | r %6.3f | R2 %6.3f | medWithinR %6.3f\n",
      nm, mean(abs(P-Y)), weighted.mean(nu, table(C)[folds]), cor(P,Y),
      1-sum((P-Y)^2)/sum((Y-mean(Y))^2), median(wr, na.rm=TRUE)))
  print(round(sort(wr), 3)); cat("\n"); invisible(wr)
}
cat("=== genus level (same samples/folds) ===\n"); wg <- run(clr(G), "genus")
cat("=== ASV level ===\n");                        wa <- run(clr(A), "ASV")
cat("per-fold change in within-cohort r (ASV - genus):\n")
print(round(wa[names(wg)] - wg, 3))
