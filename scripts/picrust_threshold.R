# --------------------------------------------------------------------------
# picrust_threshold.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
md <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
KEEPS <- md$run_accession
tx <- readRDS("taxonomy/taxa_oriented.rds")$tax
bact <- !is.na(tx[,"Phylum"]) & !is.na(tx[,"Kingdom"]) & tx[,"Kingdom"]=="Bacteria"
names(bact) <- rownames(tx)
f <- list.files(c("final_split/train","final_split/valid"), pattern="rds$", full.names=TRUE)

glob <- list(); perc <- list(); tot <- list(); rd <- list()
for (p in f) {
  st <- readRDS(p); st <- st[rownames(st) %in% KEEPS, , drop=FALSE]
  if (!nrow(st)) next
  co <- md$cohort_id[match(rownames(st), md$run_accession)]
  glob[[p]] <- colSums(st > 0)
  perc[[p]] <- do.call(pmax, c(lapply(split(seq_len(nrow(st)), co),
                 function(i) colSums(st[i, , drop=FALSE] > 0)), list(0)))
  tot[[p]] <- colSums(st); rd[[p]] <- list(st=st, co=co)
}
agg <- function(L) { n <- unlist(lapply(L, names)); v <- unlist(L, use.names=FALSE)
  tapply(v, n, function(z) if (length(z)) max(z) else 0) }
sumagg <- function(L) { n <- unlist(lapply(L, names)); v <- unlist(L, use.names=FALSE)
  tapply(v, n, sum) }
G <- sumagg(glob); P <- agg(perc); b <- bact[names(G)]

retention <- function(sel) {
  k <- c(); cc <- c()
  for (p in names(rd)) { st <- rd[[p]]$st
    cm <- intersect(colnames(st), sel)
    k <- c(k, rowSums(st[, cm, drop=FALSE]) / rowSums(st)); cc <- c(cc, rd[[p]]$co) }
  q <- tapply(100*k, cc, median)
  c(nasv=length(sel), med=median(100*k), min=min(q), worst=names(which.min(q)))
}
cat(sprintf("%-34s %8s %9s %9s %s\n","rule","ASVs","median %","worst %","worst cohort"))
for (nm in c("global >=2","global >=5","within-cohort >=2","within-cohort >=3")) {
  sel <- switch(nm,
    "global >=2"        = names(G)[b & G >= 2],
    "global >=5"        = names(G)[b & G >= 5],
    "within-cohort >=2" = names(P)[b & P >= 2],
    "within-cohort >=3" = names(P)[b & P >= 3])
  r <- retention(sel)
  cat(sprintf("%-34s %8s %8.1f%% %8.1f%% %s\n", nm, r[1], as.numeric(r[2]),
      as.numeric(r[3]), r[4]))
}
