# --------------------------------------------------------------------------
# prep_picrust.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
args  <- commandArgs(TRUE)
BUILD <- "--build" %in% args
MINS  <- if (any(grepl("^--minsamples=", args)))
            as.integer(sub("^--minsamples=","", grep("^--minsamples=", args, value=TRUE))) else 10L

tx   <- readRDS("taxonomy/taxa_oriented.rds")
seqs <- tx$oriented_seq; names(seqs) <- rownames(tx$tax)
bact <- !is.na(tx$tax[,"Phylum"]) & !is.na(tx$tax[,"Kingdom"]) & tx$tax[,"Kingdom"]=="Bacteria"
names(bact) <- rownames(tx$tax)

f <- list.files(c("final_split/train","final_split/valid"), pattern="rds$", full.names=TRUE)
KEEPS <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)$run_accession
nmv <- character(0); prs <- integer(0); tot <- numeric(0); nsamp <- 0L
for (p in f) {
  st <- readRDS(p); st <- st[rownames(st) %in% KEEPS, , drop=FALSE]
  if (!nrow(st)) next
  nsamp <- nsamp + nrow(st)
  nmv <- c(nmv, colnames(st))
  prs <- c(prs, as.integer(colSums(st > 0)))
  tot <- c(tot, as.numeric(colSums(st)))
}
pv <- tapply(prs, nmv, sum); tv <- tapply(tot, nmv, sum)
b  <- bact[names(pv)]
cat(sprintf("samples %d | distinct ASVs %d | bacterial %d\n", nsamp, length(pv), sum(b)))
cat(sprintf("\n%12s %10s %14s %16s\n","min samples","ASVs","%% of reads","matrix cells"))
for (m in c(2,5,10,20,50,100)) {
  k <- b & pv >= m
  cat(sprintf("%12d %10d %13.2f%% %15.0fM\n", m, sum(k), 100*sum(tv[k])/sum(tv[b]),
      sum(k)*nsamp/1e6))
}
if (!BUILD) {
  cat(sprintf("\nDRY RUN. Chosen threshold would be %d samples.\n", MINS))
  cat("Re-run with --build (optionally --minsamples=N) to write the files.\n")
  quit(status=0)
}

sel <- names(pv)[b & pv >= MINS]
dir.create("picrust", showWarnings=FALSE)
ids <- paste0("ASV", seq_along(sel)); names(ids) <- sel
writeLines(as.vector(rbind(paste0(">", ids), unname(seqs[sel]))), "picrust/asv.fna")

A <- matrix(0L, length(sel), nsamp, dimnames=list(ids, NULL))
sn <- character(0); j <- 0L
for (p in f) {
  st <- readRDS(p); st <- st[rownames(st) %in% KEEPS, , drop=FALSE]
  if (!nrow(st)) next
  cm <- intersect(colnames(st), sel)
  A[ids[cm], (j+1):(j+nrow(st))] <- t(st[, cm, drop=FALSE])
  sn <- c(sn, rownames(st)); j <- j + nrow(st)
}
colnames(A) <- sn
write.table(cbind(`#OTU ID`=rownames(A), A), "picrust/asv_table.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(asv_id=ids, sequence=sel, n_samples=pv[sel], total_reads=tv[sel]),
            "picrust/asv_map.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("\nwritten: picrust/asv.fna (%d seqs, %d-%d bp)\n", length(sel),
    min(nchar(seqs[sel])), max(nchar(seqs[sel]))))
cat(sprintf("         picrust/asv_table.tsv (%d x %d, %.0f MB)\n", length(sel), nsamp,
    file.size("picrust/asv_table.tsv")/1e6))

cat("\n--- reads retained per cohort (uniform = filter adds no batch effect) ---\n")
md <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
md <- md[match(colnames(A), md$run_accession), ]
kept <- colSums(A); tt <- numeric(ncol(A)); j <- 0L
for (p in f) { st <- readRDS(p); st <- st[rownames(st) %in% KEEPS, , drop=FALSE]
  if (!nrow(st)) next
  tt[(j+1):(j+nrow(st))] <- rowSums(st); j <- j + nrow(st) }
fr <- 100*kept/tt
q <- tapply(fr, md$cohort_id, median)
cat(sprintf("median retained %.1f%% | range across cohorts %.1f%% - %.1f%%\n",
    median(fr), min(q), max(q)))
print(round(sort(q)[1:8], 1))
