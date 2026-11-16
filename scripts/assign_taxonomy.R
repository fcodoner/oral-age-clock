#!/usr/bin/env Rscript
# --------------------------------------------------------------------------
# assign_taxonomy.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages(library(dada2))
args  <- commandArgs(TRUE)
RUN   <- "--run" %in% args
CHUNK <- 20000L
THREADS <- max(1L, parallel::detectCores() - 1L)

rf    <- list.files("ref", full.names=TRUE)
train <- "ref/eHOMD_RefSeq_dada2_V15.22.fasta.gz"
spec  <- "ref/eHOMD_RefSeq_dada2_assign_species_V15.22.fasta.gz"
stopifnot(file.exists(train), file.exists(spec))

cat("--- reference files in ref/ ---\n"); for (f in rf) cat("  ", f, "\n")
cat("\nDETECTED training set :", if (length(train)) train else "NONE", "\n")
cat("DETECTED species file :", if (length(spec)) spec else "NONE", "\n")
if (length(train) != 1) stop("could not pick exactly one training set -- tell me the filenames")

asv <- unique(unlist(lapply(
  list.files(c("final_split/train","final_split/valid"), pattern="\\.rds$", full.names=TRUE),
  function(x) colnames(readRDS(x)))))
cat(sprintf("\npooled unique ASVs: %d\nlength range: %d-%d bp\nchunks of %d: %d\nthreads: %d\n",
    length(asv), min(nchar(asv)), max(nchar(asv)), CHUNK, ceiling(length(asv)/CHUNK), THREADS))

if (!RUN) { cat("\nDRY RUN -- nothing classified. Re-run with --run to start.\n"); quit(status=0) }

dir.create("taxonomy", showWarnings=FALSE)
saveRDS(asv, "taxonomy/pooled_asvs.rds")
idx <- split(seq_along(asv), ceiling(seq_along(asv)/CHUNK))
out <- vector("list", length(idx))
for (i in seq_along(idx)) {
  t0 <- Sys.time()
  set.seed(100 + i)
  out[[i]] <- assignTaxonomy(asv[idx[[i]]], train, minBoot=50,
                             multithread=THREADS, outputBootstraps=TRUE, verbose=FALSE)
  saveRDS(out[[i]], sprintf("taxonomy/chunk_%03d.rds", i))
  cat(sprintf("chunk %d/%d  n=%d  %.1f min\n", i, length(idx), length(idx[[i]]),
      as.numeric(difftime(Sys.time(), t0, units="mins")))); flush.console()
}
tax <- do.call(rbind, lapply(out, function(x) x$tax))
boo <- do.call(rbind, lapply(out, function(x) x$boot))
rownames(tax) <- rownames(boo) <- asv
if (length(spec) == 1) {
  cat("adding species-level assignment...\n")
  tax <- addSpecies(tax, spec, allowMultiple=TRUE, n=50000, verbose=TRUE)
}
saveRDS(list(tax=tax, boot=boo), "taxonomy/taxa_pooled.rds")
write.table(data.frame(asv=rownames(tax), tax, check.names=FALSE),
            "taxonomy/taxa_pooled.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("\ndone. ranks:", paste(colnames(tax), collapse=", "), "\n")
cat("unassigned at phylum:", sum(is.na(tax[,2])), "/", nrow(tax),
    sprintf("(%.1f%%)\n", 100*mean(is.na(tax[,2]))))
