# --------------------------------------------------------------------------
# orient_pass2.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages(library(dada2)); suppressMessages(library(Biostrings))
TRAIN <- "ref/eHOMD_RefSeq_dada2_V15.22.fasta.gz"
SPEC  <- "ref/eHOMD_RefSeq_dada2_assign_species_V15.22.fasta.gz"
CHUNK <- 20000L; THREADS <- 3L
dir.create("taxonomy/rev", showWarnings=FALSE, recursive=TRUE)

obj <- readRDS("taxonomy/taxa_pooled.rds")
tax <- obj$tax; boot <- obj$boot
asv <- rownames(tax)
need <- which(is.na(tax[,"Phylum"]))
cat(sprintf("unassigned at phylum: %d / %d\n", length(need), nrow(tax)))

rc  <- as.character(reverseComplement(DNAStringSet(asv[need])))
idx <- split(seq_along(need), ceiling(seq_along(need)/CHUNK))
rt <- vector("list", length(idx)); rb <- vector("list", length(idx))
for (i in seq_along(idx)) {
  cf <- sprintf("taxonomy/rev/chunk_%03d.rds", i)
  if (file.exists(cf)) { o <- readRDS(cf)
    if (nrow(o$tax) == length(idx[[i]])) {
      rt[[i]] <- o$tax; rb[[i]] <- o$boot
      cat(sprintf("rev chunk %d/%d cached\n", i, length(idx))); flush.console(); next } }
  t0 <- Sys.time(); set.seed(500 + i)
  o <- assignTaxonomy(rc[idx[[i]]], TRAIN, minBoot=50,
                      multithread=THREADS, outputBootstraps=TRUE, verbose=FALSE)
  saveRDS(o, cf); rt[[i]] <- o$tax; rb[[i]] <- o$boot
  cat(sprintf("rev chunk %d/%d n=%d %.1f min\n", i, length(idx), length(idx[[i]]),
      as.numeric(difftime(Sys.time(), t0, units="mins")))); flush.console()
}
rtax <- do.call(rbind, rt); rboo <- do.call(rbind, rb)

ok  <- !is.na(rtax[,"Phylum"])
hit <- need[ok]
cat(sprintf("\nrecovered by reverse complement: %d (%.1f%% of unassigned)\n",
            sum(ok), 100*mean(ok)))

rk <- intersect(colnames(rtax), colnames(tax))
bk <- intersect(colnames(rboo), colnames(boot))
cat("ranks copied:", paste(rk, collapse=", "), "\n")
tax[hit, rk]  <- rtax[ok, rk, drop=FALSE]
boot[hit, bk] <- rboo[ok, bk, drop=FALSE]

orient <- rep("fwd", nrow(tax)); names(orient) <- asv
orient[hit] <- "rev"
rcseq <- asv; rcseq[hit] <- rc[ok]

if ("Species" %in% colnames(tax)) {
  cat("adding species for recovered ASVs...\n")
  sp <- addSpecies(rtax[ok, , drop=FALSE], SPEC, allowMultiple=TRUE, n=50000, verbose=TRUE)
  tax[hit, "Species"] <- if ("Species" %in% colnames(sp)) sp[, "Species"] else NA
}

saveRDS(list(tax=tax, boot=boot, orient=orient, oriented_seq=rcseq),
        "taxonomy/taxa_oriented.rds")
write.table(data.frame(asv=rownames(tax), orient=orient, tax, check.names=FALSE),
            "taxonomy/taxa_oriented.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("\nfinal unassigned at phylum: %d / %d (%.1f%%)\n",
    sum(is.na(tax[,"Phylum"])), nrow(tax), 100*mean(is.na(tax[,"Phylum"]))))
cat("orientation:", sum(orient=="fwd"), "fwd /", sum(orient=="rev"), "rev\n")
