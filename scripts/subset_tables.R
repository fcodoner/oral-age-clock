#!/usr/bin/env Rscript
# --------------------------------------------------------------------------
# subset_tables.R
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
WRITE <- "--write" %in% args
MINDEPTH <- 3000

part <- read.delim("metadata/partition.tsv", stringsAsFactors=FALSE)
cls  <- setNames(part$class, part$run_accession)

files <- list.files("final", pattern="\\.rds$", full.names=TRUE)
if (!length(files)) stop("no .rds files in final/ -- list the directory and tell me what is there")

if (WRITE) { dir.create("final_split/train", recursive=TRUE, showWarnings=FALSE)
             dir.create("final_split/valid", recursive=TRUE, showWarnings=FALSE) }

norm <- function(x) sub("_(F|R)?_?filt.*$","", sub("\\.fastq.*$","", sub("_1$","", x)))

tot <- c(train=0, valid=0, unmatched=0, shallow=0); rowsall <- 0L
cat(sprintf("%-26s %6s %6s %6s %6s %7s %7s\n",
            "group","rows","match","train","valid","unmatch","shallow"))
for (f in files) {
  g  <- sub("\\.rds$","", basename(f))
  st <- readRDS(f)
  if (is.null(dim(st))) { cat(sprintf("%-26s  NOT A MATRIX\n", g)); next }
  rn <- norm(rownames(st))
  k  <- cls[rn]
  dep <- rowSums(st)
  low <- !is.na(k) & dep < MINDEPTH
  nt <- sum(k=="train" & !low, na.rm=TRUE); nv <- sum(k=="validate" & !low, na.rm=TRUE)
  nu <- sum(is.na(k))
  cat(sprintf("%-26s %6d %6d %6d %6d %7d %7d%s\n", g, nrow(st), sum(!is.na(k)),
              nt, nv, nu, sum(low), if (nu==nrow(st)) "   <-- NO NAMES MATCHED" else ""))
  tot <- tot + c(nt, nv, nu, sum(low)); rowsall <- rowsall + nrow(st)
  if (WRITE) {
    for (cl in c("train","validate")) {
      sel <- which(k==cl & !low)
      if (!length(sel)) next
      sub <- st[sel, , drop=FALSE]
      sub <- sub[, colSums(sub) > 0, drop=FALSE]
      rownames(sub) <- rn[sel]
      saveRDS(sub, file.path(if (cl=="train") "final_split/train" else "final_split/valid",
                             paste0(g, ".rds")))
    }
  }
}
cat(sprintf("\nfiles %d | table rows %d\nTRAIN %d | VALID %d | unmatched %d | below %d reads %d\n",
            length(files), rowsall, tot[["train"]], tot[["valid"]],
            tot[["unmatched"]], MINDEPTH, tot[["shallow"]]))
if (WRITE) cat("written: final_split/train/ and final_split/valid/\n")
