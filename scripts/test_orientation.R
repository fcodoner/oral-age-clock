# --------------------------------------------------------------------------
# test_orientation.R
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
N <- 2000L
grp <- c("USA_39492__single__V1-V3","USA_74379__single__V1-V3","LEE17_CN__single__V4",
         "SPAIN_74796__paired__V1-V2","PORTUGAL_74247__paired__V1-V2",
         "CHINA_48248__single__V4","SOUTHKORE_04331__paired__V3-V4","YU20_CN__single__V4",
         "CHINA_48248__paired__V4")
rate <- function(s) { set.seed(1)
  t <- assignTaxonomy(s, TRAIN, minBoot=50, multithread=3, verbose=FALSE)
  100*mean(!is.na(t[,"Phylum"])) }
cat(sprintf("%-32s %6s %10s %10s %8s\n","group","n","as-is%","revcomp%","verdict"))
for (g in grp) {
  p <- file.path("final_split/train", paste0(g,"_seqtab.rds"))
  if (!file.exists(p)) p <- file.path("final_split/valid", paste0(g,"_seqtab.rds"))
  if (!file.exists(p)) { cat(sprintf("%-32s  missing\n", g)); next }
  st <- readRDS(p); a <- colSums(st)
  s <- names(sort(a, decreasing=TRUE))[seq_len(min(N, length(a)))]
  f <- rate(s)
  r <- rate(as.character(reverseComplement(DNAStringSet(s))))
  cat(sprintf("%-32s %6d %9.1f%% %9.1f%% %8s\n", g, length(s), f, r,
      if (r > f + 20) "FLIPPED" else if (f > 90) "ok" else "neither"))
}
