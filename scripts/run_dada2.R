#!/usr/bin/env Rscript
# --------------------------------------------------------------------------
# run_dada2.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# run_dada2.R - denoise one processing group with DADA2 and write its ASV table.
#
#   Rscript scripts/run_dada2.R <group name> \
#           --params metadata/dada2/dada2_params.tsv \
#           --out processed
#
#   # or every group, one after another:
#   Rscript scripts/run_dada2.R ALL --params metadata/dada2/dada2_params.tsv --out processed
#
# Each group is denoised on its own, which is what DADA2 requires: the error model
# is learned per sequencing run, not across studies. The per-group ASV tables are
# merged later, in a separate step.
#
# Needs: dada2 (>= 1.26). Install with:
#   conda install -c bioconda bioconductor-dada2

suppressPackageStartupMessages(library(dada2))

# ---------------------------------------------------------------------------
# Binned quality scores (NextSeq / NovaSeq).
#
# These instruments report only a handful of distinct quality values instead of
# the full 0-40 range. DADA2's standard error model fits a loess curve through
# those values; with 4 points the fit fails and learnErrors() returns NULL,
# which surfaces as "Error matrix is NULL". The replacement below widens the
# loess span, fills the gaps at both ends of the range, and forces the error
# rate to decline monotonically with quality - the published workaround for
# binned data. Everything else about the run is unchanged.
# ---------------------------------------------------------------------------
loessErrfun_binned <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) for (ntj in c("A", "C", "G", "T")) if (nti != ntj) {
    errs  <- trans[paste0(nti, "2", ntj), ]
    tot   <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
    rlogp <- log10((errs + 1) / tot)
    rlogp[is.infinite(rlogp)] <- NA
    df <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
    mod.lo <- tryCatch(loess(rlogp ~ q, df, weights = log10(tot), span = 2),
                       error = function(e) NULL)
    pred <- if (is.null(mod.lo)) rep(NA_real_, length(qq)) else predict(mod.lo, qq)
    if (all(is.na(pred))) pred <- rep(mean(rlogp, na.rm = TRUE), length(qq))
    maxrli <- max(which(!is.na(pred))); minrli <- min(which(!is.na(pred)))
    pred[seq_along(pred) > maxrli] <- pred[maxrli]
    pred[seq_along(pred) < minrli] <- pred[minrli]
    est <- rbind(est, 10^pred)
  }
  MAX_ERROR_RATE <- 0.25; MIN_ERROR_RATE <- 1e-7
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
  est <- matrix(t(apply(est, 1, function(x) {
    if (length(x) > 1) for (i in 2:length(x)) if (x[i] > x[i - 1]) x[i] <- x[i - 1]
    x
  })), nrow = 12)
  s <- function(i) est[i, , drop = FALSE]
  err <- rbind(1 - colSums(s(1:3)), s(1:3),
               s(4), 1 - colSums(s(4:6)), s(5:6),
               s(7:8), 1 - colSums(s(7:9)), s(9),
               s(10:12), 1 - colSums(s(10:12)))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  err
}

# ---------------------------------------------------------------------------
# No quality information at all.
#
# Some deposits carry a single quality character for every base: the quality
# strings were stripped or rebuilt before submission. No quality-based model can
# be fitted to that, so this one ignores quality and estimates a single rate per
# substitution type (A->C, A->G, ...) from the observed transition counts, held
# constant across the quality range. Denoising then rests on sequence abundance
# rather than per-base quality - a limitation to report, not an error.
# ---------------------------------------------------------------------------
noqualErrfun <- function(trans) {
  nts <- c("A", "C", "G", "T")
  err <- matrix(0, nrow = 16, ncol = ncol(trans),
                dimnames = list(paste0(rep(nts, each = 4), "2", nts), colnames(trans)))
  for (nti in nts) {
    tot <- sum(trans[paste0(nti, "2", nts), , drop = FALSE])
    for (ntj in setdiff(nts, nti)) {
      e <- sum(trans[paste0(nti, "2", ntj), , drop = FALSE])
      r <- min(max((e + 1) / (tot + 4), 1e-7), 0.25)
      err[paste0(nti, "2", ntj), ] <- r
    }
    err[paste0(nti, "2", nti), ] <- 1 - sum(err[paste0(nti, "2", setdiff(nts, nti)), 1])
  }
  err
}

# how many distinct quality characters a FASTQ actually uses
n_quality_levels <- function(path, n_reads = 2000) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  lines <- readLines(con, n = n_reads * 4, warn = FALSE)
  if (length(lines) < 4) return(NA_integer_)
  quals <- lines[seq(4, length(lines), by = 4)]
  length(unique(unlist(strsplit(paste(quals, collapse = ""), ""))))
}

# learnErrors with an automatic fallback for binned data
learn <- function(files, threads, label = "") {
  nq <- suppressWarnings(n_quality_levels(files[1]))
  # three models, tried in the order the data call for
  if (!is.na(nq) && nq <= 1) {
    order <- c("noqual", "binned")
    message("    ", label, "a single quality value throughout - quality carries no ",
            "information here, using the quality-independent error model")
  } else if (!is.na(nq) && nq <= 8) {
    order <- c("binned", "noqual")
    message("    ", label, "only ", nq, " distinct quality values - using the binned error model")
  } else {
    order <- c("standard", "binned", "noqual")
  }
  fns <- list(standard = NULL, binned = loessErrfun_binned, noqual = noqualErrfun)
  err <- NULL
  for (nm in order) {
    if (!is.null(err)) break
    if (nm != order[1])
      message("    ", label, "falling back to the ", nm, " error model")
    err <- tryCatch({
      e <- if (is.null(fns[[nm]]))
             learnErrors(files, multithread = threads, verbose = FALSE)
           else
             learnErrors(files, multithread = threads, verbose = FALSE,
                         errorEstimationFunction = fns[[nm]])
      if (is.null(e) || anyNA(getErrors(e))) NULL else e
    }, error = function(e) NULL)
  }
  if (is.null(err)) stop("could not learn an error model for ", label, "reads")
  err
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_dada2.R <group|ALL> [--params FILE] [--out DIR]")
group_arg <- args[1]
getopt <- function(flag, default) {
  i <- match(flag, args); if (!is.na(i) && length(args) > i) args[i + 1] else default
}
params_file <- getopt("--params", "metadata/dada2/dada2_params.tsv")
outdir      <- getopt("--out", "processed")
threads     <- as.integer(getopt("--threads", "0"))      # 0 = all available
if (threads == 0) threads <- TRUE

params <- read.delim(params_file, stringsAsFactors = FALSE)
groups <- if (group_arg == "ALL") params$group else group_arg
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "qc"), showWarnings = FALSE)

process_group <- function(gname) {
  p <- params[params$group == gname, ]
  if (nrow(p) == 0) { message("No parameters for group ", gname); return(invisible(NULL)) }
  p <- p[1, ]
  done <- file.path(outdir, paste0(gname, "_seqtab.rds"))
  if (file.exists(done)) { message("[", gname, "] already done, skipping"); return(invisible(NULL)) }

  message("\n=== ", gname, " | ", p$layout, " | ", p$region, " | ", p$n_runs, " runs ===")
  man <- read.delim(p$manifest, stringsAsFactors = FALSE)
  man <- man[file.exists(man$fwd), ]
  if (nrow(man) == 0) { message("  no readable files, skipped"); return(invisible(NULL)) }
  forward_only <- isTRUE(!is.na(p$forward_only) && p$forward_only == "yes")
  paired <- (p$layout == "paired") && !forward_only &&
            all(nzchar(man$rev)) && all(file.exists(man$rev))
  if (forward_only) message("  forward reads only (the pair cannot overlap for this amplicon)")
  minLen <- if (!is.null(p$minLen) && !is.na(p$minLen)) as.integer(p$minLen) else 50
  # on an already-merged amplicon the REVERSE primer sits at the 3' end, and nothing
  # else removes it; left in, it inflates chimeras and splits ASVs
  trimRight <- if (!is.null(p$trimRight) && !is.na(p$trimRight)) as.integer(p$trimRight) else 0
  if (trimRight > 0) message("  trimming ", trimRight, " bp from the 3' end (reverse primer)")

  filtdir <- file.path(outdir, "filtered", gname)
  dir.create(filtdir, showWarnings = FALSE, recursive = TRUE)
  fF <- file.path(filtdir, paste0(man$run_accession, "_F_filt.fastq.gz"))
  fR <- if (paired) file.path(filtdir, paste0(man$run_accession, "_R_filt.fastq.gz")) else NULL

  # quality profile of the first few samples, for the record
  pdf(file.path(outdir, "qc", paste0(gname, "_quality.pdf")), width = 9, height = 5)
  try(print(plotQualityProfile(man$fwd[1:min(4, nrow(man))])), silent = TRUE)
  if (paired) try(print(plotQualityProfile(man$rev[1:min(4, nrow(man))])), silent = TRUE)
  dev.off()

  # 454 and other long reads need a wider alignment band and no homopolymer penalty
  long_read <- p$layout == "single_long"
  if (long_read) {
    setDadaOpt(HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32)
    message("  long-read settings applied (HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)")
  }

  message("  filtering and trimming...")
  if (paired) {
    out <- filterAndTrim(man$fwd, fF, man$rev, fR,
                         truncLen = c(p$truncLen_fwd, p$truncLen_rev),
                         trimLeft = c(p$trimLeft, p$trimLeft),
                         maxEE = c(as.numeric(p$maxEE_fwd), as.numeric(p$maxEE_rev)),
                         truncQ = p$truncQ, minLen = minLen, maxN = 0, rm.phix = TRUE,
                         compress = TRUE, multithread = threads, verbose = FALSE)
  } else {
    # truncLen 0 means "do not truncate": used for 454 and already-merged amplicons,
    # where reads legitimately vary in length and minLen does the filtering instead
    tl <- as.integer(p$truncLen_fwd)
    out <- filterAndTrim(man$fwd, fF,
                         truncLen = tl, trimLeft = as.integer(p$trimLeft),
                         trimRight = trimRight,
                         maxEE = as.numeric(p$maxEE_fwd), truncQ = p$truncQ,
                         minLen = minLen, maxN = 0,
                         rm.phix = TRUE, compress = TRUE,
                         multithread = threads, verbose = FALSE)
  }

  keep <- file.exists(fF)
  if (sum(keep) == 0) { message("  every sample lost all reads - check truncLen"); return(invisible(NULL)) }
  if (sum(keep) < length(fF))
    message("  ", length(fF) - sum(keep), " sample(s) lost all reads and were dropped")
  fF <- fF[keep]; if (paired) fR <- fR[keep]
  sample_names <- man$run_accession[keep]

  message("  learning the error model...")
  errF <- learn(fF, threads, "forward: ")
  if (paired) errR <- learn(fR, threads, "reverse: ")
  pdf(file.path(outdir, "qc", paste0(gname, "_errors.pdf")), width = 8, height = 6)
  try(print(plotErrors(errF, nominalQ = TRUE)), silent = TRUE)
  dev.off()

  message("  denoising...")
  ddF <- dada(fF, err = errF, multithread = threads, verbose = FALSE)
  if (paired) {
    ddR <- dada(fR, err = errR, multithread = threads, verbose = FALSE)
    mrg <- mergePairs(ddF, fF, ddR, fR, minOverlap = 12, verbose = FALSE)
    seqtab <- makeSequenceTable(mrg)
  } else {
    seqtab <- makeSequenceTable(ddF)
  }
  rownames(seqtab) <- sample_names

  message("  removing chimeras...")
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = threads, verbose = FALSE)

  # how many reads survived each step: the table reviewers ask for
  getN <- function(x) sum(getUniques(x))
  track <- data.frame(
    run_accession = sample_names,
    input     = out[keep, 1],
    filtered  = out[keep, 2],
    denoisedF = sapply(ddF, getN),
    merged    = if (paired) sapply(mrg, getN) else NA,
    nonchim   = rowSums(seqtab_nochim),
    stringsAsFactors = FALSE)
  track$pct_kept <- round(100 * track$nonchim / track$input, 1)

  saveRDS(seqtab_nochim, done)
  write.table(track, file.path(outdir, paste0(gname, "_track.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (forward_only) attr(seqtab_nochim, "forward_only") <- TRUE
  message("  ", ncol(seqtab_nochim), " ASVs across ", nrow(seqtab_nochim), " samples; ",
          "median ", round(median(track$pct_kept, na.rm = TRUE), 1), "% of reads kept")
  if (median(track$pct_kept, na.rm = TRUE) < 40)
    message("  WARNING: under 40% of reads survived - revisit truncLen for this group")
  invisible(NULL)
}

for (g in groups) {
  tryCatch(process_group(g),
           error = function(e) message("[", g, "] FAILED: ", conditionMessage(e)))
}
message("\nDone. ASV tables in ", outdir, "/, read-tracking tables alongside them.")
