# --------------------------------------------------------------------------
# revision_round2.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# revision_round2.R
# Re-runs the paragraph [29] sensitivity analyses under the MAIN pipeline
# (118 genera, CLR pseudocount 0.5, eligible-fold rule) so that they are
# comparable with the headline leave-one-cohort-out result.
# F. M. Codoner
# ---------------------------------------------------------------------------
suppressMessages({library(glmnet); library(vegan)})

src <- readLines("scripts/reported_values.R")
j <- grep("^folds *<-", src)[1]
k <- grep("^\\s*R\\$", src)
i <- k[k > j][1] - 1          # stop just before the first result assignment
stopifnot(!is.na(j), !is.na(i), i > j)
cat(sprintf("evaluating preamble: lines 1-%d of reported_values.R\n", i))
eval(parse(text = paste(src[1:i], collapse = "\n")), envir = globalenv())
stopifnot("loco() was not defined by the preamble" = is.function(loco),
          "fit() was not defined by the preamble"  = is.function(fit))

cat("###########################################################\n")
cat("# 0. METADATA DIAGNOSTICS                                 #\n")
cat("###########################################################\n\n")
cat("columns in sm:\n"); print(names(sm))

cat("\nsearching every text column for 'tongue':\n")
hit <- FALSE
for (cl in names(sm)) {
  v <- sm[[cl]]
  if (is.character(v) || is.factor(v)) {
    v <- as.character(v)
    ix <- grepl("tongue", v, ignore.case = TRUE)
    if (any(ix, na.rm = TRUE)) {
      hit <- TRUE
      cat(sprintf("  '%s' -> cohorts: %s\n", cl,
                  paste(sort(unique(sm$cohort_id[ix])), collapse = ", ")))
      cat(sprintf("     values: %s\n", paste(sort(unique(v[ix])), collapse = " | ")))
    }
  }
}
if (!hit) cat("  no column mentions 'tongue'\n")

sexcol <- intersect(c("sex", "gender", "host_sex"), names(sm))
cat("\nsex column:", if (length(sexcol)) sexcol[1] else "NOT FOUND", "\n")
if (length(sexcol)) {
  print(table(sm[[sexcol[1]]][tr], useNA = "ifany"))
}

cat("\nspecimen route by cohort:\n")
if ("route" %in% names(sm)) {
  print(table(route = sm$route[tr], useNA = "ifany"))
  cat("\ncohorts per route:\n")
  for (rt in sort(unique(as.character(sm$route[tr])))) {
    cat(sprintf("  %-18s %s\n", rt,
        paste(sort(unique(co[as.character(sm$route[tr]) == rt])), collapse = ", ")))
  }
} else cat("  no route column\n")

if ("region" %in% names(sm)) {
  cat("\namplicon region by cohort count:\n")
  print(table(region = sm$region[tr], useNA = "ifany"))
}

cat("\ncohort sizes among training samples (largest first):\n")
print(head(sort(table(co), decreasing = TRUE), 8))

# --------------------------------------------------------------------------
sens <- function(label, keepvec, cohortvec = NULL, extracol = NULL) {
  X <- XG[tr, , drop = FALSE][keepvec, , drop = FALSE]
  if (!is.null(extracol)) X <- cbind(X, extra = extracol[keepvec])
  yy <- y[keepvec]
  cc <- if (is.null(cohortvec)) co[keepvec] else cohortvec[keepvec]
  a  <- tapply(yy, cc, function(v) diff(range(v)))
  n  <- table(cc)
  ff <- names(a)[a >= 40 & n[names(a)] >= 40]
  if (!length(ff)) { cat(sprintf("%-34s no eligible folds\n", label)); return(invisible(NULL)) }
  r <- loco(X, yy, cc, ff)
  cat(sprintf("%-34s n=%4d  folds=%2d  MAE %6.3f  null %6.3f  r %.3f  R2 %6.3f  medWithinR %.3f\n",
              label, length(yy), length(ff),
              r["MAE"], r["nullMAE"], r["r"], r["R2"], r["medWithinR"]))
  invisible(r)
}

cat("\n###########################################################\n")
cat("# 1. BASELINE - MUST REPRODUCE age_loco                   #\n")
cat("###########################################################\n\n")
cat("stored age_loco: MAE 12.003  null 14.444  r 0.470  R2 0.199  medWithinR 0.363\n\n")
R0 <- sens("baseline (all training samples)", rep(TRUE, length(y)))

cat("\n###########################################################\n")
cat("# 2. SEX AS A COVARIATE                                   #\n")
cat("###########################################################\n\n")
if (length(sexcol)) {
  sv <- as.character(sm[[sexcol[1]]][tr])
  has <- sv %in% c("male", "female", "Male", "Female", "M", "F", "m", "f")
  cat(sprintf("samples with sex recorded: %d of %d\n\n", sum(has), length(has)))
  male <- as.numeric(tolower(substr(sv, 1, 1)) == "m")
  sens("sex subset, without sex term", has)
  sens("sex subset, with sex term",    has, extracol = male)
} else cat("no sex column available; analysis cannot be run\n")

cat("\n###########################################################\n")
cat("# 3. COHORTS FROM ONE LABORATORY POOLED INTO ONE FOLD     #\n")
cat("###########################################################\n\n")
c2 <- co
c2[c2 %in% c("USA_39492", "USA_74379")] <- "UIC_POOLED"
cat(sprintf("pooled fold size: %d samples\n\n", sum(c2 == "UIC_POOLED")))
sens("UIC cohorts pooled", rep(TRUE, length(y)), cohortvec = c2)

cat("\n###########################################################\n")
cat("# 4. LARGEST COHORT EXCLUDED  (reproduces age_noAGP)      #\n")
cat("###########################################################\n\n")
cat("stored age_noAGP: MAE 11.805  null 14.792  r 0.520  R2 0.245  medWithinR 0.363\n\n")
sens("AGP_ORAL excluded", co != "AGP_ORAL")
big <- names(sort(table(co), decreasing = TRUE))[1]
if (big != "AGP_ORAL") {
  cat(sprintf("\nnote: the largest training cohort is %s, not AGP_ORAL\n", big))
  sens(sprintf("%s excluded", big), co != big)
}

cat("\n###########################################################\n")
cat("# 5. TONGUE-SWAB COHORT EXCLUDED                          #\n")
cat("###########################################################\n\n")
tc <- character(0)
for (cl in names(sm)) {
  v <- sm[[cl]]
  if (is.character(v) || is.factor(v)) {
    ix <- grepl("tongue", as.character(v), ignore.case = TRUE)
    if (any(ix, na.rm = TRUE)) tc <- union(tc, unique(sm$cohort_id[ix]))
  }
}
if (length(tc)) {
  cat("tongue-associated cohorts:", paste(tc, collapse = ", "), "\n\n")
  sens(paste("excluding", paste(tc, collapse = "+")), !(co %in% tc))
} else cat("no tongue-swab cohort identifiable from the metadata\n")

cat("\n###########################################################\n")
cat("# 6. PROVENANCE OF THE 0.219 FIGURE                       #\n")
cat("###########################################################\n\n")
f <- "results/reviewer_checks.rds"
if (file.exists(f)) {
  RK <- readRDS(f)
  G2 <- readRDS("features/genus_counts.rds")
  SM2 <- read.delim("features/sample_metadata.tsv", stringsAsFactors = FALSE)
  rownames(G2) <- unname(rownames(G2))
  SM2 <- SM2[match(rownames(G2), SM2$run_accession), ]
  H <- !nzchar(SM2$label) & !is.na(SM2$age)
  y2 <- SM2$age[H]; coh2 <- SM2$cohort_id[H]
  p <- RK$p_base
  cat(sprintf("reviewer_checks pipeline: %d samples, %d cohorts\n", length(y2), length(unique(coh2))))
  if (length(p) == length(y2)) {
    wr <- sapply(split(seq_along(y2), coh2), function(ii)
      if (length(ii) >= 3) suppressWarnings(cor(p[ii], y2[ii])) else NA_real_)
    cat(sprintf("median within-cohort r over ALL %d cohorts : %.3f\n",
                sum(!is.na(wr)), median(wr, na.rm = TRUE)))
    el <- names(wr) %in% folds
    cat(sprintf("median within-cohort r over the %d eligible : %.3f\n",
                sum(el & !is.na(wr)), median(wr[el], na.rm = TRUE)))
    cat("\n-> if the first value is near 0.219 the manuscript figure came from\n")
    cat("   this pipeline across all cohorts, not from the main pipeline.\n")
  } else {
    cat(sprintf("length mismatch (p_base %d vs metadata %d); cannot align\n",
                length(p), length(y2)))
  }
} else cat("results/reviewer_checks.rds not found\n")

cat("\ndone\n")
