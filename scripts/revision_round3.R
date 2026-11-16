# --------------------------------------------------------------------------
# revision_round3.R
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
# revision_round3.R
# Specimen type: recovers the curated saliva/rinse vs swab classification and
# tests whether pooling specimen types accounts for the limited transfer.
# F. M. Codoner
# ---------------------------------------------------------------------------
suppressMessages({library(glmnet); library(vegan)})

src <- readLines("scripts/reported_values.R")
j <- grep("^folds *<-", src)[1]
k <- grep("^\\s*R\\$", src)
i <- k[k > j][1] - 1
stopifnot(!is.na(j), !is.na(i), i > j)
cat(sprintf("evaluating preamble: lines 1-%d\n", i))
eval(parse(text = paste(src[1:i], collapse = "\n")), envir = globalenv())
stopifnot(is.function(loco), is.function(fit))

sens <- function(label, keepvec, cohortvec = NULL, extracol = NULL) {
  X <- XG[tr, , drop = FALSE][keepvec, , drop = FALSE]
  if (!is.null(extracol)) X <- cbind(X, extra = extracol[keepvec])
  yy <- y[keepvec]
  cc <- if (is.null(cohortvec)) co[keepvec] else cohortvec[keepvec]
  a <- tapply(yy, cc, function(v) diff(range(v))); n <- table(cc)
  ff <- names(a)[a >= 40 & n[names(a)] >= 40]
  if (!length(ff)) {
    cat(sprintf("%-34s n=%4d  NO ELIGIBLE FOLDS\n", label, length(yy))); return(invisible(NULL)) }
  r <- loco(X, yy, cc, ff)
  cat(sprintf("%-34s n=%4d  folds=%2d  MAE %6.3f  null %6.3f  r %.3f  medWithinR %.3f\n",
              label, length(yy), length(ff), r["MAE"], r["nullMAE"], r["r"], r["medWithinR"]))
  invisible(r)
}

cat("\n=== 1. RECOVERING THE CURATED SPECIMEN CLASSIFICATION ===\n\n")
CS <- read.delim("metadata/curated_core/curated_samples.tsv", stringsAsFactors = FALSE)
cat("columns in curated_samples.tsv:\n"); print(names(CS))

idcol <- intersect(c("run_accession", "run", "accession", "sample_accession", "run_id"), names(CS))
stopifnot("no accession column found" = length(idcol) > 0)
cat("\njoining on:", idcol[1], "\n")

st_all <- CS$sample_type[match(sm$run_accession, CS[[idcol[1]]])]
cat(sprintf("matched %d of %d samples in sm\n", sum(!is.na(st_all)), nrow(sm)))
st <- st_all[tr]
cat("\nspecimen type among training samples:\n")
print(table(st, useNA = "ifany"))

cat("\nspecimen type by cohort (training samples):\n")
tb <- table(cohort = co, type = st, useNA = "ifany")
print(tb)
mix <- rowSums(tb > 0)
cat(sprintf("\ncohorts with more than one specimen type: %d of %d\n", sum(mix > 1), length(mix)))
if (any(mix > 1)) print(tb[mix > 1, , drop = FALSE])

cat("\nage by specimen type:\n")
for (lv in sort(unique(st[!is.na(st)]))) {
  v <- y[which(st == lv)]
  cat(sprintf("  %-14s n=%4d  mean age %5.1f  sd %5.1f  range %2.0f-%2.0f\n",
              lv, length(v), mean(v), sd(v), min(v), max(v)))
}

cat("\n=== 2. STRATIFIED LEAVE-ONE-COHORT-OUT BY SPECIMEN TYPE ===\n\n")
cat("primary analysis for reference: MAE 12.003  null 14.444  r 0.469  medWithinR 0.363\n\n")
for (lv in sort(unique(st[!is.na(st)]))) {
  sens(sprintf("%s only", lv), !is.na(st) & st == lv)
}

cat("\n=== 3. SPECIMEN TYPE AS A COVARIATE ===\n\n")
ok <- !is.na(st)
if (length(unique(st[ok])) > 1) {
  isswab <- as.numeric(st == "swab")
  sens("typed subset, no type term", ok)
  sens("typed subset, with type term", ok, extracol = isswab)
} else cat("only one specimen type present; covariate test not applicable\n")

cat("\n=== 4. SPECIMEN TYPE IN THE PERMANOVA CONTEXT ===\n\n")
if (exists("genus_relab")) {
  RM <- genus_relab
} else {
  RM <- readRDS("features/genus_relab.rds")
}
rn <- rownames(RM)
idx <- match(sm$run_accession[tr], rn)
good <- !is.na(idx) & ok
cat(sprintf("samples usable for distance-based test: %d\n", sum(good)))
set.seed(1)
sub <- which(good)
if (length(sub) > 1500) sub <- sort(sample(sub, 1500))
D <- vegdist(RM[idx[sub], , drop = FALSE], method = "bray")
md <- data.frame(type = factor(st[sub]), cohort = factor(co[sub]), age = y[sub])
cat(sprintf("subsample used: %d samples, %d cohorts\n\n", length(sub), nlevels(md$cohort)))
cat("specimen type alone:\n")
print(adonis2(D ~ type, data = md, permutations = 999))
cat("\nspecimen type after study of origin:\n")
print(adonis2(D ~ cohort + type, data = md, permutations = 999, by = "terms"))

saveRDS(list(sample_type = st_all, table = tb), "results/revision_round3.rds")
cat("\nsaved: results/revision_round3.rds\n")
