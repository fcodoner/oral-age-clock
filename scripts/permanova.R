# --------------------------------------------------------------------------
# permanova.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
if (!requireNamespace("vegan", quietly=TRUE))
  install.packages("vegan", repos="https://cloud.r-project.org")
suppressMessages(library(vegan))

R  <- readRDS("features/genus_relab.rds")
rownames(R) <- unname(rownames(R))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(R), sm$run_accession), ]
stopifnot(identical(sm$run_accession, rownames(R)))

prev <- colMeans(R > 0)
R <- R[, prev > 0.10, drop=FALSE]
R <- R / rowSums(R)
cat(sprintf("samples %d | genera after prevalence>10%%: %d\n", nrow(R), ncol(R)))

sm$sex2 <- ifelse(is.na(sm$sex) | sm$sex=="" , "unknown", tolower(sm$sex))
sm$age_grp <- cut(sm$age, c(17,30,45,60,120), labels=c("18-30","31-45","46-60","61+"))

cat("\n--- age distribution by cohort (10 widest, 10 narrowest) ---\n")
ag <- do.call(rbind, lapply(split(sm$age, sm$cohort_id), function(x)
      data.frame(n=length(x), min=min(x), med=median(x), max=max(x), span=diff(range(x)))))
ag <- ag[order(-ag$span), ]
print(head(ag, 10)); print(tail(ag, 10))

cat("\ncomputing Bray-Curtis distance...\n")
D <- vegdist(R, method="bray")

cat("\n--- age alone ---\n")
print(adonis2(D ~ age, data=sm, permutations=199))

cat("\n--- sequential: study terms first, then age ---\n")
A1 <- adonis2(D ~ region + route + cohort_id + age + sex2,
              data=sm, permutations=199, by="terms")
print(A1)

cat("\n--- sequential: age first, then study terms ---\n")
A2 <- adonis2(D ~ age + sex2 + region + route + cohort_id,
              data=sm, permutations=199, by="terms")
print(A2)

dir.create("results", showWarnings=FALSE)
saveRDS(list(study_first=A1, age_first=A2, n=nrow(R), genera=ncol(R)),
        "results/permanova.rds")
cat("\nsaved: results/permanova.rds\n")
