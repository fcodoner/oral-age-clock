# --------------------------------------------------------------------------
# export_S6.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
# Exports Supplementary Table S6 from the meta-analysis result object.
# Writes results/S6_meta.csv. Run from the project root.
M <- readRDS("results/meta_detail.rds")

out <- do.call(rbind, lapply(names(M), function(k) {
  d <- M[[k]]
  d <- d[!is.na(d$q) & d$q < 0.05, , drop = FALSE]
  d <- d[order(d$q), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  tot <- as.integer(sub(".*/", "", d$agree))
  agr <- as.integer(sub("/.*",  "", d$agree))
  data.frame(
    Type                          = if (k == "genus") "Genus" else "Inferred pathway",
    Feature                       = if ("feature_full" %in% names(d)) d$feature_full else d$feature,
    `Pooled Spearman rho`         = sprintf("%.3f", d$rho),
    Direction                     = ifelse(d$rho > 0, "increases with age", "decreases with age"),
    `95% CI lower`                = sprintf("%.3f", d$lo),
    `95% CI upper`                = sprintf("%.3f", d$hi),
    `Cohorts contributing`        = tot,
    `Cohorts in pooled direction` = agr,
    `q value (BH)`                = signif(d$q, 3),
    check.names = FALSE, stringsAsFactors = FALSE)
}))

write.csv(out, "results/S6_meta.csv", row.names = FALSE)
cat(sprintf("wrote results/S6_meta.csv: %d rows (%d genera, %d pathways)\n",
            nrow(out), sum(out$Type == "Genus"), sum(out$Type == "Inferred pathway")))
cat(sprintf("cohorts contributing: range %d to %d\n",
            min(out$`Cohorts contributing`), max(out$`Cohorts contributing`)))
print(out, row.names = FALSE)
