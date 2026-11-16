# --------------------------------------------------------------------------
# build_release.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages(library(glmnet))
V <- "1.0"; OUT <- file.path("release", paste0("oral-health-score-v", V))
for (d in c("models","scripts","reference","examples","docs"))
  dir.create(file.path(OUT,d), recursive=TRUE, showWarnings=FALSE)

B  <- readRDS("models/oralclock_final.rds")
M  <- readRDS("features/genus_counts.rds"); rownames(M) <- unname(rownames(M))
sm <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
sm <- sm[match(rownames(M), sm$run_accession), ]

file.copy("models/oralclock_final.rds", file.path(OUT,"models"), overwrite=TRUE)
file.copy("scripts/score.R",            file.path(OUT,"scripts"), overwrite=TRUE)

write.table(data.frame(genus=B$genera), file.path(OUT,"reference","genus_template.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cf <- B$coefs[names(B$coefs)!="(Intercept)"]
write.table(data.frame(genus=names(cf), coefficient=round(cf,4),
            direction=ifelse(cf>0,"raises score","lowers score")),
            file.path(OUT,"reference","coefficients.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(score=round(B$global_ref,4)),
            file.path(OUT,"reference","global_healthy_scores.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(genus=rownames(B$loadings), PC2_periodontal=round(B$loadings[,1],4),
            PC3_commensal=round(B$loadings[,2],4)),
            file.path(OUT,"reference","axis_loadings.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

rv  <- which(sm$cohort_id=="RELVAS21_EU")
lab <- ifelse(sm$class[rv]=="train","healthy", sub(".*Grade_","grade_",sm$label[rv]))
ctl <- rv[lab=="healthy"]; pat <- rv[lab!="healthy"]
wr <- function(idx, f) { D <- as.data.frame(M[idx,,drop=FALSE])
  D <- cbind(sample_id=rownames(M)[idx], D)
  write.table(D, file.path(OUT,"examples",f), sep="\t", quote=FALSE, row.names=FALSE) }
wr(ctl, "example_controls.tsv"); wr(pat, "example_patients.tsv")
write.table(data.frame(sample_id=rownames(M)[pat], clinical_grade=lab[lab!="healthy"],
            age=sm$age[pat]), file.path(OUT,"examples","example_patients_truth.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

writeLines(c(
'# Generates a one-page HTML report. Usage:',
'#   source("scripts/score.R"); source("scripts/report_html.R")',
'#   ref <- lab_reference(controls); report_html(patient_counts, ref, "report.html", "ID-001")',
'report_html <- function(cts, ref, file="report.html", sample_id="") {',
'  r <- report(cts, ref)',
'  col <- c(LOW="#2E7D32", INTERMEDIATE="#E8871A", ELEVATED="#C62828")[r$band]',
'  bar <- function(v) sprintf("<div style=\'background:#eee;height:10px;border-radius:5px\'>',
'    <div style=\'background:%s;height:10px;width:%.0f%%;border-radius:5px\'></div></div>", col, 100*v)',
'  tb <- function(x, up) paste(sprintf(',
'    "<tr><td style=\'padding:4px 10px\'><i>%s</i></td><td style=\'padding:4px 10px;text-align:right;color:%s\'>%+.2f</td></tr>",',
'    names(x), ifelse(up,"#C62828","#2E7D32"), x), collapse="")',
'  writeLines(sprintf("<!doctype html><html><head><meta charset=\'utf-8\'>',
'  <title>Oral-health report</title></head>',
'  <body style=\'font-family:system-ui,sans-serif;max-width:720px;margin:40px auto;color:#222\'>',
'  <h1 style=\'margin:0;font-size:24px\'>Salivary oral-health report</h1>',
'  <p style=\'color:#666;margin:4px 0 24px\'>Sample %s &middot; %s</p>',
'  <div style=\'border:1px solid #ddd;border-radius:8px;padding:20px;margin-bottom:20px\'>',
'    <div style=\'font-size:13px;color:#666\'>RESULT</div>',
'    <div style=\'font-size:30px;font-weight:600;color:%s;margin:6px 0\'>%s</div>',
'    <div style=\'margin:14px 0 6px;font-size:14px\'>%.0fth percentile of %d healthy controls</div>',
'    %s</div>',
'  <h3 style=\'font-size:15px\'>Community axes</h3>',
'  <table style=\'width:100%%;font-size:14px;border-collapse:collapse\'>',
'  <tr><td style=\'padding:4px 10px\'>Periodontal pathogen axis</td>',
'      <td style=\'padding:4px 10px;text-align:right\'>%+.2f</td></tr>',
'  <tr><td style=\'padding:4px 10px\'>Commensal / aerobe axis</td>',
'      <td style=\'padding:4px 10px;text-align:right\'>%+.2f</td></tr></table>',
'  <h3 style=\'font-size:15px\'>Raising this score</h3>',
'  <table style=\'width:100%%;font-size:14px;border-collapse:collapse\'>%s</table>',
'  <h3 style=\'font-size:15px\'>Depleted relative to reference</h3>',
'  <table style=\'width:100%%;font-size:14px;border-collapse:collapse\'>%s</table>',
'  <p style=\'margin-top:28px;font-size:12px;color:#888;border-top:1px solid #eee;padding-top:12px\'>',
'  Wellness and triage indicator. Not a diagnostic test and not a cancer screen.',
'  Scores are valid only against the reference panel of this laboratory.',
'  Genus coverage %.1f%%.</p></body></html>",',
'  sample_id, format(Sys.Date()), col, r$band, 100*r$percentile, ref$n, bar(r$percentile),',
'  r$axes[1], r$axes[2], tb(r$drivers_up,TRUE), tb(r$drivers_down,FALSE), 100*r$coverage), file)',
'  cat("written:", file, "\\n"); invisible(r) }'),
  file.path(OUT,"scripts","report_html.R"))

writeLines(c(
'# Runs the full example end to end.',
'source("scripts/score.R"); source("scripts/report_html.R")',
'ctl <- read.delim("examples/example_controls.tsv", row.names=1, check.names=FALSE)',
'pat <- read.delim("examples/example_patients.tsv", row.names=1, check.names=FALSE)',
'tru <- read.delim("examples/example_patients_truth.tsv")',
'ref <- lab_reference(as.matrix(ctl))',
'cat(sprintf("reference: %d controls\\n\\n", ref$n))',
's <- apply(as.matrix(pat), 1, function(r) report(r, ref)$percentile)',
'print(round(tapply(s, tru$clinical_grade, median), 2))',
'k <- which(tru$clinical_grade=="grade_3")[1]',
'report_html(as.matrix(pat)[k,], ref, "example_report.html", tru$sample_id[k])'),
  file.path(OUT,"examples","run_example.R"))

writeLines(c(
sprintf("# Salivary Oral-Health Score v%s", V), "",
"A saliva-based triage and wellness score derived from 16S rRNA genus profiles.",
"**Not a diagnostic test. Not a cancer screen.**", "",
"## Quick start", "```r",
'setwd("<this folder>")',
'source("scripts/score.R"); source("scripts/report_html.R")',
'ctl <- read.delim("examples/example_controls.tsv", row.names=1, check.names=FALSE)',
'ref <- lab_reference(as.matrix(ctl))      # YOUR laboratory\'s healthy controls',
'report(patient_counts, ref)               # named vector of genus counts', "```", "",
"## The local calibration requirement", "",
"Scores are **not comparable between laboratories**. Study of origin accounts for",
"27.7% of compositional variance. Each site must build its own reference panel from",
"healthy controls processed identically. Minimum 30 recommended; validation used 17.",
"Against a global reference only 27% of severe cases were flagged; against local",
"controls, 82%.", "",
"## Performance (external cohort, never seen in training)", "",
"| Metric | Value |", "|---|---|",
"| AUC healthy vs moderate-severe | 0.844 |",
"| Sensitivity at 90th-percentile threshold | 70% |",
"| False positives among healthy | 12% |", "",
"Roughly one patient in five bands against their clinical grade. Not a per-patient",
"diagnosis.", "",
"## Input requirements", "",
"- 16S rRNA, any of V1-V2 / V1-V3 / V3-V4 / V4",
"- Genus-level counts classified against eHOMD V15.22",
"- Minimum 3,000 bacterial reads after host removal",
"- Genus coverage >= 80% or the score is withheld", "",
"## Contents", "",
"- `models/oralclock_final.rds` — model, axes, global reference",
"- `scripts/score.R` — scoring engine",
"- `scripts/report_html.R` — one-page HTML report",
"- `reference/genus_template.tsv` — the expected genera",
"- `reference/coefficients.tsv` — model weights and direction",
"- `reference/axis_loadings.tsv` — the two interpretable axes",
"- `examples/` — 17 controls and 60 graded patients with ground truth", "",
"## Citation", "",
"Codoñer FM et al. (in preparation). IATA-CSIC."),
  file.path(OUT,"README.md"))

f <- list.files(OUT, recursive=TRUE, full.names=TRUE)
writeLines(c(sprintf("oral-health-score v%s   built %s", V, Sys.Date()),
  sprintf("R %s | glmnet %s", getRversion(), packageVersion("glmnet")),
  sprintf("model: %d genera, %d selected", length(B$genera), length(cf)), "",
  sprintf("%-52s %10s", "file", "bytes"),
  sprintf("%-52s %10d", sub(paste0(OUT,"/"),"",f), file.size(f))),
  file.path(OUT,"MANIFEST.txt"))
cat(sprintf("\npackage built: %s\n%d files, %.1f MB\n", OUT, length(f), sum(file.size(f))/1e6))
