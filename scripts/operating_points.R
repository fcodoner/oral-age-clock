# --------------------------------------------------------------------------
# operating_points.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
b <- readRDS("models/oralclock_v1.rds")
pv <- function(sens, spec, prev) c(
  PPV = sens*prev / (sens*prev + (1-spec)*(1-prev)),
  NPV = spec*(1-prev) / (spec*(1-prev) + (1-sens)*prev))

cat("=== ORAL HEALTH / periodontitis ===\n")
op <- data.frame(thr=c(0.30,0.40,0.50,0.60,0.70),
                 sens=c(0.959,0.905,0.824,0.635,0.486),
                 spec=c(0.280,0.460,0.560,0.660,0.860))
for (prev in c(0.60,0.45,0.30,0.15)) {
  cat(sprintf("\n-- clinic prevalence %.0f%% --\n%8s %7s %7s %7s %7s\n",
      100*prev,"thresh","sens","spec","PPV","NPV"))
  for (i in 1:nrow(op)) { r <- pv(op$sens[i], op$spec[i], prev)
    cat(sprintf("%8.2f %7.3f %7.3f %7.3f %7.3f\n", op$thr[i], op$sens[i], op$spec[i], r[1], r[2])) }
}

cat("\n\n=== OSCC ===\n")
oc <- data.frame(thr=c(0.30,0.50,0.70),
                 sens=c(0.987,0.735,0.329),
                 spec=c(0.130,0.668,0.880))
for (prev in c(0.56,0.05,0.001)) {
  lab <- if (prev>0.5) "study design (case-control)" else
         if (prev>0.01) "high-risk referral clinic" else "general-population screening"
  cat(sprintf("\n-- prevalence %.1f%%  (%s) --\n%8s %7s %7s %9s %7s\n",
      100*prev, lab, "thresh","sens","spec","PPV","NPV"))
  for (i in 1:nrow(oc)) { r <- pv(oc$sens[i], oc$spec[i], prev)
    cat(sprintf("%8.2f %7.3f %7.3f %9.4f %7.4f\n", oc$thr[i], oc$sens[i], oc$spec[i], r[1], r[2])) }
}
