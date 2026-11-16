# --------------------------------------------------------------------------
# revision_round1.R
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
# revision_round1.R
# Sensitivity and robustness analyses for the first revision round.
# F. M. Codoner
# ---------------------------------------------------------------------------
suppressMessages({library(glmnet); library(vegan)})

src <- readLines("scripts/reported_values.R")
i <- grep("^folds *<-", src)[1]
stopifnot(!is.na(i))
eval(parse(text = paste(src[1:i], collapse = "\n")), envir = globalenv())

mae <- function(p, a) mean(abs(p - a))
auc <- function(s, d) {
  r <- rank(s); n1 <- sum(d == 1); n0 <- sum(d == 0)
  if (n1 < 2 || n0 < 2) return(NA_real_)
  (sum(r[d == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
wc <- function(p, yy, cc)
  sapply(split(seq_along(yy), cc),
         function(ii) if (length(ii) >= 20) cor(p[ii], yy[ii]) else NA_real_)

cat("###########################################################\n")
cat("# A. TRAINING-SET ACCOUNTING AND SIZE-MATCHED COMPARISON  #\n")
cat("###########################################################\n\n")

keep <- co %in% folds
Xs <- XG[tr, ][keep, , drop = FALSE]; ys <- y[keep]; cs <- co[keep]
cat(sprintf("full training set : %d samples, %d cohorts\n", length(y), length(unique(co))))
cat(sprintf("matched set       : %d samples, %d cohorts\n\n", nrow(Xs), length(unique(cs))))

tz_full <- sapply(folds, function(f) sum(co != f))
tz_mtch <- sapply(folds, function(f) sum(cs != f))
cat("training size per fold:\n")
cat(sprintf("  primary LOCO  (trains on all training cohorts) : median %d, range %d-%d\n",
            as.integer(median(tz_full)), min(tz_full), max(tz_full)))
cat(sprintf("  matched LOCO  (trains on eligible cohorts only): median %d, range %d-%d\n",
            as.integer(median(tz_mtch)), min(tz_mtch), max(tz_mtch)))
cat(sprintf("  random  5-fold on matched set                  : %d\n", as.integer(round(nrow(Xs) * 4 / 5))))
cat(sprintf("  random 18-fold on matched set                  : %d\n\n", as.integer(round(nrow(Xs) * 17 / 18))))

cat("per-fold detail (matched set):\n")
print(data.frame(cohort = folds, n_test = as.integer(table(cs)[folds]),
                 n_train_matched = as.integer(tz_mtch),
                 n_train_primary = as.integer(tz_full)), row.names = FALSE)

PL <- rep(NA_real_, length(ys))
for (f in folds) {
  te <- cs == f
  cv <- fit(Xs[!te, , drop = FALSE], ys[!te])
  PL[te] <- as.numeric(predict(cv, Xs[te, , drop = FALSE], s = "lambda.min"))
}
rcv <- function(k) {
  set.seed(1); kk <- sample(rep(1:k, length.out = length(ys))); P <- rep(NA_real_, length(ys))
  for (j in 1:k) {
    cv <- fit(Xs[kk != j, , drop = FALSE], ys[kk != j])
    P[kk == j] <- as.numeric(predict(cv, Xs[kk == j, , drop = FALSE], s = "lambda.min"))
  }
  P
}
P5 <- rcv(5); P18 <- rcv(18)

cat(sprintf("\nmatched LOCO     MAE %.3f  r %.3f\n", mae(PL, ys), cor(PL, ys)))
cat(sprintf("random  5-fold   MAE %.3f  r %.3f   gap %+.3f\n",
            mae(P5, ys),  cor(P5, ys),  mae(PL, ys) - mae(P5, ys)))
cat(sprintf("random 18-fold   MAE %.3f  r %.3f   gap %+.3f   <- size-matched\n",
            mae(P18, ys), cor(P18, ys), mae(PL, ys) - mae(P18, ys)))

CU <- unique(cs); set.seed(7); N <- 2000; g5 <- g18 <- numeric(N)
for (b in 1:N) {
  s <- sample(CU, length(CU), replace = TRUE)
  ix <- unlist(lapply(s, function(k) which(cs == k)))
  g5[b]  <- mae(PL[ix], ys[ix]) - mae(P5[ix],  ys[ix])
  g18[b] <- mae(PL[ix], ys[ix]) - mae(P18[ix], ys[ix])
}
qq <- function(v) sprintf("%.2f [%.2f, %.2f]", mean(v), quantile(v, .025), quantile(v, .975))
cat("\npaired cluster bootstrap over cohorts (2000 draws):\n")
cat("  LOCO - random  5-fold :", qq(g5),  "\n")
cat("  LOCO - random 18-fold :", qq(g18), "\n")

cat("\n###########################################################\n")
cat("# B. LEAKAGE SENSITIVITY: FILTERS INSIDE TRAINING FOLDS   #\n")
cat("###########################################################\n\n")
cat("primary pipeline computes the 10% prevalence filter and the contaminant\n")
cat("set on all samples; this refits both inside each training fold.\n\n")

Gt <- G[tr, ][keep, , drop = FALSE]
PLin <- rep(NA_real_, length(ys)); nk <- integer(0)
for (f in folds) {
  te <- cs == f
  Kin <- setdiff(colnames(Gt)[colMeans(Gt[!te, , drop = FALSE] > 0) > 0.10], CT)
  nk <- c(nk, length(Kin))
  Xi <- clr(Gt[, Kin, drop = FALSE])
  cv <- fit(Xi[!te, , drop = FALSE], ys[!te])
  PLin[te] <- as.numeric(predict(cv, Xi[te, , drop = FALSE], s = "lambda.min"))
}
cat(sprintf("genera retained per fold: median %d, range %d-%d (primary: %d)\n\n",
            as.integer(median(nk)), min(nk), max(nk), ncol(Xs)))
cat(sprintf("filter on all samples    MAE %.3f  r %.3f  medWithinR %.3f\n",
            mae(PL, ys),   cor(PL, ys),   median(wc(PL, ys, cs), na.rm = TRUE)))
cat(sprintf("filter inside folds      MAE %.3f  r %.3f  medWithinR %.3f\n",
            mae(PLin, ys), cor(PLin, ys), median(wc(PLin, ys, cs), na.rm = TRUE)))

cat("\n###########################################################\n")
cat("# C. DISEASE AUC CONFIDENCE INTERVALS AND MATCHED COUNTS  #\n")
cat("###########################################################\n\n")

RC <- readRDS("results/disease_roc.rds")
set.seed(1)
for (nm in names(RC)) {
  s <- RC[[nm]]$score; d <- RC[[nm]]$label; ag <- RC[[nm]]$age
  bs <- replicate(2000, {
    ii <- sample(seq_along(d), replace = TRUE)
    if (length(unique(d[ii])) < 2) NA_real_ else auc(s[ii], d[ii])
  })
  ca <- which(d == 1); ct <- which(d == 0); mn <- c(); ma <- c()
  for (rp in 1:200) {
    pool <- ct; pc <- c()
    for (k in sample(ca)) {
      cand <- pool[abs(ag[pool] - ag[k]) <= 5]
      if (!length(cand)) next
      j <- if (length(cand) == 1) cand else sample(cand, 1)
      pc <- c(pc, k, j); pool <- setdiff(pool, j)
    }
    if (length(pc) >= 10) { mn <- c(mn, length(pc)); ma <- c(ma, auc(s[pc], d[pc])) }
  }
  cat(sprintf("%-28s n=%3d cases=%3d ctrl=%3d\n", nm, length(d), sum(d), sum(d == 0)))
  cat(sprintf("%-28s AUC %.3f [%.3f, %.3f]\n", "", auc(s, d),
              quantile(bs, .025, na.rm = TRUE), quantile(bs, .975, na.rm = TRUE)))
  if (length(ma)) {
    cat(sprintf("%-28s age-matched n=%.0f  AUC %.3f  IQR [%.3f, %.3f]  2.5-97.5%% [%.3f, %.3f]  (%d/200 draws)\n\n",
                "", median(mn), median(ma), quantile(ma, .25), quantile(ma, .75),
                quantile(ma, .025), quantile(ma, .975), length(ma)))
  } else cat(sprintf("%-28s age-matched: too few matchable pairs\n\n", ""))
}

cat("###########################################################\n")
cat("# D. CONTAMINANT-RETAINED SENSITIVITY                     #\n")
cat("###########################################################\n\n")
cat(sprintf("contaminant genera excluded (%d): %s\n\n", length(CT), paste(CT, collapse = ", ")))

Kall <- colnames(G)[colMeans(G > 0) > 0.10]
Xall <- clr(G[, Kall, drop = FALSE])[tr, ][keep, , drop = FALSE]
PC <- rep(NA_real_, length(ys))
for (f in folds) {
  te <- cs == f
  cv <- fit(Xall[!te, , drop = FALSE], ys[!te])
  PC[te] <- as.numeric(predict(cv, Xall[te, , drop = FALSE], s = "lambda.min"))
}
cat(sprintf("contaminants removed  (%3d genera)  MAE %.3f  r %.3f\n", ncol(Xs),     mae(PL, ys), cor(PL, ys)))
cat(sprintf("contaminants retained (%3d genera)  MAE %.3f  r %.3f\n", length(Kall), mae(PC, ys), cor(PC, ys)))

cat("\n###########################################################\n")
cat("# E. ORAL-HEALTH GRADIENT: THRESHOLD AND SAMPLE PROVENANCE#\n")
cat("###########################################################\n\n")

GD <- readRDS("results/gradient_detail.rds")
sc <- GD$score; gr <- GD$grade
cat("stored object fields:", paste(names(GD), collapse = ", "), "\n")
cat("grade levels:", paste(levels(gr), collapse = " | "), "\n")
cat("grade counts:", paste(sprintf("%s=%d", levels(gr), as.integer(table(gr))), collapse = "  "), "\n\n")

rv <- which(sm$cohort_id == "RELVAS21_EU")
cat(sprintf("RELVAS21_EU samples: %d\n", length(rv)))
cat("class assignment within RELVAS21_EU:\n"); print(table(sm$class[rv]))
cat("\ncross-tabulation of class against grade label:\n")
print(table(class = sm$class[rv], grade = gr))
cat("\n-> if grade 0 falls entirely in 'train', those scores are IN-SAMPLE\n")
cat("   fitted values while grades 1-3 are out-of-sample predictions.\n\n")

g0lab <- levels(gr)[1]
loc <- sc[gr == g0lab]
thr_in <- quantile(loc, .90)
cat("--- as currently reported (threshold from the same 17 controls) ---\n")
cat(sprintf("threshold %.3f\n", thr_in))
for (g in levels(gr))
  cat(sprintf("  %-9s flagged %5.1f%% (%d/%d)%s\n", g, 100 * mean(sc[gr == g] >= thr_in),
              sum(sc[gr == g] >= thr_in), sum(gr == g),
              if (g == g0lab) "   <- CIRCULAR" else ""))

g0 <- which(gr == g0lab)
loo <- sapply(g0, function(k) sc[k] >= quantile(sc[setdiff(g0, k)], .90))
cat(sprintf("\n--- leave-one-out threshold for grade 0 ---\n  %-9s flagged %5.1f%% (%d/%d)\n",
            g0lab, 100 * mean(loo), sum(loo), length(loo)))

cat("\n--- clean refit: every RELVAS21_EU sample held out of training ---\n")
trn <- setdiff(tr, rv)
cat(sprintf("training n = %d (was %d)\n", length(trn), length(tr)))
cv2 <- fit(XG[trn, , drop = FALSE], sm$age[trn])
s2 <- as.numeric(predict(cv2, XG[rv, , drop = FALSE], s = "lambda.min"))
cat(sprintf("correlation between stored and refit scores: r = %.3f\n", cor(sc, s2)))
loc2 <- s2[gr == g0lab]
g0b <- which(gr == g0lab)
loo2 <- sapply(g0b, function(k) s2[k] >= quantile(s2[setdiff(g0b, k)], .90))
thr2 <- quantile(loc2, .90)
cat(sprintf("\nall 77 scored out-of-sample; threshold %.3f\n", thr2))
for (g in levels(gr))
  cat(sprintf("  %-9s flagged %5.1f%% (%d/%d)\n", g, 100 * mean(s2[gr == g] >= thr2),
              sum(s2[gr == g] >= thr2), sum(gr == g)))
cat(sprintf("  %-9s flagged %5.1f%% under a leave-one-out threshold\n", g0lab, 100 * mean(loo2)))
gn <- as.integer(gr)
cat(sprintf("\nSpearman(score, grade): stored %.3f  clean refit %.3f\n",
            cor(sc, gn, method = "spearman"), cor(s2, gn, method = "spearman")))
cat(sprintf("AUC grade>=1 vs grade 0: stored %.3f  clean refit %.3f\n",
            auc(sc, as.integer(gn > 1)), auc(s2, as.integer(gn > 1))))

cat("\n###########################################################\n")
cat("# F. RECONCILING medWithinR 0.363 AND 0.219              #\n")
cat("###########################################################\n\n")

v1 <- readRDS("results/age_predictions.rds")$loco
cat(sprintf("primary LOCO  (trains on all %d training samples, tests on %d)\n",
            length(y), length(v1$pred)))
cat(sprintf("              MAE %.3f  r %.3f  median within-cohort r %.3f\n",
            mae(v1$pred, v1$actual), cor(v1$pred, v1$actual),
            median(wc(v1$pred, v1$actual, v1$cohort), na.rm = TRUE)))
cat(sprintf("matched LOCO  (trains on the %d eligible samples only, tests on %d)\n",
            nrow(Xs), length(PL)))
cat(sprintf("              MAE %.3f  r %.3f  median within-cohort r %.3f\n",
            mae(PL, ys), cor(PL, ys), median(wc(PL, ys, cs), na.rm = TRUE)))
cat("\nper-cohort within-cohort r, both arms:\n")
a <- wc(v1$pred, v1$actual, v1$cohort); b <- wc(PL, ys, cs)
k <- intersect(names(a), names(b))
print(data.frame(cohort = k, n = as.integer(table(cs)[k]),
                 r_primary = round(a[k], 3), r_matched = round(b[k], 3)), row.names = FALSE)

cat("\n###########################################################\n")
cat("# G. HETEROGENEITY AND DIRECTIONAL AGREEMENT              #\n")
cat("###########################################################\n\n")

M <- readRDS("results/meta_detail.rds")
for (nm in names(M)) {
  d <- M[[nm]]; d <- d[d$q < 0.05, ]
  d <- d[order(d$q), ]
  cat(sprintf("-- %s: %d features at q < 0.05 --\n", nm, nrow(d)))
  print(d, row.names = FALSE); cat("\n")
}

cat("--- recomputed I-squared (method cross-check against stored rho) ---\n")
RM <- if (exists("genus_relab")) genus_relab else G
cat(sprintf("matrix used: %s  (%d x %d)\n", if (exists("genus_relab")) "genus_relab" else "G",
            nrow(RM), ncol(RM)))
dl <- function(z, v) {
  w <- 1 / v; m <- sum(w * z) / sum(w); Q <- sum(w * (z - m)^2); df <- length(z) - 1
  C <- sum(w) - sum(w^2) / sum(w); t2 <- max(0, (Q - df) / C)
  I2 <- if (Q > 0) max(0, 100 * (Q - df) / Q) else 0
  ws <- 1 / (v + t2); ms <- sum(ws * z) / sum(ws)
  c(rho = tanh(ms), tau2 = t2, I2 = I2, k = length(z))
}
GS <- RM[tr, , drop = FALSE]; feats <- M$genus$feature
out <- do.call(rbind, lapply(feats, function(ft) {
  if (!ft %in% colnames(GS)) return(NULL)
  zz <- c(); vv <- c()
  for (f in unique(co)) {
    ii <- which(co == f); if (length(ii) < 20) next
    x <- GS[ii, ft]; if (length(unique(x)) < 3) next
    r <- suppressWarnings(cor(x, y[ii], method = "spearman"))
    if (is.na(r) || abs(r) >= 1) next
    zz <- c(zz, atanh(r)); vv <- c(vv, 1 / (length(ii) - 3))
  }
  if (length(zz) < 5) return(NULL)
  data.frame(feature = ft, t(dl(zz, vv)))
}))
st <- M$genus[match(out$feature, M$genus$feature), c("feature", "rho", "agree", "q")]
cmp <- data.frame(feature = out$feature, rho_stored = round(st$rho, 3),
                  rho_recomputed = round(out$rho, 3),
                  diff = round(out$rho - st$rho, 3),
                  I2 = round(out$I2, 1), tau2 = round(out$tau2, 4),
                  k = out$k, agree = st$agree, q = signif(st$q, 3))
cmp <- cmp[order(cmp$q), ]
cat(sprintf("\nmax |difference| between stored and recomputed rho: %.4f\n", max(abs(cmp$diff))))
cat("if that is near zero the method matches and I2 can be reported.\n\n")
print(head(cmp, 25), row.names = FALSE)
cat(sprintf("\nI2 across the %d significant genera: median %.1f, range %.1f-%.1f\n",
            sum(cmp$q < 0.05), median(cmp$I2[cmp$q < 0.05]),
            min(cmp$I2[cmp$q < 0.05]), max(cmp$I2[cmp$q < 0.05])))

saveRDS(list(matched = list(loco = PL, r5 = P5, r18 = P18, y = ys, cohort = cs),
             fold_internal = PLin, contam_retained = PC,
             gradient_clean = s2, heterogeneity = cmp),
        "results/revision_round1.rds")
cat("\nsaved: results/revision_round1.rds\n")
