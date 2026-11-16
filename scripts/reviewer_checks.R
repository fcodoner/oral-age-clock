# --------------------------------------------------------------------------
# reviewer_checks.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
suppressMessages({library(glmnet)})
G  <- readRDS("features/genus_counts.rds")
SM <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
rownames(G) <- unname(rownames(G))
SM <- SM[match(rownames(G), SM$run_accession), ]
stopifnot(identical(SM$run_accession, rownames(G)))

H <- !nzchar(SM$label) & !is.na(SM$age)
G <- G[H,]; SM <- SM[H,]
cat("healthy samples:", nrow(G), "| cohorts:", length(unique(SM$cohort_id)), "\n")

clr <- function(M){ M <- M+1; L <- log(M); L - rowMeans(L) }
X <- clr(G); y <- SM$age; coh <- SM$cohort_id
FF <- function(n) rep(1:5, length.out=n)

loco <- function(X, y, folds){
  p <- rep(NA_real_, length(y))
  for (f in unique(folds)) {
    te <- folds==f; tr <- !te
    fit <- cv.glmnet(X[tr,,drop=FALSE], y[tr], alpha=0.5, foldid=FF(sum(tr)))
    p[te] <- as.numeric(predict(fit, X[te,,drop=FALSE], s="lambda.min"))
  }
  p
}
mae <- function(a,b) mean(abs(a-b))

cat("\n== 1. baseline LOCO ==\n")
p_base <- loco(X, y, coh)
cat(sprintf("MAE %.3f | null %.3f | r %.3f\n", mae(p_base,y), mae(mean(y),y), cor(p_base,y)))

cat("\n== 2. random 5-fold CV on the same data ==\n")
set.seed(1); rf <- sample(FF(length(y)))
p_rand <- loco(X, y, rf)
cat(sprintf("MAE %.3f | r %.3f\n", mae(p_rand,y), cor(p_rand,y)))

cat("\n== 3. cluster bootstrap over cohorts (2000 draws) ==\n")
CU <- unique(coh); set.seed(7)
bl <- bn <- br <- bg <- numeric(2000)
for (i in 1:2000) {
  s  <- sample(CU, length(CU), replace=TRUE)
  ix <- unlist(lapply(s, function(c) which(coh==c)))
  bl[i] <- mae(p_base[ix], y[ix])
  bn[i] <- mae(mean(y[ix]), y[ix])
  br[i] <- mae(p_rand[ix], y[ix])
  bg[i] <- bl[i] - br[i]
}
q <- function(v) sprintf("%.2f [%.2f, %.2f]", mean(v), quantile(v,.025), quantile(v,.975))
cat("LOCO MAE          ", q(bl), "\n")
cat("null MAE          ", q(bn), "\n")
cat("LOCO improvement  ", q(bn-bl), "  (above 0 = better than null)\n")
cat("randomCV MAE      ", q(br), "\n")
cat("LOCO - randomCV   ", q(bg), "  (above 0 = LOCO worse)\n")

cat("\n== 4. sex as covariate ==\n")
hs <- SM$sex %in% c("male","female")
cat("samples with sex:", sum(hs), "of", nrow(SM), "\n")
Xs <- cbind(X[hs,,drop=FALSE], sexM = as.numeric(SM$sex[hs]=="male"))
p_sex <- loco(Xs, y[hs], coh[hs])
p_nos <- loco(X[hs,,drop=FALSE], y[hs], coh[hs])
cat(sprintf("without sex  MAE %.3f  r %.3f\n", mae(p_nos,y[hs]), cor(p_nos,y[hs])))
cat(sprintf("with sex     MAE %.3f  r %.3f\n", mae(p_sex,y[hs]), cor(p_sex,y[hs])))

cat("\n== 5. UIC cohorts pooled into one fold ==\n")
c2 <- coh; c2[c2 %in% c("USA_39492","USA_74379")] <- "UIC_POOLED"
if (sum(c2=="UIC_POOLED")>0) {
  p_uic <- loco(X, y, c2)
  cat(sprintf("UIC pooled   MAE %.3f  r %.3f   (baseline %.3f / %.3f)\n",
      mae(p_uic,y), cor(p_uic,y), mae(p_base,y), cor(p_base,y)))
} else cat("UIC cohorts not present among healthy samples\n")

cat("\n== 6. AGP tongue swabs excluded ==\n")
k <- coh!="AGP_ORAL"
p_nag <- loco(X[k,,drop=FALSE], y[k], coh[k])
cat(sprintf("without AGP  MAE %.3f  r %.3f   n=%d\n", mae(p_nag,y[k]), cor(p_nag,y[k]), sum(k)))

saveRDS(list(p_base=p_base,p_rand=p_rand,boot=list(loco=bl,null=bn,rand=br,gap=bg)),
        "results/reviewer_checks.rds")
cat("\nsaved results/reviewer_checks.rds\n")
