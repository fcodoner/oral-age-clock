# --------------------------------------------------------------------------
# boot_aligned.R
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
src <- readLines("scripts/reported_values.R")
i <- grep("^XG *<-", src)[1]
stopifnot(!is.na(i))
cat("evaluating reported_values.R lines 1 to", i, "\n")
eval(parse(text=paste(src[1:i], collapse="\n")), envir=globalenv())
cat("XG:", nrow(XG), "samples x", ncol(XG), "genera | K =", length(K), "\n")

SM <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
SM <- SM[match(rownames(XG), SM$run_accession), ]
stopifnot(identical(SM$run_accession, rownames(XG)))
H  <- !nzchar(SM$label) & !is.na(SM$age)
X  <- XG[H,,drop=FALSE]; y <- SM$age[H]; coh <- SM$cohort_id[H]
cat("healthy:", nrow(X), "| cohorts:", length(unique(coh)), "\n")

FF <- function(n) rep(1:5, length.out=n)
loco <- function(X,y,f){ p <- rep(NA_real_, length(y))
  for (k in unique(f)) { te <- f==k
    cv <- cv.glmnet(X[!te,,drop=FALSE], y[!te], alpha=0.5, foldid=FF(sum(!te)))
    p[te] <- as.numeric(predict(cv, X[te,,drop=FALSE], s="lambda.min")) }
  p }
mae <- function(a,b) mean(abs(a-b))

p_loco <- loco(X,y,coh)
set.seed(1); p_rand <- loco(X,y,sample(FF(length(y))))
cat(sprintf("\nLOCO      MAE %.3f  null %.3f  r %.3f\n", mae(p_loco,y), mae(mean(y),y), cor(p_loco,y)))
cat(sprintf("randomCV  MAE %.3f           r %.3f\n", mae(p_rand,y), cor(p_rand,y)))

CU <- unique(coh); set.seed(7); N <- 2000
bl<-bn<-br<-numeric(N)
for (b in 1:N){ s <- sample(CU,length(CU),replace=TRUE)
  ix <- unlist(lapply(s, function(c) which(coh==c)))
  bl[b]<-mae(p_loco[ix],y[ix]); bn[b]<-mae(mean(y[ix]),y[ix]); br[b]<-mae(p_rand[ix],y[ix]) }
q <- function(v) sprintf("%.2f [%.2f, %.2f]", mean(v), quantile(v,.025), quantile(v,.975))
cat("\n-- cluster bootstrap over cohorts --\n")
cat("LOCO MAE           ", q(bl), "\n")
cat("improvement vs null", q(bn-bl), "\n")
cat("LOCO - randomCV    ", q(bl-br), "\n")

cat("\n-- leave-TRIAD-out --\n")
k <- coh!="TRIAD_ES"
if (sum(!k)) { pt <- loco(X[k,,drop=FALSE], y[k], coh[k])
  cat(sprintf("without TRIAD  MAE %.3f  r %.3f  n=%d  (with: %.3f / %.3f)\n",
      mae(pt,y[k]), cor(pt,y[k]), sum(k), mae(p_loco,y), cor(p_loco,y)))
  wc <- function(p,yy,cc) median(sapply(split(seq_along(yy),cc),
        function(i) if(length(i)>=20) cor(p[i],yy[i]) else NA), na.rm=TRUE)
  cat(sprintf("median within-cohort r: with %.3f | without %.3f\n",
      wc(p_loco,y,coh), wc(pt,y[k],coh[k]))) } else cat("TRIAD_ES absent\n")

saveRDS(list(loco=p_loco, rand=p_rand, boot=list(l=bl,n=bn,r=br)), "results/boot_aligned.rds")
cat("\nsaved results/boot_aligned.rds\n")
