# --------------------------------------------------------------------------
# fdr_and_kumpitsch.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
## ---- 1. FDR on genus-age associations -------------------------------
G  <- readRDS("features/genus_counts.rds")
SM <- read.delim("features/sample_metadata.tsv", stringsAsFactors=FALSE)
rownames(G) <- unname(rownames(G)); SM <- SM[match(rownames(G), SM$run_accession), ]
H <- !nzchar(SM$label) & !is.na(SM$age); G <- G[H,]; SM <- SM[H,]
clr <- function(M){ M <- M+1; L <- log(M); L - rowMeans(L) }
X <- clr(G); coh <- SM$cohort_id; age <- SM$age

z <- function(r) 0.5*log((1+r)/(1-r))
res <- t(sapply(colnames(X), function(g){
  st <- do.call(rbind, lapply(split(seq_along(age), coh), function(ix){
    if (length(ix) < 20) return(NULL)
    r <- suppressWarnings(cor(X[ix,g], age[ix], method="spearman"))
    if (is.na(r)) return(NULL)
    c(z(max(min(r,0.999),-0.999)), 1/(length(ix)-3))
  }))
  if (is.null(st) || nrow(st) < 5) return(c(NA,NA,NA,NA))
  yi <- st[,1]; vi <- st[,2]
  w <- 1/vi; mu <- sum(w*yi)/sum(w)
  Q <- sum(w*(yi-mu)^2); df <- length(yi)-1
  tau2 <- max(0, (Q-df)/(sum(w)-sum(w^2)/sum(w)))
  w2 <- 1/(vi+tau2); mu2 <- sum(w2*yi)/sum(w2); se <- sqrt(1/sum(w2))
  p <- 2*pnorm(-abs(mu2/se))
  c(tanh(mu2), p, nrow(st), tau2)
}))
colnames(res) <- c("pooled_rho","p","k_cohorts","tau2")
R <- data.frame(genus=rownames(res), res, row.names=NULL)
R <- R[!is.na(R$p), ]
R$q_BH <- p.adjust(R$p, "BH")
R <- R[order(R$q_BH), ]
write.csv(R, "results/genus_age_fdr.csv", row.names=FALSE)
cat(sprintf("genera tested %d | p<0.05: %d | q<0.05: %d | q<0.10: %d\n",
    nrow(R), sum(R$p<0.05), sum(R$q_BH<0.05), sum(R$q_BH<0.10)))
print(head(R[R$q_BH<0.10, c("genus","pooled_rho","p","q_BH","k_cohorts")], 25), row.names=FALSE)

## ---- 2. disease result without Kumpitsch ----------------------------
cat("\n---- disease cohorts, leave-Kumpitsch-out ----\n")
rv <- tryCatch(readRDS("results/reported_values.rds"), error=function(e) NULL)
if (is.null(rv)) { cat("reported_values.rds not found; skipping\n") } else {
  d <- rv$disease
  if (is.null(d)) { cat("no $disease element; names: ", paste(names(rv),collapse=", "), "\n") }
  else { print(d)
    a <- suppressWarnings(as.numeric(d$auc %||% d$AUC))
    nm <- d$cohort %||% rownames(d)
    k <- !grepl("KUMPITSCH", nm, ignore.case=TRUE)
    cat(sprintf("\nmean AUC all %d cohorts: %.3f | without Kumpitsch (%d): %.3f\n",
        length(a), mean(a,na.rm=TRUE), sum(k), mean(a[k],na.rm=TRUE))) }
}
