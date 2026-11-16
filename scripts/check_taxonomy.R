# --------------------------------------------------------------------------
# check_taxonomy.R
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
tx <- readRDS("taxonomy/taxa_oriented.rds")$tax
ph <- is.na(tx[,"Phylum"]); gn <- is.na(tx[,"Genus"])
names(ph) <- names(gn) <- rownames(tx)

f  <- list.files(c("final_split/train","final_split/valid"), pattern="rds$", full.names=TRUE)
b  <- sub("_seqtab.rds$", "", basename(f))
cl <- ifelse(grepl("/train/", f), "T", "V")
o  <- order(sub(".*__","",b), b)

cat(sprintf("%-34s %2s %9s %6s %10s %10s %9s\n",
    "group","cl","reads","modal","unasgPhy%","unasgGen%","ASVunas%"))
res <- list()
for (i in o) {
  st <- readRDS(f[i]); a <- colSums(st); s <- colnames(st)
  md <- as.integer(names(which.max(tapply(a, nchar(s), sum))))
  p <- 100*sum(a[ph[s]])/sum(a); g <- 100*sum(a[gn[s]])/sum(a)
  cat(sprintf("%-34s %2s %9.0f %6d %9.2f%% %9.2f%% %8.1f%%%s\n",
      b[i], cl[i], sum(a), md, p, g, 100*mean(ph[s]),
      if (p > 15) "   <-- GATE" else ""))
  res[[b[i]]] <- c(p=p, g=g)
}
m <- do.call(rbind, res)
cat(sprintf("\nread-weighted unassigned at phylum: median %.2f%%  max %.2f%%\n",
    median(m[,"p"]), max(m[,"p"])))
cat(sprintf("read-weighted unassigned at genus : median %.2f%%  max %.2f%%\n",
    median(m[,"g"]), max(m[,"g"])))
cat("\n--- is unassignment a length artefact? ---\n")
asv <- rownames(tx); L <- nchar(asv)
print(round(tapply(ph, cut(L, c(0,180,220,260,320,380,440,600)), mean)*100, 1))
