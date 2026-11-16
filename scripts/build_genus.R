# --------------------------------------------------------------------------
# build_genus.R
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
MIN <- 3000
dir.create("features", showWarnings=FALSE)

lab <- tx[,"Genus"]
for (r in c("Family","Order","Class","Phylum")) {
  i <- is.na(lab) & !is.na(tx[,r])
  lab[i] <- paste0(tx[i,r], " (unclassified)")
}
keep <- !is.na(tx[,"Phylum"]) & !is.na(tx[,"Kingdom"]) & tx[,"Kingdom"]=="Bacteria"
names(lab) <- names(keep) <- rownames(tx)
cat(sprintf("bacterial ASVs %d / %d | distinct genus labels %d\n",
    sum(keep), nrow(tx), length(unique(lab[keep]))))

FWDONLY <- c("PORTUGAL_74247__paired__V1-V2","SPAIN_74796__paired__V1-V2",
             "CHAUDHARI20_IN__paired__V3-V4","ANTARCTIC_19274__paired__V4")
part <- read.delim("metadata/partition.tsv", stringsAsFactors=FALSE)
meta <- read.delim("metadata/FINAL_sample_set_M1.tsv", stringsAsFactors=FALSE,
                   quote="", comment.char="")
age <- setNames(meta$age_years, meta$run_accession)
sex <- setNames(meta$sex,       meta$run_accession)

f <- list.files(c("final_split/train","final_split/valid"), pattern="rds$", full.names=TRUE)
tabs <- list(); rows <- list()
for (i in seq_along(f)) {
  g  <- sub("_seqtab.rds$","",basename(f[i]))
  st <- readRDS(f[i]); s <- colnames(st)
  st <- st[, keep[s], drop=FALSE]
  bd <- rowSums(st)
  st <- st[bd >= MIN, , drop=FALSE]
  if (!nrow(st)) next
  cs <- t(rowsum(t(st), lab[colnames(st)]))
  key <- paste(if (grepl("/train/", f[i])) "T" else "V", g, sep="|")
  tabs[[key]] <- cs
  rte <- if (g %in% FWDONLY) "single" else strsplit(g,"__")[[1]][2]
  rows[[key]] <- data.frame(run_accession=rownames(cs), group=g,
    cohort_id=strsplit(g,"__")[[1]][1], region=sub(".*__","",g), route=rte,
    bact_depth=rowSums(cs), stringsAsFactors=FALSE)
}
gen <- sort(unique(unlist(lapply(tabs, colnames))))
M <- matrix(0L, nrow=sum(sapply(tabs,nrow)), ncol=length(gen),
            dimnames=list(unlist(lapply(tabs,rownames)), gen))
r0 <- 0L
for (cs in tabs) { M[(r0+1):(r0+nrow(cs)), colnames(cs)] <- cs; r0 <- r0+nrow(cs) }

sm <- do.call(rbind, rows); rownames(sm) <- NULL
sm <- merge(sm, part[,c("run_accession","class","label")], by="run_accession", all.x=TRUE)
sm$age <- age[sm$run_accession]; sm$sex <- sex[sm$run_accession]
sm <- sm[match(rownames(M), sm$run_accession), ]

R <- M / rowSums(M)
saveRDS(M, "features/genus_counts.rds"); saveRDS(R, "features/genus_relab.rds")
write.table(sm, "features/sample_metadata.tsv", sep="\t", quote=FALSE, row.names=FALSE)

cat(sprintf("\nmatrix: %d samples x %d genera | sparsity %.1f%% zeros\n",
    nrow(M), ncol(M), 100*mean(M==0)))
cat(sprintf("train %d | valid %d | cohorts %d | ages %d-%d (median %d)\n",
    sum(sm$class=="train"), sum(sm$class=="validate"),
    length(unique(sm$cohort_id)), min(sm$age,na.rm=TRUE), max(sm$age,na.rm=TRUE),
    median(sm$age,na.rm=TRUE)))
cat("missing age:", sum(is.na(sm$age)), "| missing sex:", sum(is.na(sm$sex)|sm$sex==""), "\n")
cat("\nprevalence: genera in >10% of samples:", sum(colMeans(M>0)>0.10),
    "| >50%:", sum(colMeans(M>0)>0.50), "\n")
tg <- sort(colMeans(R), decreasing=TRUE)[1:15]
cat("\ntop 15 genera by mean relative abundance:\n")
print(round(100*tg, 2))
