# --------------------------------------------------------------------------
# diag_iberian.R
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
ph <- is.na(tx[,"Phylum"]); names(ph) <- rownames(tx)
gen <- tx[,"Genus"]; names(gen) <- rownames(tx)

for (g in c("SPAIN_74796__paired__V1-V2","PORTUGAL_74247__paired__V1-V2")) {
  for (d in c("train","valid")) {
    p <- sprintf("final_split/%s/%s_seqtab.rds", d, g)
    if (!file.exists(p)) next
    st <- readRDS(p); a <- colSums(st); s <- colnames(st)
    u <- s[ph[s]]; au <- a[u]
    cat(sprintf("\n===== %s [%s] =====\n", g, d))
    cat(sprintf("reads %.0f | ASVs %d | unassigned ASVs %d carrying %.1f%% of reads\n",
        sum(a), length(a), length(u), 100*sum(au)/sum(a)))
    o <- order(au, decreasing=TRUE)
    cat(sprintf("top 10 unassigned carry %.1f%% of ALL reads; top 50 carry %.1f%%\n",
        100*sum(au[o][1:min(10,length(o))])/sum(a),
        100*sum(au[o][1:min(50,length(o))])/sum(a)))
    cat("\n  rank  reads   %all  first40bp\n")
    for (i in seq_len(min(8, length(o)))) {
      cat(sprintf("  %4d %7.0f %6.2f  %s\n", i, au[o][i], 100*au[o][i]/sum(a),
          substr(u[o][i], 1, 40)))
    }
    pre <- table(substr(u, 1, 12))
    cat("\n  commonest 12bp prefixes among UNASSIGNED:\n")
    print(head(sort(pre, decreasing=TRUE), 5))
    pa <- table(substr(s[!ph[s]], 1, 12))
    cat("  commonest 12bp prefixes among ASSIGNED:\n")
    print(head(sort(pa, decreasing=TRUE), 5))
    cat("  top assigned genera:", paste(head(names(sort(tapply(a[!ph[s]],
        gen[s[!ph[s]]], sum), decreasing=TRUE)), 6), collapse=", "), "\n")
  }
}
