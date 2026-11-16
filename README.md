# Oral microbiome age clock: cross-study validation of salivary 16S models

Analysis code and derived data for the manuscript *Oral disease signatures
transfer between studies while chronological-age prediction is attenuated:
4,426 public salivary 16S rRNA samples*.

Francisco M. Codoner, Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
ORCID [0000-0003-3671-8148](https://orcid.org/0000-0003-3671-8148)

## What is here

```
scripts/      analysis code
features/     genus count and relative-abundance matrices, sample metadata
results/      every object from which a reported value is taken
figures/      the figures as submitted
metadata/     curated per-sample and per-cohort tables; the two largest are gzipped and are read directly by read.delim() in R
env/          sessionInfo and external tool versions
```

Raw sequences are **not** deposited here. All 45 BioProject accessions are
listed in Supplementary Table S4 and are available from ENA, SRA and DDBJ.
`features/` contains the small derived matrices, so the analyses below can be
re-run without reprocessing the raw data.

## Reproducing the reported values

Requires R >= 4.5 with `glmnet`, `vegan`, `ranger`. From the repository root:

```r
Rscript scripts/reported_values.R      # Tables 1 and 2, most Results values
Rscript scripts/matched_comparison.R   # LOCO vs random cross-validation
Rscript scripts/permanova.R            # variance partition
Rscript scripts/age_consistency_full.R # meta-analysis across cohorts
Rscript scripts/revision_round1.R      # leakage, contaminants, disease AUC CIs
Rscript scripts/revision_round2.R      # sex, laboratory, cohort exclusions
Rscript scripts/revision_round3.R      # specimen type
Rscript scripts/export_S6.R            # Supplementary Table S6
```

Fold assignments are fixed rather than drawn at random, so reported statistics
are exactly reproducible.

## Which script produces what

| Script | Produces |
|---|---|
| `reported_values.R` | Table 1, Table 2 and most Results values -> results/reported_values.rds |
| `matched_comparison.R` | matched LOCO vs random CV comparison -> results/matched_comparison.rds |
| `permanova.R` | variance partition (Fig. 1a) -> results/permanova.rds |
| `age_consistency_full.R` | genus and pathway meta-analysis -> results/meta_detail.rds |
| `revision_round1.R` | leakage, contaminant and disease-AUC sensitivity -> results/revision_round1.rds |
| `revision_round2.R` | sex, laboratory, cohort-exclusion sensitivity (Supplementary Results 2) |
| `revision_round3.R` | specimen-type stratification -> results/revision_round3.rds |
| `export_S6.R` | Supplementary Table S6 -> results/S6_meta.csv |
| `fig2.R` | Figure 1 |
| `fig56.R` | Figure 3 and Figure 4 |
| `figures.R` | Figure 2, Supplementary Figures 1 and 2 |
| `run_dada2.R` | DADA2 denoising per processing group |
| `assign_taxonomy.R` | eHOMD classification |
| `build_genus.R` | genus-level feature construction |
| `prep_picrust.R` | PICRUSt2 input preparation |

Other scripts in `scripts/` document the data assembly, quality control and
exploratory analyses carried out during the study; the table above is the
canonical path from raw features to the values reported in the manuscript.

## Notes on reproducibility

- `features/sample_metadata.tsv` defines the healthy and disease sets, the
  train/validate split and the cohort identifiers used throughout.
- `results/permanova.rds` holds the variance partition read directly by
  `fig2.R`; the figure contains no hand-entered values.
- `results/meta_detail.rds` carries both `feature` (the 58-character label used
  in figures) and `feature_full` (the complete MetaCyc descriptor). Use
  `feature_full` for any lookup: truncation makes a few distinct pathways share
  a label.

## License

MIT. See `LICENSE`.
