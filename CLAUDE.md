# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-only bioinformatics project studying translational regulation by eIF3d and eIF3e during hypoxia in MDA-MB-231 breast cancer cells. It integrates ribosome profiling + RNA-seq to compute translation efficiency (TE), then builds Random Forest models to predict which genes are translationally regulated by eIF3d/3e and what sequence/structural features drive that regulation.

## Running Code

All analysis lives in R Markdown notebooks (`.Rmd`). Open `hypoxia_timecourse.Rproj` in RStudio and knit notebooks individually. There is no Makefile or pipeline runner — notebooks are run in the numerical order implied by their filenames.

To run a notebook non-interactively:
```r
rmarkdown::render("code/predictive_modeling/07_rf_classification_eif3d_targets_no_boruta.Rmd")
```

The one-time setup script `code/predictive_modeling/00_build_txdb_v49.R` builds the GENCODE v49 transcript database (`accessories/human/txdb.gencode49.sqlite`). Run it once before any feature extraction.

No renv or conda — packages must be installed manually. Key packages: `DESeq2`, `tidyverse`, `here`, `GenomicFeatures`, `BSgenome.Hsapiens.UCSC.hg38`, `Biostrings`, `coRdon`, `randomForest`, `ranger`, `caret`, `clusterProfiler`, `msigdbr`, `ggrepel`, `pheatmap`.

## Repository Structure

```
code/                        # All analysis notebooks
  functions.R                # Shared utility functions — check here before writing new helpers
  predictive_modeling/       # Feature extraction + ML pipeline (numbered 01–40)
output/                      # Generated CSVs and RDS files
  predictive_modeling/       # Feature matrices and RF model outputs
  genesets/                  # Gene signature CSVs (e.g. hypoxia_3d_promotes_TE_1hr_lfc0.5.csv)
accessories/
  human/                     # GENCODE v49 GTF, txdb SQLite, gene annotation table
  translation_signatures_literature.csv
  csc_data/                  # CLIP-seq datasets
plots/                       # All generated PDFs (~01j_*, rscu_*, etc.)
counts/                      # Raw RNA-seq / ribosome profiling read counts
cluster_scripts/             # ViennaRNA scripts for secondary structure computation
```

## Predictive Modeling Pipeline

The `code/predictive_modeling/` notebooks form a sequential pipeline:

**Feature extraction (01_*.Rmd)** — each notebook reads the previous notebook's output RDS, appends new columns, and writes a new RDS to `output/predictive_modeling/`:

| Notebook | Output RDS | Features added |
|----------|-----------|----------------|
| `01_` | `feature_matrix_1hr.rds` | UTR/CDS lengths, GC%, k-mers, RBP motifs (~900 features) |
| `01b_` | `feature_matrix_1hr_targeted.rds` | eIF3 CLIP binding (7 conditions) |
| `01c_` | `feature_matrix_1hr_combined.rds` | Combined all above (~6,300 features) |
| `01d_` | `feature_matrix_1hr_structure.rds` | RNA secondary structure (RNAplfold/RNAfold) |
| `01e_` | `feature_matrix_1hr_clip_split.rds` | Split CLIP into `clip_eif3_count` / `clip_cap_binding_count` |
| `01f_` | `feature_matrix_te_lfc.rds` | TE features per condition (8 columns) |
| `01g_` | `feature_matrix_codon_optimality.rds` | CSC, GC3 codon optimality |
| `01h_` | `feature_matrix_kozak.rds` | Kozak sequence features |
| `01i_` | `feature_matrix_external_stability.rds` | External mRNA stability features |
| `01j_` | `feature_matrix_csc_positional.rds` | Positional CSC (full CDS + quarters + first 75 codons) |
| `01k_` | (appends to csc_positional) | CDS/3'UTR structure features |
| `01l_` | `feature_matrix_codon_features.rds` | bg-corrected RSCU + raw codon frequencies (121 features) |

**Negative control definition:**
- `02_negative_control_set.Rmd` → `output/predictive_modeling/negative_control_genes.csv`
- Negative controls = genes with |TE_LFC| < 0.5 AND padj > 0.3 across ALL si3d conditions

**RF classification models (03–10):**
- `07_rf_classification_eif3d_targets_no_boruta.Rmd` — primary model, run across multiple conditions
- `09_cross_condition_importance_comparison.Rmd` — compare feature importance across conditions

**Analysis notebooks (11–40):** Feature-specific deep dives, external validations, regression models.

## Key Data Objects and Conventions

**Gene identifiers** — three related columns always present together in feature matrices and output tables:
- `gene_id` — versioned ENSEMBL ID (e.g. `ENSG00000000003.14`)
- `gene_id_clean` — version-stripped: `sub("\\..*", "", gene_id)`
- `symbol` — HGNC gene symbol
- `transcript_id_clean` — version-stripped transcript ID (primary key in feature matrices)

**Feature matrices** are RDS files with one row per gene/transcript, keyed by `transcript_id_clean`. List columns may be present (e.g. nested codon data from `coRdon`); coerce before using dplyr aggregation: `as.integer(unlist(col))`.

**Translation efficiency (TE)** = log2(ribosome footprint / RNA) computed under si3d vs sictrl knockdown. Positive TE LFC = eIF3d promotes translation; negative = inhibits.

**Gene group conventions** (used across all visualization notebooks):
```r
group_colors <- c(
  "eIF3d hypoxia promotes" = "#E69F00",
  "Negative control"       = "grey40",
  "All other genes"        = "grey80"
)
```
Positive gene set: `output/genesets/hypoxia_3d_promotes_TE_1hr_lfc0.5.csv` → column `ensembl_gene`
Negative controls: `output/predictive_modeling/negative_control_genes.csv`, size-matched with `set.seed(9)`.

**Plot conventions:**
- Save to `plots/` via `here("plots", "notebooknum_description.pdf")`
- `theme_classic(base_size = 12)`, density `alpha = 0.15`, `linewidth = 0.8`
- Wilcoxon tests annotated with `x = Inf, y = Inf, hjust = 1.1, vjust = 1.4`

## Reference Files

| File | Contents |
|------|---------|
| `accessories/human/txdb.gencode49.sqlite` | GENCODE v49 transcript database (built by `00_build_txdb_v49.R`) |
| `accessories/human/gene_anno_hs_dm_v49_r111.tsv` | Gene annotation: gene_id, symbol, biotype |
| `accessories/translation_signatures_literature.csv` | Published translation signature gene sets |
| `output/translation_categories_si3d_vs_sictrl_normoxia_1and4hr.csv` | Full TE differential results |
| `output/predictive_modeling/negative_control_genes.csv` | Negative control gene pool |
