# Translational Regulation by eIF3d and eIF3e During Hypoxia

Analysis code for a ribosome profiling + RNA-seq study of translational control by eIF3d and eIF3e in breast cancer cells under normoxia and hypoxia.

## Overview

This project integrates ribosome profiling and RNA-seq data to measure translation efficiency (TE) genome-wide, then applies machine learning to identify sequence and structural features that predict eIF3d/eIF3e-dependent translational regulation. Analyses are performed in MDA-MB-231 and MCF7-SIX1 breast cancer cell lines across a hypoxia timecourse.

**General workflow:**
1. Compute translation efficiency (TE = ribosome footprints / RNA) via DESeq2 under eIF3d or eIF3e knockdown vs. control
2. Extract sequence, structural, and codon-level features for each transcript
3. Train Random Forest classifiers to distinguish eIF3-regulated from unregulated transcripts
4. Interpret feature importance to identify RNA regulatory elements associated with eIF3 sensitivity

## Requirements

- R (≥ 4.2)
- Key packages: `DESeq2`, `tidyverse`, `here`, `GenomicFeatures`, `BSgenome.Hsapiens.UCSC.hg38`, `Biostrings`, `coRdon`, `randomForest`, `ranger`, `caret`, `clusterProfiler`, `msigdbr`, `ggrepel`, `pheatmap`
- GENCODE v49 annotation (GTF) and a pre-built `TxDb` object (see setup below)

Packages must be installed manually — there is no `renv` lockfile or conda environment.

## Setup

Build the transcript database once before running any analysis notebooks:

```r
source("code/predictive_modeling/00_build_txdb_v49.R")
```

This writes `accessories/human/txdb.gencode49.sqlite`, which all downstream notebooks read.

## Running the Analysis

All analyses are R Markdown notebooks (`.Rmd`). Open `hypoxia_timecourse.Rproj` in RStudio and knit notebooks in the numerical order implied by their filenames. Notebooks in `code/predictive_modeling/` form a sequential pipeline; notebooks in `code/` cover differential expression, GSEA, and visualization.

To render a notebook non-interactively:

```r
rmarkdown::render("code/predictive_modeling/07_rf_classification_eif3d_targets_no_boruta.Rmd")
```

RNA secondary structure features were pre-computed on a compute cluster using scripts in `cluster_scripts/` (ViennaRNA). The resulting flat files are read directly by the feature extraction notebooks.

## Worked Example: Predicting eIF3d-Promoted Genes at 1 hr Hypoxia

This is the project's primary Random Forest model — it asks whether mRNA sequence, structure, and stability features can distinguish genes whose translation efficiency depends on eIF3d under acute hypoxia from genes eIF3d never regulates. Use it as the template for every other notebook 07 run.

### Before you start: `output/` is not in the repository

`output/` is listed in `.gitignore`, so a fresh clone contains **no** feature matrices, gene sets, or model results. Notebook 07 reads three files that you must obtain before it will run:

| File | What it is | Rows |
|---|---|---|
| `output/predictive_modeling/feature_matrix_external_stability.rds` | Transcript-level feature matrix (~4.4 MB) | 10,172 × 134 |
| `output/genesets/hypoxia_3d_promotes_TE_1hr_lfc0.5.csv` | Positive class: eIF3d-promoted genes | 1,604 |
| `output/predictive_modeling/negative_control_genes.csv` | Negative class: never regulated by si3d | 3,319 |

There are two ways to get them:

1. **Ask Kate for the three files** (recommended — they total under 5 MB). Drop them into the paths above, preserving the directory structure.
2. **Rebuild them.** The feature matrix requires the full extraction chain, in this order: `01_` → `01b_` → `01d_` → `01e_` → `01f_` → `01g_` → `01h_` → `22_cnot3_riboseq_slamseq` → `23_dhx29_riboseq_slamseq` → `01i_` → `01j_` → `01k_`. Note that `01j_` and `01k_` **overwrite** `feature_matrix_external_stability.rds` in place, so both must finish before notebook 07 runs. The gene sets come from `code/generate_signatures.Rmd` (positive set) and `code/predictive_modeling/02_negative_control_set.Rmd` (negative controls). This path also needs the GENCODE v49 TxDb from the Setup section and the cluster-generated ViennaRNA files.

### Running it

Open `hypoxia_timecourse.Rproj` first — every notebook resolves paths with `here()`, which needs the project root.

The defaults in `07_rf_classification_eif3d_targets_no_boruta.Rmd` are already set to this exact analysis, so knitting it in RStudio with no changes reproduces the model. To be explicit (and for scripted runs), pass the parameters:

```r
rmarkdown::render(
  "code/predictive_modeling/07_rf_classification_eif3d_targets_no_boruta.Rmd",
  params = list(
    condition               = "hypoxia",
    gene                    = "3d",
    direction               = "promotes",
    timepoint               = "1hr",
    lfc                     = "0.5",
    feature_set             = "external_stability",
    clip_mode               = "aggregate",
    exclude_own_te_features = TRUE
  )
)
```

Two parameters are worth understanding rather than copying blindly:

- **`feature_set = "external_stability"`** selects `feature_matrix_external_stability.rds`. The mapping from `feature_set` to RDS filename is hard-coded in the `derived_values` chunk; you never name the RDS directly (the `feature_matrix_path` parameter is an override reserved for cell-line-specific matrices such as MCF7-SIX1).
- **`exclude_own_te_features = TRUE`** drops the `indiv_te_*`, `delta_te_*`, and `sig_*` columns. This is **required** here, not optional: the positive and negative gene sets are themselves defined by TE, so leaving those columns in leaks the label into the features and produces a meaninglessly high AUC.

Leave `neg_geneset_csv` empty. Because `gene = "3d"`, the notebook automatically picks up `negative_control_genes.csv`; that parameter is only for non-default backgrounds such as the 3d-vs-3e contrast.

### What you get

Results are written to `output/predictive_modeling/` with the suffix
`eif33d_promotes_hypoxia_1hr_lfc0.5_external_stability_clipaggregate_no_boruta_noOwnTE`
(the doubled `eif33d` is expected — the suffix is built as `"eif3"` + the `gene` parameter):

- `rf_classifier_<suffix>.rds` — the trained model
- `rf_classification_importance_<suffix>.csv` — feature importance, the main result
- `rf_classification_performance_<suffix>.csv` and `rf_classification_cv_results_<suffix>.csv` — AUC and cross-validation
- `rf_classification_predictions_<suffix>.csv` and `rf_model_gene_sets_<suffix>.csv` — per-gene predictions and the exact gene sets used

Change any parameter and the suffix changes with it, so runs never overwrite each other. To compare importance across two runs (for example hypoxia vs. normoxia), feed both into `09_cross_condition_importance_comparison.Rmd`.

### Sanity checks

The notebook prints its resolved inputs early — confirm they read:

```
Gene set CSV: hypoxia_3d_promotes_TE_1hr_lfc0.5.csv
Feature matrix: feature_matrix_external_stability.rds
```

`set.seed(9)` is the project-wide seed and is set before every stochastic step, so a correct run is reproducible: your feature importance ranking should match Kate's run exactly. If it doesn't, the feature matrix is the likely culprit — check that both `01j_` and `01k_` were run.

### If you use Claude Code

`CLAUDE.md` in the repository root is the machine-readable version of everything above, plus the project's conventions: contrast directions (which are *not* symmetric between the two cell lines), the feature-extraction chain, validation requirements, and plotting rules. Claude Code loads it automatically on session start, so you can point an agent at this repository and ask it to run or extend notebook 07 without restating the pipeline. Read it yourself too before modifying any notebook — the DESeq2 contrast-direction section in particular documents a sign convention that is easy to get backwards.

## Repository Structure

```
code/
  functions.R                  # Shared DESeq2 wrappers and utility functions
  predictive_modeling/         # Feature extraction + ML pipeline (notebooks 01–58)
  *.Rmd                        # TE analysis, GSEA, visualization notebooks
output/
  predictive_modeling/         # Feature matrices (.rds) and RF model outputs
  genesets/                    # Gene signature CSVs
accessories/
  human/                       # GENCODE v49 annotation and gene tables
  translation_signatures_literature.csv
  csc_data/                    # External CLIP-seq datasets
plots/                         # Generated figures (PDF)
counts/                        # Raw RNA-seq and ribosome profiling count matrices
cluster_scripts/               # ViennaRNA secondary structure scripts
```

## Contact

Kate Matlin — katherine.matlin@ucdenver.edu
