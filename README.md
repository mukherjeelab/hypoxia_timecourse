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
