# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-only bioinformatics project studying translational regulation by eIF3d and eIF3e during hypoxia in MDA-MB-231 and MCF7-SIX1 breast cancer cells. It integrates ribosome profiling + RNA-seq to compute translation efficiency (TE), then builds Random Forest models to predict which genes are translationally regulated by eIF3d/3e and what sequence/structural features drive that regulation.

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
| `01c_kmer_motif_precompute` | (intermediate cache) | Pre-computes k-mers/RBP motifs; must run before `01c_` |
| `01_` | `feature_matrix_1hr.rds` | UTR/CDS lengths, GC%, k-mers, RBP motifs (~900 features) |
| `01b_` | `feature_matrix_1hr_targeted.rds` | eIF3 CLIP binding (7 conditions) |
| `01c_` | `feature_matrix_1hr_combined.rds` | Combined all above (~6,300 features) |
| `01d_` | `feature_matrix_1hr_structure.rds` | RNA secondary structure (RNAplfold/RNAfold) |
| `01e_` | `feature_matrix_1hr_clip_split.rds` | Split CLIP into `clip_eif3_count` / `clip_cap_binding_count` |
| `01f_` | `feature_matrix_te_lfc.rds` | TE features per condition (8 columns) |
| `01g_` | `feature_matrix_codon_optimality.rds` | CSC, GC3 codon optimality |
| `01h_` | `feature_matrix_kozak.rds` | Kozak sequence features |
| `22_cnot3_riboseq_slamseq` | `feature_matrix_cnot3.rds` | Zhu 2024 SLAM-seq half-lives + CNOT3 KO riboseq |
| `23_dhx29_riboseq_slamseq` | `feature_matrix_dhx29.rds` | Hia 2026 DHX29 SLAM-seq + riboseq occupancy |
| `01i_` | `feature_matrix_external_stability.rds` | Karner 2026 MDA-MB-231 SLAM-seq (reads `feature_matrix_dhx29.rds`) |
| `01j_` | `feature_matrix_csc_positional.rds` | Positional CSC (full CDS + quarters + first 75 codons) |
| `01k_` | (appends to csc_positional) | CDS/3'UTR structure features |
| `01l_` | `feature_matrix_codon_features.rds` | bg-corrected RSCU + raw codon frequencies, 121 features (reads `feature_matrix_external_stability.rds`) |
| `01m_` | `feature_matrix_codon_only.rds` | Codon features without external stability base |
| `43_` | `feature_matrix_mcf7six1_codon_features.rds` | MCF7-SIX1 codon features with cell-line-specific transcript selection |

The chain `22_cnot3` → `23_dhx29` → `01i_` is non-obvious: two non-`01_`-prefixed notebooks are required steps before `01i_` can run. RNA secondary structure inputs (from `01d_`) were pre-computed on a compute cluster using `cluster_scripts/` and read in as flat files — they cannot be recomputed locally.

**Negative control definition:**
- `02_negative_control_set.Rmd` → `output/predictive_modeling/negative_control_genes.csv`
- Negative controls = genes with |TE_LFC| < 0.5 AND padj > 0.3 across ALL si3d conditions

**RF classification models (03–10):**
- `07_rf_classification_eif3d_targets_no_boruta.Rmd` — primary model, fully parameterized via YAML; key params:
  - `condition`, `gene`, `direction`, `timepoint`, `lfc`, `feature_set`, `clip_mode`
  - `neg_geneset_csv` — custom negative control CSV from `output/genesets/` (default uses `negative_control_genes.csv`)
  - `feature_matrix_path` — override to load a named RDS from `output/predictive_modeling/` (e.g. for MCF7-SIX1)
  - `exclude_own_te_features` — removes `indiv_te_*`, `delta_te_*`, `sig_*` columns; set TRUE when pos/neg sets are defined by TE
- `09_cross_condition_importance_comparison.Rmd` — compare feature importance across two 07 runs; key params:
  - `condition_hyp`, `condition_nor` — condition strings used to construct the 07 output file suffixes (default `"hypoxia"` / `"normoxia"`; set to e.g. `"mcf7six1_hypoxia"` for MCF7-SIX1)

**Analysis notebooks (11–43):** Feature-specific deep dives, external validations, regression models, and MCF7-SIX1 analyses.

## Key Data Objects and Conventions

**Gene identifiers** — three related columns always present together in feature matrices and output tables:
- `gene_id` — versioned ENSEMBL ID (e.g. `ENSG00000000003.14`)
- `gene_id_clean` — version-stripped: `sub("\\..*", "", gene_id)`
- `symbol` — HGNC gene symbol
- `transcript_id_clean` — version-stripped transcript ID (primary key in feature matrices)

**Feature matrices** are RDS files with one row per gene/transcript, keyed by `transcript_id_clean`. List columns may be present (e.g. nested codon data from `coRdon`); coerce before using dplyr aggregation: `as.integer(unlist(col))`.

**Translation efficiency (TE)** = log2(ribosome footprint / RNA) computed under si3d vs sictrl knockdown. Positive TE LFC = eIF3d promotes translation; negative = inhibits.

**MCF7-SIX1 directionality differs:** DESeq2 contrast is log2(si3d / sictrl), so positive LFC = eIF3d **inhibits**. Always negate: `te_lfc = -log2FoldChange` so positive = promotes (matching MDA-MB-231 convention). MCF7-SIX1 TE data lives in a separate repo: `/Users/katematlin/github/2024_eIF3e_hypoxia/2024_eIF3e_hypoxia/`. Salmon RNA-seq samples for transcript selection: `NM2023_0101`, `NM2023_0111`, `NM2023_0121` (rRNADepletion/KAPA_RNAHyper). Ribosome profiling samples are `NM2023_0037–0066` — do not use for transcript selection.

**Gene group conventions** (used across all visualization notebooks):
```r
group_colors <- c(
  "eIF3d hypoxia promotes" = "#E69F00",
  "Negative control"       = "grey40",
  "All other genes"        = "grey80"
)
```
Positive gene set: `output/genesets/hypoxia_3d_promotes_TE_1hr_lfc0.5.csv` → column `ensembl_gene`; also `db_gene_symbol` for notebook 09 Venn diagrams.
Negative controls: `output/predictive_modeling/negative_control_genes.csv`, size-matched with `set.seed(9)`.

**Gene set CSV naming convention:** `{condition}_{gene}_{direction}_TE_{timepoint}_lfc{lfc}.csv` in `output/genesets/`. Notebooks 07 and 09 construct paths from params using this pattern. For MCF7-SIX1: `mcf7six1_hypoxia_3d_promotes_TE_1hr_lfc0.5.csv`; negative controls go in `output/genesets/` (not `output/predictive_modeling/`) when using the `neg_geneset_csv` param.

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
| `output/predictive_modeling/negative_control_genes.csv` | Negative control gene pool (MDA-MB-231 si3d) |
| `accessories/all_combined_sigs_1_13_26.csv` | All gene signatures across conditions and cell lines; column `gs_name` identifies set (e.g. `MCF7-SIX1_hypoxia_3d_promotes_TE`) |
