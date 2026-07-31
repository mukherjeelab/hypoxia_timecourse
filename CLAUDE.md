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

## Shared Helper Functions (`code/functions.R`)

Check here before writing new helpers. Three key functions:

- **`deFunction(mode, ...)`** — main DESeq2 wrapper; modes: `"TE"` (interaction design for translation efficiency), `"RNA"` (RNA-seq only), `"RIBO"` (ribosome profiling only). Two-stage spike-in handling: spike-ins are used to estimate size factors, then **removed from the count matrix before fitting** to avoid inflating gene counts. Contains fallback logic for DESeq2 coefficient naming convention changes.
- **`deFunction_time_normalized()`** — extends `deFunction()` by adding `hours` as a batch covariate; design becomes `~ hours + Condition + SeqType + Condition:SeqType`.
- **`categorize_translation_changes()`** — classifies genes into Forwarded / Exclusive / Buffered / TE_only categories. The function's own default parameter is `lfc_cutoff = log2(1.5) ≈ 0.585`, but this default is never actually used — every call site in the codebase (`rna_ribo_te_analysis.Rmd`, `timecourse_qc.Rmd`, `eif3e_eif3d_normoxia_and_hypoxia.Rmd`, 39 calls total) explicitly overrides it to `lfc_cutoff = 0.0`, so the `category` column reflects significance by `padj` alone with no magnitude filter. **The project's real fold-change standard is 0.5** (used for gene sets, negative controls, and everywhere else) — it's applied as a separate filter directly on `te_lfc`/`te_padj` *after* this function runs, not through its `lfc_cutoff` argument. Don't read the `log2(1.5)` default as "the project convention"; it isn't exercised anywhere.
- **`createNormMatrix()`** — DESeq2 normalization with optional spike-in support; uses the TE interaction design.

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
| `01i_` | `feature_matrix_external_stability.rds` | Karner 2026 MDA-MB-231 SLAM-seq (reads `feature_matrix_cnot3.rds` as its base, then re-joins the DHX29 columns directly from `dhx29_slamseq_halflives.rds` + `dhx29_riboseq_occupancy.rds` — so `23_` must still run first, but `feature_matrix_dhx29.rds` itself is not consumed by this pipeline) |
| `01j_` | overwrites `feature_matrix_external_stability.rds` | Positional CSC (full CDS + quarters + first 75 codons) |
| `01k_` | overwrites `feature_matrix_external_stability.rds` | CDS/3'UTR structure features (drops 5'UTR RNAfold cols) |
| `01l_` | `feature_matrix_codon_features.rds` | bg-corrected RSCU + raw codon frequencies, 121 features (reads `feature_matrix_external_stability.rds`) |
| `01m_` | `feature_matrix_codon_only.rds` | Codon features without external stability base |
| `01n_` | `feature_matrix_ripseq.rds` | eIF3d/eIF4E RIP-seq IP/IgG log2FC + IP/RNA log2 ratio, 4 timepoints, 16 features (reads `feature_matrix_external_stability.rds`; use with 07 via `feature_matrix_path`) |
| `01v_` | `feature_matrix_positional_gc.rds` | GC1/GC2/GC3 per codon position + `gc3_residual` (GC3 orthogonalized against `cds_gc`). Reads `feature_matrix_external_stability.rds` and writes a **new** file, leaving the base untouched. By default drops `cds_gc` and `transcript_gc` (param `drop_aggregate_gc`) — `cds_gc` is ~the mean of the three positions, so leaving it in lets the forest fall back on the aggregate and the positional split never gets tested |
| `43_` | `feature_matrix_mcf7six1_codon_features.rds` | MCF7-SIX1 codon features with cell-line-specific transcript selection |

The chain `22_cnot3` → `23_dhx29` → `01i_` → `01j_` → `01k_` is non-obvious. `01j_` and `01k_` overwrite `feature_matrix_external_stability.rds` in place — both must run after `01i_` before notebook 07 is called. RNA secondary structure inputs (from `01d_`, `01k_`) were pre-computed on a compute cluster using `cluster_scripts/` and read in as flat files. `01c_` (k-mer + RBP motif features) can be skipped for the `external_stability` pipeline — `01d_` is patched to read directly from `01_` + `01b_` outputs instead.

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

**`07b_` and `07c_` are downstream of `07_` — run 07 first.** Both rebuild 07's gene sets, features and split from the same params (deterministic under `set.seed(9)`), so their params must match the 07 run being analyzed. 07c additionally **loads the saved `rf_classifier_{suffix}.rds` and asserts against it**, so it fails with `notebook 07 model RDS not found` if 07 has not been knit with those params. Full order for a new feature set: **`01x_` (build matrix) → `07_` → `07c_`**.

- `07b_rf_permutation_importance_eif3d_targets.Rmd` — refits with `importance = "permutation"` as a less-biased check on the Gini ranking.
- `07c_rf_shap_eif3d_targets.Rmd` — per-gene Shapley decomposition; gives **signed contributions per gene per feature**, which 07/07b cannot. Key points:
  - `treeshap` does not support ranger **probability** forests, so 07c fits a **surrogate regression forest** on the 0/1 target with `mtry`/`min.node.size`/`num.trees` read off the saved 07 model (never hardcoded — ranger's defaults differ by tree type; `min.node.size` is 10 for probability but 1 for classification and 5 for regression). Check A asserts the config matches structurally; Check B asserts the surrogate reproduces 07's held-out predictions within `fidelity_tol` (0.02, calibrated on the 56-feature `external_stability` set against deliberate misconfigurations — treat as a sanity bound, not a tuned one, on other feature sets).
  - SHAP is additive, so summing within a feature family is exact. This is the only way to resolve credit-splitting across correlated blocks (`csc`/`csc_first75`/`csc_q1`–`q4`, the eight `struct_accessibility_*` families, GC and length families) — impurity importance structurally cannot.
  - **`shap_value_cor` is the directional quantity** (correlation of a feature's value with its own SHAP value). `mean_signed_shap` is **not** a measure of direction — it averages to ~0 for every feature and its residual offset reflects the gap between explain-set mean prediction and the training baseline. It is exported for completeness only.
  - `run_interactions` param — TreeSHAP interaction values scale with feature count (~10 min at 56 features on 200 genes, far longer at 177). Set FALSE for large feature sets. Results cache to `shap_interactions_{suffix}.rds` and reload on re-knit.
  - Sections that need the `cds_gc`/`csc`/`hia2026_dhx29_occupancy` trio, or per-codon features, are chunk-gated (`eval=run_pair_diag` / `run_codon_tests` / `run_gc3_screen`) and skip cleanly on matrices that lack them.
  - **Use a distinct `feature_set` label for a new matrix.** The output suffix is built from params, so reusing an existing label silently overwrites that run's outputs. Existing labels used with `feature_matrix_path`: `codon_rscu` (→ `feature_matrix_codon_features.rds`), `positional_gc` (→ `feature_matrix_positional_gc.rds`).

**Codon analysis: use bg-RSCU, not raw codon frequency, to test synonymous preference.** `codon_freq_AAA` and `codon_freq_AAG` both rise with lysine content, so raw-frequency features are dominated by amino-acid composition and report "same direction" whether or not a synonymous preference exists. `bg_rscu_*` is normalized within each synonymous family, so its SHAP direction isolates wobble choice. In the si3d hypoxia 1hr contrast these dissociate cleanly: raw frequency tracks codon positions 1–2 (ρ = −0.39, amino-acid identity) while bg-RSCU tracks position 3 (ρ = −0.67, synonymous choice). Aggregate `gc3` cannot substitute — it is 0.93-correlated with `cds_gc`, 0.66 with `gc1` via isochore structure, and nets out opposing within-family preferences.

## Integrating External Data from Papers

Before incorporating any dataset from a published paper into the feature matrix or analysis code, first verify that the data reproduces a key result from that paper. This protects against using the wrong file, a misread column, or a pre/post-processing mismatch that would silently corrupt downstream models.

**What to replicate:** A single quantitative result that is easy to check — a reported median half-life, a named gene's enrichment score, a summary statistic from a figure, or the gene count in a reported set. It does not need to be a full analysis.

**Where to put the check:** In the notebook section that loads the external data, immediately after reading the file and before any joins. Use `stopifnot` or an explicit comparison with `cat` reporting the expected vs. observed value.

Example pattern:
```r
# Zhu 2024: paper reports median sgNT half-life ~3.5h (Fig 2B)
stopifnot("Zhu2024 median half-life outside expected range" =
  between(median(halflives$sgNT_halflife, na.rm = TRUE), 2.5, 5.0))
cat("Zhu2024 median sgNT half-life:", round(median(halflives$sgNT_halflife, na.rm = TRUE), 2), "h\n")
```

## Feature Extraction Validation Conventions

Every feature extraction notebook (including new feature groups added to `43_`) **must** include both types of checks for each feature group added:

### 1. Computational checks (`stopifnot`)
Assert value ranges that must hold by definition. Use `stopifnot()` — never `cat(check)`, which passes silently on failure.

| Feature type | Required check |
|---|---|
| GC content (%) | `all(x >= 0 & x <= 100, na.rm = TRUE)` |
| Accessibility (RNAplfold) | `all(x >= 0 & x <= 1, na.rm = TRUE)` |
| MFE (RNAfold) | `all(x <= 0, na.rm = TRUE)` |
| CAI, tAI, FOP | `all(x >= 0 & x <= 1, na.rm = TRUE)` |
| -log10 enrichment (riboseq) | `all(x >= 0, na.rm = TRUE)` |
| Half-lives | `all(x > 0, na.rm = TRUE)` |
| Log-fold changes | `all(is.finite(x[!is.na(x)]))` |
| CLIP counts | range within [0, n_conditions] |
| Feature column presence | `stopifnot(all(expected_cols %in% colnames(df)))` |
| Leakage guard | `stopifnot(!"target_col" %in% colnames(df))` |

### 2. Biological spot-checks
Assert known invariants about specific genes or codon families. These catch implementation errors that pass range checks.

Established examples (copy/adapt for new notebooks):
- **ACTB CDS GC%**: expected 50–70% (`01_`)
- **Ribosomal proteins have above-median CAI**: `median(ribo_cai) > median(all_cai)` (`01g_`)
- **CTG is optimal Leu codon**: `w["CTG"] == 1` (`01g_`)
- **Kozak PWM monotonicity**: mean PWM increases tier 0 → 3 (`01h_`, `43_`)
- **GC3 median in [40, 65]%** (`01_`, `43_`)

Both check types must be present. `cat()` summaries are informative but are not substitutes for `stopifnot()`.

## Reproducibility: Random Seed Convention

**Always use `set.seed(9)`** immediately before any stochastic operation — GSEA/fgsea permutation testing, negative control gene sampling, Random Forest training, train/test splits, or any other call that depends on R's RNG. This is the project-wide standard seed; do not use `42` or other values in new code. (Existing `set.seed(42)` calls in `code/acute_to_long_hypoxic_response.Rmd` predate this convention and have not been retroactively changed, since doing so would shift already-reported results — flag before touching them.)

## Key Data Objects and Conventions

**Gene identifiers** — three related columns always present together in feature matrices and output tables:
- `gene_id` — versioned ENSEMBL ID (e.g. `ENSG00000000003.14`)
- `gene_id_clean` — version-stripped: `sub("\\..*", "", gene_id)`
- `symbol` — HGNC gene symbol
- `transcript_id_clean` — version-stripped transcript ID (primary key in feature matrices)

**Feature matrices** are RDS files with one row per gene/transcript, keyed by `transcript_id_clean`. List columns may be present (e.g. nested codon data from `coRdon`); coerce before using dplyr aggregation: `as.integer(unlist(col))`. Feature matrices are transcript-keyed, not gene-keyed — if you need gene-level operations, deduplicate first: `distinct(gene_id_clean, .keep_all = TRUE)`.

**Translation efficiency (TE)** = log2(ribosome footprint / RNA) computed under si3d vs sictrl knockdown. Positive TE LFC = eIF3d promotes translation; negative = inhibits.

**MDA-MB-231 DESeq2 contrast direction** — despite files being named `translation_categories_si3d_vs_sictrl_*.csv`, the DESeq2 call inside `eif3e_eif3d_normoxia_and_hypoxia.Rmd` runs **sictrl as `test_condition` and si3d as `control_condition` (reference)**. The extracted coefficient is therefore log2(sictrl / si3d). Positive te_lfc = sictrl TE > si3d TE = knockdown reduces TE = eIF3d promotes. The file naming describes the experiment ("si3d knockdown study vs sictrl"), not the DESeq2 numerator/denominator. Do not confuse with MCF7-SIX1 (below).

**MCF7-SIX1 directionality differs:** DESeq2 contrast is log2(si3d / sictrl), so positive LFC = eIF3d **inhibits**. Always negate: `te_lfc = -log2FoldChange` so positive = promotes (matching MDA-MB-231 convention). MCF7-SIX1 TE data lives in a separate repo: `/Users/katematlin/github/2024_eIF3e_hypoxia/2024_eIF3e_hypoxia/`. Salmon RNA-seq samples for transcript selection: `NM2023_0104`, `NM2023_0114`, `NM2023_0124` (sictrl normoxia, KAPA_RNAHyper). Matches MDA-MB-231 convention. Ribosome profiling samples are `NM2023_0037–0066` — do not use for transcript selection.

**Gene group conventions** (used across all visualization notebooks):
```r
group_colors <- c(
  "eIF3d hypoxia promotes" = "#E69F00",
  "Negative control"       = "grey40",
  "All other genes"        = "grey80"
)
```

**eIF3d / eIF3e color convention (required in all new plots):**
```r
kd_colors <- c(
  "si3d" = "#E69F00",  # Okabe-Ito orange       -> always eIF3d / si3d
  "si3e" = "#56B4E9"   # Okabe-Ito light blue   -> always eIF3e / si3e
)
```
- eIF3d / si3d and eIF3d-defined gene sets: Okabe-Ito orange `#E69F00`
- eIF3e / si3e and eIF3e-defined gene sets: Okabe-Ito light (sky) blue `#56B4E9`
- Negative control set: `black`
- Rest of the transcriptome / all other genes: grey (`grey70`)

These two hues are reserved for the two factors. When a plot also needs up/down or
significance colors, pick from the remaining Okabe-Ito hues (e.g. `#0072B2` dark blue,
`#CC79A7` reddish purple) so nothing else in the figure reads as eIF3d or eIF3e.
Applied in `code/predictive_modeling/66_rna_changes_si3d_si3e.Rmd`. The older
`group_colors` block above keeps `grey40`/`grey80` for backward compatibility with
already-generated figures; new plots should use `black` / `grey70`.
Positive gene set: `output/genesets/hypoxia_3d_promotes_TE_1hr_lfc0.5.csv` → column `ensembl_gene`; also `db_gene_symbol` for notebook 09 Venn diagrams.
Negative controls: `output/predictive_modeling/negative_control_genes.csv`, size-matched with `set.seed(9)`.

**Gene set CSV naming convention:** `{condition}_{gene}_{direction}_TE_{timepoint}_lfc{lfc}.csv` in `output/genesets/`. Notebooks 07 and 09 construct paths from params using this pattern. For MCF7-SIX1: `mcf7six1_hypoxia_3d_promotes_TE_1hr_lfc0.5.csv`; negative controls go in `output/genesets/` (not `output/predictive_modeling/`) when using the `neg_geneset_csv` param.

**RSCU output naming convention** (notebooks 41–47):
- Single-condition bg-RSCU: `bg_rscu_{cell_line}_{condition}_{timepoint}.csv` (e.g. `bg_rscu_mcf7six1_hypoxia_1hr.csv`)
- Delta-delta (cross-condition): `bg_rscu_delta_delta_{cond1}_vs_{cond2}.csv` (e.g. `bg_rscu_delta_delta_mcf7six1_hyp_vs_nor.csv`)
- Raw codon frequencies (intermediate): `codon_freq_{cell_line}_{condition}.csv`
- Sequential dependency: nb47 (delta-delta) requires nb46's output file. If nb46 has not been knit, nb47 will fail at the `stopifnot` file-existence check.

**Namespace conventions:**
- Always use `dplyr::select()`, never bare `select()` — avoids conflicts with `MASS`, `Biostrings`, and other packages that export `select`
- Always use `dplyr::rename()`, never bare `rename()` — same conflict risk
- Always use `dplyr::count()`, never bare `count()` — same conflict risk (e.g. `plyr`/other packages exporting `count()`); bare `count()` has failed with `Error in count(): Argument 'x' is not a vector: list` on an otherwise-valid tibble

**Plot conventions:**
- Save to `plots/` via `here("plots", "notebooknum_description.pdf")`
- `theme_classic(base_size = 12)`, density `alpha = 0.15`, `linewidth = 0.8`
- Wilcoxon tests annotated with `x = Inf, y = Inf, hjust = 1.1, vjust = 1.4`
- Do not use em dashes (`—`) in axis labels — they render as `...` in PDF output. Use plain text only.
- Diverging heatmaps (e.g. `pheatmap` NES/enrichment heatmaps): Okabe-Ito dark blue `#0072B2` (low) to white to Okabe-Ito vermillion/red-orange `#D55E00` (high) — `colorRampPalette(c("#0072B2", "white", "#D55E00"))`. Used in `gsea_si3d_vs_sictrl.Rmd` / `gsea_si3e_vs_sictrl.Rmd`.

## Reference Files

| File | Contents |
|------|---------|
| `accessories/metadata_231_timecourse.csv` | MDA-MB-231 sample metadata: columns `id`, `type` (RiboSeq / Input_Ribo_RNAseq), `rep`, `treatment`, `condition`, `hours`; used in `01_` to select siCTRL normoxia RNA-seq samples for transcript isoform selection |
| `accessories/human/txdb.gencode49.sqlite` | GENCODE v49 transcript database (built by `00_build_txdb_v49.R`) |
| `accessories/human/gene_anno_hs_dm_v49_r111.tsv` | Gene annotation: gene_id, symbol, biotype |
| `accessories/translation_signatures_literature.csv` | Published translation signature gene sets |
| `output/translation_categories_si3d_vs_sictrl_normoxia_1and4hr.csv` | Full TE differential results |
| `output/predictive_modeling/negative_control_genes.csv` | Negative control gene pool (MDA-MB-231 si3d) |
| `counts/all_combined_sigs_1_28_26.csv` | All gene signatures across conditions and cell lines; column `gs_name` identifies set (e.g. `MCF7-SIX1_hypoxia_3d_promotes_TE`) |
| `accessories/human/human_trna_gcn.csv` | tRNA gene copy numbers (GCN) per codon; used as background for RSCU computation in codon optimality notebooks |
