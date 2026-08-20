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
- `02b_negative_control_set_si3e.Rmd` → the si3e equivalent (parameterized `neg_lfc` / `neg_padj`)

**Screened control vs expressed background — two different comparators, two different claims.**
The `02`/`02b` sets are screened: genes selected for having *no* effect. They are not a
neutral background — the si3e screened pool is markedly GC-rich, so part of any GC
separation against it is a property of the screen. `02c_expressed_background_pool_si3e.Rmd`
builds the alternative: every gene DESeq2 could **test** in the contrast (non-NA `te_padj`)
minus every gene in the positive set's parent list, with **no screening on effect size**
(10,778 genes; 9,183 in `feature_matrix_dhx29_kd_feature.rds`). eIF3e-inhibited genes are
deliberately retained — excluding them would rebuild the screened set.

| Comparator | eIF3e-specific model, hypoxia 1hr | Claim |
|---|---|---|
| screened `si3e_negative_controls_padj0.05.csv` | test AUC 0.907 | targets vs eIF3e-**protected** mRNAs |
| expressed background (`02c`) | test AUC 0.769 (seed 9); 0.740 median over 20 draws | targets vs a **typical expressed** mRNA |

The two are not comparable and must not be plotted on the same axis. The pool CSV is a
**definition**, not pre-intersected with any matrix — the models use
`feature_matrix_dhx29_kd_feature.rds` (10,172 genes) and the nb89 GC figure uses
`external_stability` joined to `positional_gc`, so each consumer intersects and reports its
own n.

**Negative-subsample replication (notebooks 90/91).** Notebook 07 keeps every positive and
draws a size-matched negative class from the pool — 399 of 9,183 for the eIF3e-specific
model — so a single run reports the AUC and feature ranking of *one arbitrary draw*.
`90_replicate_neg_subsample_driver.R` renders 07 once per `neg_seed` (20 seeds, ~4 s each);
`91_neg_subsample_stability.Rmd` reports the spread. On the eIF3e-specific model: test AUC
median 0.740, sd 0.029, range 0.683–0.803, with **seed 9 sitting at the 90th percentile** —
the default draw is on the lucky side. Top-20 feature agreement between replicate pairs is
median Jaccard 0.67, and only **7 features are top-20 in every replicate**: `csc_q4`, `csc`,
`log2_cds_length`, `csc_q3`, `log2_utr5_length`, `csc_q2`, `karner2026_mda231_log2_ct`.
Rank a feature off a single 07 run only if it clears that bar.

**Null-band GC figures (`89_gc_position_density_by_geneset.Rmd`, section 4b).** Against the
02c background the comparator is drawn from the same population as the grey curve, so the
black curves become a **null band** — what N random expressed mRNAs look like — rather than
a contrast. `gc_null_band_figure()` plots 50 draws and computes an empirical p on 1,000
(`(k+1)/(n+1)`, so the floor is 1/1001, not 0). Draws are cheap here and are deliberately
**decoupled from the 20 model replicates**, which are limited by forest-fitting cost. On the
eIF3e-specific set only GC3 clears the band decisively (median 62.8% vs null 57.4%, z = 3.6,
p < 0.001); CDS GC is marginal (53.1 vs 51.5, p = 0.019) and GC1/GC2 are inside it
(p = 0.069 / 0.159). Against the *screened* pool all four looked separated — that difference
is the screen, not the biology.

**RF classification models (03–10):**
- `07_rf_classification_eif3d_targets_no_boruta.Rmd` — primary model, fully parameterized via YAML; key params:
  - `condition`, `gene`, `direction`, `timepoint`, `lfc`, `feature_set`, `clip_mode`
  - `neg_geneset_csv` — custom negative control CSV from `output/genesets/` (default uses `negative_control_genes.csv`)
  - `feature_matrix_path` — override to load a named RDS from `output/predictive_modeling/` (e.g. for MCF7-SIX1)
  - `exclude_own_te_features` — removes `indiv_te_*`, `delta_te_*`, `sig_*` columns; set TRUE when pos/neg sets are defined by TE
  - `neg_seed` — seeds **only** the negative draw in `create_target`. The train/test split,
    the forest and the CV stay pinned at 9, so across a seed sweep every difference is
    attributable to negative-class composition. Default 9 = the historical draw.
  - `neg_tag` — appended to the output suffix. **Required whenever you change the negative
    POOL** (not just the draw): `neg_geneset_csv` is *not* part of the suffix, so a new pool
    at `neg_seed = 9` silently overwrites the model, importance and predictions of the run
    that used the old pool. `"exprbg"` is the tag for the 02c background pool. Both
    `neg_seed` and `neg_tag` are mirrored in 07b/07c, which rebuild the draw themselves and
    would otherwise explain a model trained on different genes.
  - `save_plots` — FALSE suppresses the 8 `ggsave` calls; use for seed sweeps.
  - Notebook 07 now **removes positives from the negative pool before sampling**. Previously
    an overlapping gene was sampled, then relabelled positive by the `case_when`, silently
    shrinking the negative class below `number_to_match`. It also reports when the pool is
    no larger than the positive set, where every seed draws the same genes and a sweep
    produces N copies of one model.
  - `readout` — `"TE"` (default) or `"RNA"`. Selects the gene set filename (`..._TE_...` vs `..._RNA_...`) and is appended to the output suffix **only when not "TE"**, so every pre-existing TE run keeps its filenames and 07c still resolves its saved models. Mirrored in 07b/07c/08/09. `"RNA"` requires an explicit `neg_geneset_csv` (the default `negative_control_genes.csv` is TE-defined and 307 of its genes satisfy the RNA positive definition) and is rejected with `timepoint="alltimepoints"/"1and4hr"` only for TE, since that shortcut filters `te_lfc` inline.
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

**RNA-level (abundance) readout — notebooks 70–73.** Counterpart to the TE pipeline, asking which mRNAs *lose abundance* on si3d rather than which lose TE.

- `70_rna_genesets_si3d.Rmd` — builds the RNA positive + negative gene sets. YAML-parameterized (`pos_condition`, `pos_timepoint`, `pos_lfc`, `neg_padj`); defaults reproduce the hypoxia 1hr configuration byte-identically. Excludes **EIF3D** from the positive set (the siRNA depleting its own target is not eIF3d-dependent biology) while asserting it *would* have qualified, which doubles as the direction check.
- `71_gc_content_vs_rna_changes_transcriptome.Rmd` — GC/length vs RNA and TE change, whole transcriptome.
- `72_rna_geneset_overlap_conditions.Rmd` — hypoxia vs normoxia RNA gene set overlap (Venn), plus a power/threshold decomposition.
- `73_global_shift_rna_vs_te_by_timepoint.Rmd` — RNA vs TE LFC densities per timepoint.

**The 0.5 LFC convention does NOT transfer to raw RNA.** These contrasts are spike-in normalized and si3d lowers total mRNA per cell, so `rna_lfc` is not centred on zero and the median drifts with time:

| median `rna_lfc` | 1hr | 4hr | 24hr | 1+4hr regressed |
|---|---|---|---|---|
| normoxia | 0.368 | 0.587 | 0.838 | 0.505 |
| hypoxia | 0.409 | 0.176 | 0.209 | — |

TE medians stay near zero (normoxia 0.109 → 0.164) because the ratio cancels the per-cell loss. Consequences:

- `rna_lfc > 0.5` at hypoxia 1hr selects 4,764 genes (~40% of the transcriptome) — essentially "above median". RNA gene sets use **`lfc > 1`**, and even that is only 0.50–0.63 above the median depending on contrast.
- **Set sizes are not comparable across contrasts at a fixed raw threshold.** Normoxia 1+4hr yields 1,282 genes vs normoxia 1hr's 705 purely because its median is 0.137 higher; matched on a centred cut they are 888 vs 841. The `padj < 0.05` filter is entirely non-binding at `lfc > 1` — set size is set by the LFC cut alone.
- **Notebook 07 keeps ALL positives and downsamples only the negative pool** to `min(n_pos, n_neg)`. Every TE run has negatives > positives so this balances silently; RNA runs can invert it. Check the negative pool clears the positive count (nb70 asserts this), and relax `neg_padj` (0.2 → 0.1) rather than the LFC bound, which is never binding.
- The RNA negative set is **eIF3d-protected, not unaffected** — its median `rna_lfc` (~0.03) sits well *below* the transcriptome median, so those transcripts retained mRNA while the typical one lost 25–30%. Describe results as "eIF3d-dependent vs eIF3d-protected".

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

**MDA-MB-231 DESeq2 contrast direction — the filenames are backwards. Read this before labelling any plot.**

The DESeq2 call inside `eif3e_eif3d_normoxia_and_hypoxia.Rmd` runs **sictrl as `test_condition` and si3d as `control_condition` (reference)**. In `code/functions.R`:

```r
243:  Condition = factor(comparison_variable, levels = c(control_condition, test_condition))  # control = reference
320:  coef_name <- paste0("Condition_", test_condition, "_vs_", control_condition)             # -> "Condition_sictrl_vs_si3d"
```

So the coefficient DESeq2 actually emits is named **`sictrl_vs_si3d`**, and every column in these files is **log2(sictrl / si3d)**. The files are nonetheless written as `translation_categories_si3d_vs_sictrl_*.csv` — **the filename states the ratio in the opposite order from the coefficient it contains.** This is a known wart, kept because ~350 reference sites across ~54 notebooks depend on the current names; it is not a claim about the direction.

- Positive `te_lfc` / `rna_lfc` = sictrl > si3d = knockdown **reduces** it = eIF3d **promotes**.
- **Plot titles and axis labels must say "sictrl vs si3d" (or "log2(sictrl / si3d)"), never "si3d vs sictrl".** Copying the filename into a caption produces a figure whose title contradicts its own axis. This was caught in notebooks 71/72/73 after the fact.
- Verify empirically rather than trusting either name: **EIF3D's own `rna_lfc` must be strongly positive** (~+2.9 at hypoxia 1hr, ~87% knockdown), because the siRNA depletes its target. If it is negative, the direction has flipped. The same check on the si3e files puts EIF3E at ~+4.9. New notebooks reading these files should assert this (see nb70/71/72/73 for the pattern).

Do not confuse with MCF7-SIX1 (below).

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
