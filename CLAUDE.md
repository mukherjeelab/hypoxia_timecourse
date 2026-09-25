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

## Translational Control Overview (`code/translational_control_overview/`)

A **fork** of the notebook-07 modelling question, built to run across many datasets.
`code/predictive_modeling/07*` is **frozen** and remains the reproducible baseline; nothing
here writes to `output/predictive_modeling/`. Outputs go to
`output/translational_control_overview/`.

| Notebook | Role |
|---|---|
| `00_feature_reconciliation.Rmd` | audit of our features vs the collaborator repo; writes `feature_inventory_comparison.csv` |
| `01_build_feature_matrix.Rmd` | single builder -> `feature_matrix_tco.rds` + **`feature_annotation.csv`** |
| `02_rf_model.Rmd` | fork of 07, driven by `include_blocks` instead of eight `feature_set` branches |
| `01b_feature_dictionary.R` | regenerates **`FEATURES.md`**, the committed feature dictionary |
| `01c_mutation_test_checks.R` | renders `01_` once per deliberate corruption and asserts the matching check FIRES |
| `03_shap.Rmd` | SHAP decomposition; family attribution is the reason it exists |
| `02_sweep_presets.R` | the sweep arms, shared by `02b_` and `02c_` so a reported result stays re-runnable |
| `02b_neg_subsample_sweep.R` | `Rscript 02b_... <preset>`; 20 `neg_seed` x N arms (~5 min for 3 arms). No argument lists the presets |
| `02c_sweep_stability.Rmd` | reports any preset (`params$preset`); **the null band any claim must clear** |

**`fidelity_tol = 0.02` from `07c_` does not transfer, and `max_abs` is the wrong statistic.**
`03_shap.Rmd` derives its own bound instead of inheriting one, and measuring it shows why: on
this 72-feature model, `max |deviation|` ranks deliberately misconfigured surrogates as
*better* than the correct one (`min.node.size = 25` scores 0.069 vs the correct 0.111), because
it is a single worst-case point. Worse, **`min.node.size` is not detectable by fidelity here at
all** - every value from 1 to 25 lands within a few percent of the correct RMSE, where on the
56-feature `external_stability` set `min.node.size = 1` gave a 12x separation. Only `mtry`
moves fidelity. So Check B asserts on RMSE and **only against the misconfigurations it can
discriminate**, and `min.node.size` is guarded by Check A's structural equality instead.

**SHAP subgroups are real, but the statistic decides the answer - and the obvious one is
wrong.** `03_` section 6 clusters the 287 held-out positive-class genes on family-level SHAP
profiles (exact, since SHAP is additive), normalised within gene so clusters reflect *which*
families drove a gene rather than how strongly it was called. The null permutes each family
independently across genes, destroying co-occurrence while preserving every marginal.

Mean silhouette falls monotonically with k, so the null's best-over-k sits almost always at
k = 2. Comparing that against an observed best at k = 3 compares two different quantities, and
it returns **p = 0.74 and a confident "no structure"**. Standardising **within k** first and
only then maximising over k - which is what pays for having been free to choose k - gives
observed z = 10.9 against a null max of 3.08, **p = 0.005**. Every k >= 3 exceeds the null's
*maximum* over 200 permutations; only k = 2 sits below it, so the structure is real but is not
a two-way split.

| k | observed | null median | null max | z |
|---|---|---|---|---|
| 2 | 0.2301 | 0.2554 | 0.2724 | -4.56 |
| 3 | 0.2524 | 0.2013 | 0.2184 | +7.44 |
| 4 | 0.2128 | 0.1513 | 0.1987 | +3.20 |
| 5 | 0.1974 | 0.1403 | 0.1618 | +6.00 |
| 6 | 0.1964 | 0.1339 | 0.1504 | +10.90 |

`03_subgroup_null.R` holds the one implementation both notebooks use, and **keeps two `best_k`
apart**: `best_k_z` for testing (peaks at k = 6, where the null variance is tightest - a
property of the null, not of the structure) and `best_k_sil` for describing (k = 3, the raw
silhouette peak). At k = 3 the groups are led by `gc` (n = 135), `length` (n = 105) and a weak
`clip`/`codon_optimality` group (n = 47), and the PCA shows a continuous cloud, so treat
membership as a gradient rather than a hard class. `03b_shap_subgroup_null.Rmd` carries the
figures, including the raw statistic as a worked counter-example.

**The error was caught by plotting the per-k null band**, not by the analysis - the single
summary number hid that observed sat above the null at every k >= 3. Plot the band.

**`02_` saves the imputed split** (`rf_model_data_{suffix}.rds`, gated on `save_model_data`,
off during sweeps) so `03_` explains *that* model's data rather than re-deriving it from
params as `07c_` must. That removes a whole class of drift rather than asserting its way back.

**`FEATURES.md` is the committed feature dictionary** - one row per column with variant,
origin, coverage and definition, plus which `validate_*` chunk exercises it. Regenerate with
`01b_feature_dictionary.R` after any change to `01_`; never edit it by hand. It exists because
`feature_annotation.csv` lives in `output/`, which is gitignored, so the matrix's contents are
otherwise invisible to anyone reading the repo. Definitions are written **only where `01_` states one deliberately** - 61 of 171 columns. A
blank is not an oversight: the column is inherited, and paraphrasing the notebook that built it
would be inventing documentation. Those rows instead carry a derived **assigned in** column
listing every notebook that assigns that name. Its **verification-coverage table is the point**:
32 of 106 modelled features are checked in this fork, and the other 74 are inherited, checked
only in their original `code/predictive_modeling/` notebook. The reproduction gates would not
notice if one of those had always been wrong - which is exactly how the `01h_` Kozak and `01g_`
`csc` bugs survived. Blank in the `checks` column means *not checked here*, not *correct*.

**`feature_annotation.csv` is the central artifact.** One row per column, with `block`
(`id` / `label` / `qc` / `intrinsic` / `external` / `eif4e`), `family` (for exact SHAP
summation, since SHAP is additive), and `variant`. Feature selection is a join against this
table, not a regex cascade. Editing the CSV changes what the model sees.

**Two gates, both asserted, both green.** `01_` checks that the block join reproduces
notebook 07's 51 features exactly; `02_` checks that the resulting importance matches the
frozen baseline CSV (max deviation 3.6e-15). The second runs only under
`variant_set = "current"`, which is why superseded columns are **kept rather than deleted** -
dropping them would make the baseline permanently unreproducible. Reproduction requires
matching 07's **RNG call sequence and column order**, not just its logic: `sample()` must be
the only RNG consumer after `set.seed(neg_seed)`, and ranger sees `train_data[, feature_cols]`,
so reordering columns changes the forest even with identical seeds.

**A gate run sees baseline-era columns only.** Both gates filter `added_in == "baseline"`.
Without that, every family added to the fork afterwards enters the `variant_set = "current"`
run and fails a gate that is about *reproduction*, not about what the fork models today - the
`termination` family did exactly that and left the `02_` gate red until it was fixed. New
feature families therefore need no gate bookkeeping.

**Variant pairs.** A corrected feature is added *alongside* the one it supersedes, never over
it, and `02_`'s `variant_set` param (`"current"` | `"corrected"`, default `"corrected"`) picks
one member of each pair. The two members are near-duplicates and must never enter the same
forest, where they would split credit and hide the effect being measured.

| current | corrected | change |
|---|---|---|
| `kozak_pwm_score` | `kozak_pwm_score_v2` | published Kozak 1987 PWM; see the bug note below |
| `gc3` | `gc3_internal` | initiator and stop codon removed |
| `cai` / `fop` | `cai_internal` / `fop_internal` | initiator AUG removed |
| `csc` | `csc_internal` | iCodon given the internal-codon CDS, not the stop-containing one |
| `tai` | `tai_gtrnadb` | collaborator's GtRNAdb copy number + explicit wobble table |

**`struct_*` is renamed `tco_struct_*` in this fork** (the inventory's one NAME COLLISION row).
The collaborator has `struct_accessibility_*` too, but theirs is a whole-transcript fold while
ours is region-isolated RNAplfold (`-W 150 -L 100`) - same name, different quantity, and
nothing is refolded. `code/predictive_modeling/` is untouched and keeps the old names. Both
gates compare on the baseline's spelling, and the rename has to be applied to `baseline_cols`
as well, or all 18 columns are tagged `added_in = "tco"` and silently dropped from both gates.

**G4mer (`g4` family) is ours alone** - `00_`'s inventory has no collaborator counterpart
(`ours_only`), so unlike `struct_*` there is no name to collide. Three things govern its use:
inference is **already complete** in `output/g4mer/` (not the `g4mer/output/` in `.gitignore`);
**stride differs by region** (5'UTR 10, CDS/3'UTR 20), so the three columns are not comparable
to each other and are never differenced - `01u_` holds the harmonised version for cross-region
work; and `g4mer_max` rises with region length by construction and tracks GC because rG4 needs
G-runs, so the modelled column is `g4mer_max_resid_*`, residualised on log2 region length +
region GC (verified r = 0.000 against both). `g4mer_mean_*` and `g4mer_frac_above_*` are
per-window averages and fractions with no length bias and enter raw.

**`initiation_score` is a separate axis from `variant_set`** (`"both"` (default) | `"pwm"` |
`"noderer"`). `kozak_pwm_score_v2` and `noderer_tis_efficiency` are both start-context
strength - one estimated from Kozak's 1987 frequency table, the other measured by FACS-seq.
`00_`'s inventory says never to model both, on the assumption they are near-duplicates. They
are not: on these transcripts they correlate at **Spearman 0.154**, swapping one for the other
moves test AUC by -0.0001 (p = 0.29, 20 paired seeds), and together they hold 1.60 + 1.42 pp
against the 1.71 pp the PWM holds alone - so they are not splitting credit. Both are kept by
default; the exclusive modes are retained for that comparison. The suffix carries
`initiation_score` whenever it is not `"pwm"`, so a `"both"` run cannot overwrite a `"pwm"` one.

**Always judge a change against the sweep, never against a single run.** Notebook 02 keeps
every positive and draws a size-matched negative class, so one run reports the AUC and ranking
of *one arbitrary draw*. On this model the draw alone moves test AUC by **0.069**
(sd 0.0165, median 0.807). `neg_seed` seeds **only** the draw - split, forest and CV stay
pinned at 9 - so the sweep is **paired**: within a seed both variant sets see identical genes,
and the paired-difference sd (0.0029) is ~6x smaller than the draw-to-draw sd (0.0166). Effects
invisible in a distribution comparison are readable when paired.

**A check that has never been shown to FAIL is uninformative.** `01c_mutation_test_checks.R`
renders `01_` once per named corruption (`mutate_input`) and asserts the matching check fires;
a mutation run never writes the matrix. Currently **6 of 6 caught**. The suite earned its place
immediately: `g4_region_swap` (exchange the 5'UTR and 3'UTR G4mer summaries) originally passed
**every** check, because the motif spot-check only asserted an *ordering* within each region
and G-richness is correlated across regions of the same transcript. The CDS/5'UTR-vs-3'UTR
result rested on an unverified join. The fix is a **discriminative** check - each region's
score must track its *own* motif better than any other region's - and re-running the mutation
confirms it now fires. Add a mutation whenever you add a check.

**Run every new family against its own shuffled null, not against zero.** `02_`'s
`shuffle_families` permutes a family's values across transcripts, destroying its link to the
target while preserving every marginal, the NA pattern and the within-family correlation
structure. Impurity importance is non-negative, so an uninformative family never scores zero -
the shuffled arm is the floor a real one has to clear. Three arms per question
(`exclude_families` / default / `shuffle_families`), 20 seeds each, ~5 min.

**Arms live in `02_sweep_presets.R`, never inline.** `family_arms()` builds the three-arm shape
and gives each arm its own `feature_set_label`, since two arms sharing a label overwrite each
other's outputs. Presets sit side by side: `02b_`'s config list used to be *overwritten* by
each new question, so the termination arms survived only in one commit and the peptide arms in
another, and every number below needed a different checkout to reproduce.

**Importance shares are denominator-dependent.** A family's share falls as unrelated features
are added (this matrix went 57 -> 63 -> 72 features across three additions), so shares from
different builds are not comparable. The table below was re-run end to end on the current
matrix for that reason. The real-vs-shuffled *difference within one sweep* is unaffected, since
both arms see the same matrix.

Results so far on si3d / hypoxia / 1hr, **all of which are per-dataset and must be re-run on a
new mRNA set rather than inherited**:

All three sweeps below were re-run end to end on the **same 72-feature matrix**, so the shares
are comparable to each other.

| family | real vs shuffled share | verdict |
|---|---|---|
| `termination` (11) | 4.53 vs 4.41 pp (+0.08, 16/20, p = 0.005) | at the floor; every per-feature gain within +/-4 ranks |
| `nascent_peptide` (3) | 3.03 vs 2.34 pp (+0.67, 20/20, p < 1e-4) | **real, but carried entirely by `proline_fraction`** (gain **+24** ranks) |
| `polya_track` (2) | 0.497 vs 0.483 pp (+0.03, 12/20, **p = 0.37**) | **null.** Across three matrix versions this went p = 0.003 -> 0.058 -> 0.37 |
| `g4` (9) | 10.24 vs 8.07 pp (**+2.16**, 20/20, p < 1e-4) | **clears the floor by the widest margin**; `g4mer_mean_cds` gains **+26.5** ranks, `g4mer_max_resid_cds` +20.5 |

Re-run after the F8 phantom-zero fix, the effects got **larger**, not smaller: `g4` +1.52 -> +2.16 pp
and `proline_fraction` +12.5 -> +24 ranks. Cleaning the contaminated structure columns stopped
them absorbing credit the real features had earned.

## SHAP noise floor (`03d_shap_seed_sweep.R`)

**One run per experiment is defensible only against a measured floor.** Twenty draws with
*nothing biological varying* - same 1,436 positives every time, a different size-matched
negative draw each time (52% overlap), forest and CV pinned at 9 - give the spread a family's
share shows for free.

| family | median share | sd | 5-95% | min detectable difference |
|---|---|---|---|---|
| `length` | 19.85 | 1.15 | 18.3-21.5 | **3.21** |
| `gc` | 17.28 | 0.82 | 16.3-18.8 | 2.28 |
| `stability_external` | 13.95 | 0.98 | 12.2-15.3 | 2.73 |
| `structure_cds` | 8.36 | 0.74 | 7.2-9.5 | 2.08 |
| `g4` | 7.90 | 0.62 | 7.4-9.0 | 1.73 |
| `codon_optimality` | 7.76 | 0.87 | 6.5-9.1 | 2.45 |
| `clip` | 6.85 | 0.95 | 5.4-8.6 | 2.65 |
| `structure_utr5` | 4.30 | 0.28 | 3.8-4.7 | 0.79 |
| `termination` | 2.46 | 0.30 | 2.2-3.0 | 0.84 |
| `nascent_peptide` | 2.28 | 0.36 | 2.1-3.1 | 1.00 |
| `polya_track` | 0.33 | 0.06 | 0.2-0.4 | 0.17 |

The threshold is **per family, not one number**: `length` needs a 3.2-point difference to mean
anything, `structure_utr5` only 0.8. Two single-run experiments each carry their own draw
noise, so the relevant scale is ~sqrt(2) x sd; the last column is 2.8 x sd. Shares are bounded
at 0 and skewed for the small families, so prefer the empirical 5-95% interval over a
normal-theory one.

**SHAP is more rank-stable than impurity.** 14 of 20 features are top-20 in *every* draw
(impurity: 7 of 20) and the top-20 Jaccard is 0.818 (impurity: 0.67).

**But 14 of 69 features have a `shap_value_cor` sign that FLIPS across draws** - including
`tco_struct_accessibility_stop_proximal_utr3` at exactly 50/50, the three 3'UTR G4 columns and
`term_interval_codons`. Direction is uninterpretable for those, and they are the same features
the shuffled-null calls noise.

**The floor does NOT transfer between experiments** (`03e_experiment_pool_ratios.R`). It is
driven by the pool ratio - negative pool size over positive set size - and across 51 gene sets
that ratio spans **0.96 to 46.5**. The reference experiment sits at 1.94, predicting 51% draw
overlap against 52% measured, so the ratio forecasts the floor without fitting anything. Seven
experiments have a pool *no larger than* their positive set: every seed draws the same genes, a
sweep returns N copies of one model, and it looks reassuringly stable while measuring nothing.
Most of those are si3e sets scored against the si3d pool and need their own comparator.

**The per-feature `gain` (shuffled rank - real rank) is the readable summary**, and `02c_`
section 4 prints it. Positive means the real values earn the place; zero or negative means the
column is indistinguishable from noise.

| earns its place | gain | indistinguishable from noise | gain |
|---|---|---|---|
| `g4mer_mean_cds` | +15.5 | `max_net_charge_30aa`, `ppp_motif_density`, `max_consecutive_aaa`/`aag` | **0.0** each |
| `proline_fraction` | +12.5 | every `term_*` column | -5 to +3.5 |
| `g4mer_mean_utr5` | +9.0 | `g4mer_max_resid_utr3` | **-7.5** |
| `g4mer_max_resid_cds` | +8.5 | `g4mer_mean_utr3` | **-12.0** |
| `g4mer_max_resid_utr5` | +6.0 | | |

Four of the five nascent-peptide/poly(A) columns have a gain of **exactly 0.0** - identical
median rank to their own permuted copies. Both 3'UTR G4 columns score *worse* than their
shuffled versions.

`max_net_charge_30aa`, `ppp_motif_density`, `max_consecutive_aag` and `max_consecutive_aaa`
are **rank-identical to their shuffled copies** (44/46/48/49 of 63). No family beats the
*no-family* arm on AUC: `g4` recovers the dilution cost that 9 noise columns impose
(`shuffled - none` = -0.0015, p = 0.03) without adding beyond it.

**Within `g4` the regions separate, and that is the result.** `g4mer_mean_cds` (rank 21 vs 36.5
shuffled), `g4mer_max_resid_cds` (25 vs 33.5), `g4mer_mean_utr5` (28.5 vs 37.5) and
`g4mer_max_resid_utr5` (31 vs 37) sit well above their nulls; both 3'UTR columns sit *below*
theirs (44/45.5 vs 36.5/33.5). CDS and 5'UTR are what a scanning or elongating ribosome
traverses. This is measured on residualised columns, so it is not GC in disguise.
`g4mer_frac_above_*` is bottom-ranked in both arms in all three regions (median 0 across the
transcriptome), as is `ppp_motif_density`. `proline_fraction` carries genuine signal that the forest can mostly get
elsewhere: it correlates with `cds_gc` at rho = 0.44 (proline codons are CCN) and `cds_gc` is
already the model's strongest single feature. **All of these features are deliberately kept** -
they are validated and cheap, and the fork exists to run on other mRNA sets where they may
matter.

**Independent audit findings, all fixed in the fork (`01_`).** An independent agent recomputed
the matrix from primary sources and found seven real defects that every existing check passed.
Each is now also a standing structural invariant in `01_`, because each was a *class* of bug:

| # | defect | fix |
|---|---|---|
| F6 | six `tco_struct_*` CDS/3'UTR columns filed under `structure_utr5` - `assign_family()` anchored the region tag with `$`, but those names end in `_mean`/`_min`/`_sd` | regex now matches the tag anywhere; asserted 6/6/6 |
| F1 | 1,689 transcripts had every structure column **zero-filled instead of NA** (`replace_na(., 0)` in `01d_`/`01k_`), so 0 read as "maximally structured" and every NA-based guard downstream missed them | sentinel (mean == 0 **and** sd == 0) restored to NA |
| F2 | `tco_struct_num_structured_regions*` is region length in disguise (Spearman ~1.000, and the "count" exceeds the length) | deprecated |
| F5 | `tco_struct_accessibility_*_min` identically 0 everywhere - `min()` taken after NA became 0 | deprecated; `02_` also drops zero-variance columns generally |
| F3 | `hia2026_dhx29_occupancy` joined at **gene** level from a transcript-keyed source, so ~1,013 transcripts carry another isoform's occupancy | `hia2026_dhx29_occupancy_tx` added as the corrected variant |
| F4 | `clip_lee_complete_media_bound` is bit-identical to `clip_lee_glu_dep_bound` (already identical upstream), so `clip_cap_binding_count` counted one dataset twice | duplicate deprecated; `clip_cap_binding_count_v2` counts the 2 distinct conditions |
| F7 | `te_lfc` / `te_padj` / `te_lfc_bin` are the **siCTRL hypoxia-vs-normoxia** contrast, not a knockdown contrast (bit-for-bit `translation_categories_1hr.csv`, max dev 0) | no model impact (`block = label`, asserted out), but the identity is now pinned by an assertion and stated in the definition |

**The repairs make bit-exact baseline reproduction impossible, deliberately.** F1 changes values,
so `01_` takes `apply_audit_repairs` (default TRUE): build with it **FALSE** to reproduce the
historical matrix, then run `02_`'s gate. `02_` skips its zero-variance drop on a gate run for
the same reason. Both gates are green on that path.

**Consequences for results already reported.** F6 means any `structure_utr5` number from before
the fix was a mixture of all three regions: corrected, the split is CDS 5.78%, 5'UTR 2.62%,
3'UTR 1.49% of SHAP. The model is 66 features after the deprecations.

## Findings that affect the EXISTING pipeline

These are defects in `code/predictive_modeling/`, found by the reconciliation above. They are
fixed only in the fork; the original notebooks are unchanged.

**`kozak_pwm_score` is wrong (`01h_`).** The hand-typed Kozak 1987 matrix has `G = 0.13` at
position -3; the published table and our own transcripts both give 0.36-0.41. Position -3 is
the **purine** position (the R in GCCRCCATGG), so the PWM penalises G by 0.93 bits while
`match_kozak()` in the same notebook scores `[AG]` as equivalent - the two Kozak features
contradict each other for ~40% of transcripts. Position -6 is also wrong (argmax T, should be
G). **The existing tier-monotonicity check cannot detect this** and passes under both matrices,
because tiers are set mostly by +4 and by `[CT]`-vs-`[AG]`, so a within-tier error cancels in
the tier means. Any PWM/tier feature needs a **within-tier** check. Practical impact on the
si3d hypoxia 1hr model: none, because the feature ranks 31 of 51.

**A phantom zero contaminates every `struct_*` column (`01d_`, `01k_`).** Both notebooks find
the `_lunp` header with `str_starts(lines, "#")`, but the second line of an RNAplfold `_lunp`
file is `" #i$\tl=1\t..."` - it begins with a **space**. The test misses it, the column-header
line is passed to `read.table` as data, it parses to an all-NA row, and
`acc_vals[is.na(acc_vals)] <- 0` turns that into **0**, which means *maximally structured*. So
every transcript carries a phantom position at the 5' end with accessibility zero.

One bug, three symptoms:

| symptom | why |
|---|---|
| `struct_accessibility_*_min` is identically 0 for every transcript | the phantom zero is always the minimum |
| `struct_num_structured_regions*` equals region length **+ 1** | the phantom row is the +1 (confirmed for 7,653 of 8,572 CDS) |
| every accessibility mean is diluted and every sd inflated | one extra zero in the vector |

The positional features are the worst affected: `cap_proximal`, `start_proximal_cds` and
`stop_proximal_utr3` average exactly 30 positions, so the phantom is **3.3% of the window
regardless of transcript length**. Measured on the 5'UTR, correcting it moves
`accessibility_cap_proximal` by **+0.0164** and `accessibility_utr5_mean` by **+0.0038**
(medians) - the feature you would interrogate for 5'-end biology was distorted four times as
much as the aggregate.

**`code/predictive_modeling/` is deliberately NOT fixed**, per the freeze: every `struct_*`
value there, and anything downstream that consumed them, carries the phantom. The fork
recomputes them correctly from the same stored `_lunp` files
(`01f_recompute_structure.R`), which reproduces the buggy values bit-for-bit first as proof
the defect is characterised. Treat any published `struct_*` number from the frozen pipeline
with this in mind.

**`csc` partly encodes stop-codon identity (`01g_`).** iCodon is given the stop-containing CDS.
Removing the stop gives r = 0.939 and shifts the median 0.116 -> -0.015. Impact is **real and
verified over 20 draws**: `csc` falls from median rank 3.5 to 10.0 with **non-overlapping**
rank ranges (2-6 vs 7-12) and a -1.35 pp share drop that is negative in 20/20 seeds. `tai`
falls similarly (-0.30 pp, 20/20). Test AUC is unaffected. CAI and FOP already excluded stops
and barely move (r ~ 0.999). Treat published `csc` rankings with this in mind.

**118 transcripts carry truncated CDS.** GENCODE `cds_end_NF`: the CDS starts at ATG but the
annotation runs out, so there is no terminal stop. 65 are not divisible by 3 and `01g_`'s
`%%3` filter already drops them, but **53 are in-frame and pass every existing filter**,
receiving full `cai`/`tai`/`csc`/`fop` on an incomplete CDS (their `gc3` median is 48.9 vs
57.2). The right mask is a **missing terminal stop**, not a bad frame. A separate 25
transcripts start with a non-AUG codon (24 tagged `non_ATG_start`, 17 MANE Select) - these are
**real biology, keep them**; only their Kozak score is meaningless.

**Never transcribe reference values.** Import the file and parse constants from source, then
assert the parse reproduces the already-saved feature column. `00_` does this for the `01h_`
matrix (deviation exactly 0 over 10,140 transcripts) and `01_` does it for CAI/FOP against
`01g_` (max deviation 0). Collaborator reference tables are mirrored byte-for-byte in
`accessories/collaborator_reference/`.

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
- **Kozak within-tier purine check**: inside tier 2, mean PWM for A-at-−3 vs G-at-−3 must
  differ by < 0.75 bits. Between-tier monotonicity passes even on a wrong matrix; this is
  the check that catches it
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
| `output/g4mer/g4mer_summary_{utr5,cds,utr3}.tsv` | G4mer rG4 scores per transcript (`g4mer_max`, `_mean`, `_frac_above`, `_argmax_start`). **Inference is already done for all three regions** (9.5-10k transcripts each) - note the path is `output/g4mer/`, not the `g4mer/output/` in `.gitignore`. The 4.5 h GPU figures in `g4mer/README.md` are for re-running at stride 1, which is only needed for exact reproduction of published values. Loaded and range-checked by `code/predictive_modeling/01u_feature_extraction_g4mer.Rmd` |
