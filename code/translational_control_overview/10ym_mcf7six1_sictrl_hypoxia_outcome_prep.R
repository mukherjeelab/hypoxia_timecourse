# 10ym_mcf7six1_sictrl_hypoxia_outcome_prep.R
# MCF7-SIX1 siCTRL hypoxia vs normoxia TE change at 1hr as a regression outcome - the MCF7-SIX1
# counterpart of 10y_sictrl_hypoxia_outcome_prep.R (MDA-MB-231, 1/4/24hr), so the cells' own
# hypoxic TE response can be modelled in both cell lines (requested by Kate, 2026-11-18).
# Run 10_/11_ with outcome_kind = "te_condition", cell_line = "MCF7-SIX1", gene = "ctrl",
# condition = "mcf7six1_hypvsnor".
#
# Source: counts/all_norm_to_hyp_changes.csv in the 2024_eIF3e_hypoxia repo, column LFC_te_sictrl,
# written by code/5_hypoxicresponse.Rmd (interaction coefficient, normoxia and RNA as reference
# levels, lfcShrink-shrunk, NOT spike-in normalised). Only 1hr exists in MCF7-SIX1.
# Orientation: log2 TE(hypoxia / normoxia), positive = higher TE in hypoxia. Nothing is negated -
# unlike the si3d/si3e MCF7-SIX1 files (10m_), this contrast already has normoxia as reference.
#
# Direction is asserted two ways (the same rulers as 45_hypoxia_te_grids_sictrl.Rmd):
#   1. the producing source is parsed for the reference levels and the coefficient;
#   2. the TE change correlates positively with MDA-MB-231's siCTRL 1hr hypoxia response.
# The HIF-target RNA ruler is not usable here: MCF7-SIX1 shows no early RNA induction of them.
#
# Gene sets, defined as for the MDA-MB-231 runs (relative to the median over modelled genes,
# since the test is against the global shift, not zero). They label genes and drive 10_'s
# spot-check only; the regression uses every gene.
#   positive: te_lfc - median > 0.5                      ("resists the hypoxic TE drop")
#   negative: |te_lfc - median| < 0.5 and te_padj > 0.3  (unchanged relative to the median)
#
#   Rscript code/translational_control_overview/10ym_mcf7six1_sictrl_hypoxia_outcome_prep.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here) })
eif3_repo <- "/Users/katematlin/github/2024_eIF3e_hypoxia/2024_eIF3e_hypoxia"
O <- here("output", "translational_control_overview")

# ---- 1. Direction ruler: parse the producing source ------------------------------------------
src <- readLines(file.path(eif3_repo, "code", "5_hypoxicresponse.Rmd"))
stopifnot(
  "MCF7-SIX1 condition reference is not normoxia" =
    any(grepl('metadata_acute_sictrl\\$condition <- relevel\\(factor\\(metadata_acute_sictrl\\$condition\\), ref = "normoxia"\\)', src)),
  "MCF7-SIX1 library reference is not rna" =
    any(grepl('metadata_acute_sictrl\\$library <- relevel\\(factor\\(metadata_acute_sictrl\\$library\\), ref = "rna"\\)', src)),
  "MCF7-SIX1 coefficient is not the interaction" =
    any(grepl("results\\(ddsMat_acute_sictrl, name=resultsNames\\(ddsMat_acute_sictrl\\)\\[4\\]", src)),
  "LFC_te_sictrl is not taken from res_acute_sictrl" =
    any(grepl("sictrl_te <- res_acute_sictrl", src)) && any(grepl('colnames\\(sictrl_te\\)\\[3\\] ="LFC_te_sictrl"', src)))
cat("source parse: normoxia and RNA are the reference levels; LFC_te_sictrl = hypoxia x ribo interaction\n")

# ---- 2. Load ------------------------------------------------------------------------------------
d <- read_csv(file.path(eif3_repo, "counts", "all_norm_to_hyp_changes.csv"), show_col_types = FALSE) %>%
  transmute(gene_id, symbol, gene_id_clean = sub("[.].*", "", gene_id),
            te_lfc = LFC_te_sictrl, te_padj = padj_te_sictrl)
stopifnot("duplicate gene_id" = !anyDuplicated(d$gene_id_clean),
          "non-finite te_lfc" = all(is.finite(d$te_lfc[!is.na(d$te_lfc)])),
          "padj outside [0, 1]" = all(d$te_padj >= 0 & d$te_padj <= 1, na.rm = TRUE))

# ---- 3. Direction ruler: agreement with MDA-MB-231's own hypoxic TE response ------------------
ref <- read_csv(here("output", "translation_categories_1hr.csv"), show_col_types = FALSE) %>%
  transmute(gene_id_clean = sub("[.].*", "", gene_id), te_231 = te_lfc)
j <- inner_join(d, ref, by = "gene_id_clean") %>% filter(!is.na(te_lfc), !is.na(te_231))
rho <- cor(j$te_lfc, j$te_231, method = "spearman")
cat(sprintf("Spearman vs MDA-MB-231 siCTRL 1hr hypoxia TE change: %+.3f (n = %d); negated it would be %+.3f\n",
            rho, nrow(j), -rho))
stopifnot("MCF7-SIX1 and MDA-MB-231 hypoxic TE responses disagree in direction" = rho > 0.2)

# ---- 4. Write the outcome and the gene sets -----------------------------------------------------
fm <- readRDS(file.path(O, "feature_matrix_tco_mcf7six1.rds"))
out <- d %>% filter(!is.na(te_lfc)) %>% dplyr::select(gene_id, symbol, te_lfc, te_padj)
write_csv(out, here("output", "mcf7six1_sictrl_hypoxia_vs_normoxia_te_1hr.csv"))

modelled <- d %>% filter(!is.na(te_lfc), gene_id_clean %in% fm$gene_id_clean)
med <- median(modelled$te_lfc)
pos <- modelled %>% filter(te_lfc - med > 0.5)
neg <- modelled %>% filter(abs(te_lfc - med) < 0.5, !is.na(te_padj), te_padj > 0.3)
stopifnot("positive and negative sets overlap" = !length(intersect(pos$gene_id_clean, neg$gene_id_clean)),
          "positive set implausibly small" = nrow(pos) >= 50,
          "negative pool does not clear the positive set" = nrow(neg) > nrow(pos))

write_csv(pos %>% transmute(db_gene_symbol = symbol, ensembl_gene = gene_id_clean, te_lfc,
                            gs_name = "mcf7six1_sictrl_hypoxia_1hr_resists_te_drop",
                            gs_description = "MCF7-SIX1 siCTRL TE log2(hypoxia/normoxia) minus the median over modelled genes > 0.5; spot-check / annotation set"),
          here("output", "genesets", "mcf7six1_hypvsnor_ctrl_promotes_TE_1hr_lfc0.5.csv"))
write_csv(neg %>% dplyr::select(gene_id, symbol, te_lfc, te_padj),
          here("output", "genesets", "mcf7six1_hypvsnor_ctrl_negative_controls_lfc0.5.csv"))

cat(sprintf("1hr: %d genes (%d modelled) | median TE %+.3f over modelled genes | sd %.3f\n",
            nrow(out), nrow(modelled), med, sd(modelled$te_lfc)))
cat(sprintf("resists-drop set: %d | negative pool: %d (pool ratio %.1f)\n",
            nrow(pos), nrow(neg), nrow(neg) / nrow(pos)))
