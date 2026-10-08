# 10y_sictrl_hypoxia_outcome_prep.R
# MDA-MB-231 siCTRL hypoxia vs normoxia TE change at 1, 4 and 24hr as regression outcomes, so the
# cells' own hypoxic TE response can be modelled and SHAP-clustered exactly as the si3d contrasts
# are (requested by Kate, 2026-10-04). No knockdown: run with outcome_kind = "te_condition".
#
# Source: output/translation_categories_{1,4,24}hr.csv, written by code/rna_ribo_te_analysis.Rmd
# (deFunction, control_condition = "normoxia", test_condition = "hypoxia", sictrl, spike-in
# normalised). Nothing is refitted; this re-checks the direction and writes the outcome tables.
# Orientation: log2(hypoxia / normoxia), positive = higher TE in hypoxia.
#
# Spike-in normalisation makes TE fall across the whole transcriptome by 4hr (median about -0.8),
# so te_padj < 0.05 marks almost every gene at 4/24hr and "te_lfc > 0.5" almost none. The spot-check
# / annotation set is therefore defined RELATIVE to the transcriptome: te_lfc minus the median over
# the modelled genes > 0.5 ("resists the hypoxic TE drop"). No padj filter, because padj tests
# against zero, not against the global shift. The set only labels genes and checks the join; the
# regression and the clustering use every gene. The forest is unaffected by the global shift
# (a constant offset).
#
# rna_lfc is deliberately NOT written to the outcome table: 10_ would read it as a knockdown ruler.
#   Rscript code/translational_control_overview/10y_sictrl_hypoxia_outcome_prep.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here) })
fm  <- readRDS(here("output", "translational_control_overview", "feature_matrix_tco.rds"))
ann <- read_tsv(here("accessories", "human", "gene_anno_hs_dm_v49_r111.tsv"), col_names = FALSE,
                show_col_types = FALSE) %>%
  transmute(gene_id_clean = sub("[.].*", "", X1), symbol = X2) %>% distinct(gene_id_clean, .keep_all = TRUE)

for (t in c("1hr", "4hr", "24hr")) {
  d <- read_csv(here("output", paste0("translation_categories_", t, ".csv")), show_col_types = FALSE) %>%
    mutate(gene_id_clean = sub("[.].*", "", gene_id)) %>% left_join(ann, by = "gene_id_clean")
  stopifnot("duplicate gene_id" = !any(duplicated(d$gene_id_clean)),
            "non-finite te_lfc" = all(is.finite(d$te_lfc[!is.na(d$te_lfc)])))
  # Direction: HIF targets' RNA must rise under hypoxia. At 1hr RNA has barely moved (shrunken
  # rna_lfc IQR ~0.001-0.013), so the ruler is asserted at 24hr and only reported earlier.
  hif <- d %>% filter(symbol %in% c("CA9", "VEGFA", "EGLN3"))
  cat(sprintf("%s: HIF-target rna_lfc %s\n", t, paste(sprintf("%s %+.2f", hif$symbol, hif$rna_lfc), collapse = ", ")))
  if (t == "24hr") stopifnot("hypoxia contrast direction flipped: HIF targets' RNA not up at 24hr" =
                               nrow(hif) == 3 && all(hif$rna_lfc > 0.5))
  # The matrix's te_lfc label column IS this 1hr contrast (CLAUDE.md audit finding F7): pin it.
  if (t == "1hr") {
    m <- match(fm$gene_id_clean, d$gene_id_clean)
    dev <- max(abs(fm$te_lfc - d$te_lfc[m]), na.rm = TRUE)
    cat(sprintf("1hr te_lfc vs the matrix's te_lfc label: max deviation %g\n", dev))
    stopifnot("1hr outcome does not reproduce the matrix's te_lfc (F7)" = dev == 0)
  }
  out <- d %>% filter(!is.na(te_lfc)) %>% transmute(gene_id, symbol, te_lfc, te_padj)
  write_csv(out, here("output", paste0("sictrl_hypoxia_vs_normoxia_te_", t, ".csv")))
  modelled <- out %>% filter(sub("[.].*", "", gene_id) %in% fm$gene_id_clean)
  med <- median(modelled$te_lfc)
  pos <- modelled %>% filter(te_lfc - med > 0.5)
  write_csv(pos %>% transmute(db_gene_symbol = symbol, ensembl_gene = sub("[.].*", "", gene_id), te_lfc,
                              gs_name = paste0("sictrl_hypoxia_", t, "_resists_te_drop"),
                              gs_description = "siCTRL TE log2(hypoxia/normoxia) minus the median over modelled genes > 0.5; spot-check / annotation set"),
            here("output", "genesets", paste0("hypvsnor_ctrl_promotes_TE_", t, "_lfc0.5.csv")))
  cat(sprintf("%s: %d genes (%d modelled) | median TE %+.3f over modelled genes | resists-drop set: %d\n",
              t, nrow(out), nrow(modelled), med, nrow(pos)))
}
