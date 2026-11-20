# Hypoxia-specific eIF3d-promoted mRNAs, consistent in both cell lines.
# A gene qualifies if, in BOTH MDA-MB-231 and MCF7-SIX1:
#   hypoxia 1hr : te_padj < PADJ and te_lfc >= MIN_LFC  (same floor as 04_)
#   normoxia    : te_padj >= PADJ (not significant; NA = not tested -> excluded, counted)
# Strict tier additionally requires |normoxia te_lfc| < NOR_LFC (the project's 0.5 standard),
# so "not significant" cannot be a large but noisy normoxia effect.
# Ranked by mean hypoxia TE LFC across the two lines.
#
# Contrasts: 231 hypoxia 1hr, 231 normoxia 1+4hr time-regressed; MCF7-SIX1 hypoxia/normoxia 1hr
# (its only timepoint). All LFCs log2(sictrl / si3d): positive = eIF3d promotes.
# TE is gene-level: a hit still needs the isoform / 5'UTR check (02_) before reporter use.
#
# Usage: Rscript code/finding_positive_controls/06_hypoxia_specific_eif3d_targets.R [min_lfc]  (default 1)
# Writes output/finding_positive_controls/06_hypoxia_specific_eif3d_targets_lfc{min_lfc}.csv
suppressPackageStartupMessages({ library(tidyverse); library(here) })
args <- commandArgs(trailingOnly = TRUE)
PADJ <- 0.05; NOR_LFC <- 0.5
MIN_LFC <- if (length(args)) as.numeric(args[1]) else 1  # hypoxia TE LFC floor

read_te <- function(f) read_csv(
  here("output", paste0("translation_categories_si3d_vs_sictrl_", f, ".csv")),
  show_col_types = FALSE) %>% mutate(gene_id_clean = sub("\\..*", "", gene_id))
te <- map(c(h231 = "hypoxia_1hr", n231 = "normoxia_1and4hr",
            hmcf = "mcf7six1_hypoxia_1hr", nmcf = "mcf7six1_normoxia_1hr"), read_te)

# Direction checks (as in 04_)
for (k in c("h231", "n231")) {
  e <- te[[k]]$rna_lfc[te[[k]]$symbol == "EIF3D"]
  stopifnot("231 direction flipped" = length(e) == 1 && e > 1.5)
}
for (k in c("hmcf", "nmcf")) {
  d <- te[[k]]
  stopifnot("MCF7-SIX1 te_lfc is not the negated source LFC" =
              isTRUE(all.equal(d$te_lfc, -d$source_log2FoldChange)),
            "MCF7-SIX1 direction check failed" = d$te_lfc[d$symbol == "EIF3D"] < -0.5)
}
walk(te, ~ stopifnot(!anyDuplicated(.x$gene_id_clean)))

pick <- function(k, s) te[[k]] %>% dplyr::select(gene_id_clean, !!paste0("lfc_", s) := te_lfc,
                                                  !!paste0("padj_", s) := te_padj)
d <- te$h231 %>% dplyr::select(gene_id_clean, symbol) %>%
  inner_join(pick("h231", "hyp_231"), by = "gene_id_clean") %>%
  inner_join(pick("n231", "nor_231"), by = "gene_id_clean") %>%
  inner_join(pick("hmcf", "hyp_mcf"), by = "gene_id_clean") %>%
  inner_join(pick("nmcf", "nor_mcf"), by = "gene_id_clean")
cat("genes in all four contrasts:", nrow(d), "\n")

hyp_hit <- d %>% filter(padj_hyp_231 < PADJ, padj_hyp_mcf < PADJ,
                        lfc_hyp_231 >= MIN_LFC, lfc_hyp_mcf >= MIN_LFC)
cat("hypoxia hits in both lines (padj <", PADJ, ", TE LFC >=", MIN_LFC, "):", nrow(hyp_hit), "\n")
cat("  of which normoxia padj NA in either line (not tested, excluded):",
    sum(is.na(hyp_hit$padj_nor_231) | is.na(hyp_hit$padj_nor_mcf)), "\n")

cat("  normoxia significant in 231 only / MCF7-SIX1 only / both:",
    sum(hyp_hit$padj_nor_231 < PADJ & !(hyp_hit$padj_nor_mcf < PADJ), na.rm = TRUE),
    sum(!(hyp_hit$padj_nor_231 < PADJ) & hyp_hit$padj_nor_mcf < PADJ, na.rm = TRUE),
    sum(hyp_hit$padj_nor_231 < PADJ & hyp_hit$padj_nor_mcf < PADJ, na.rm = TRUE), "\n")
out <- hyp_hit %>%
  filter(padj_nor_231 >= PADJ, padj_nor_mcf >= PADJ) %>%
  mutate(mean_hyp_lfc = (lfc_hyp_231 + lfc_hyp_mcf) / 2,
         mean_nor_lfc = (lfc_nor_231 + lfc_nor_mcf) / 2,
         hyp_minus_nor = mean_hyp_lfc - mean_nor_lfc,
         strict = abs(lfc_nor_231) < NOR_LFC & abs(lfc_nor_mcf) < NOR_LFC) %>%
  arrange(desc(mean_hyp_lfc)) %>% mutate(rank = row_number()) %>%
  relocate(rank, symbol, strict, mean_hyp_lfc, mean_nor_lfc, hyp_minus_nor)
cat("... and normoxia NOT significant in both lines:", nrow(out), "\n")
cat("... and |normoxia TE LFC| <", NOR_LFC, "in both (strict):", sum(out$strict), "\n")
stopifnot(all(out$lfc_hyp_231 >= MIN_LFC), all(out$padj_nor_mcf >= PADJ))

out_dir <- here("output", "finding_positive_controls")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(out, file.path(out_dir, paste0("06_hypoxia_specific_eif3d_targets_lfc", MIN_LFC, ".csv")))
print(out %>% mutate(across(where(is.double), ~ signif(.x, 2))) %>%
        dplyr::select(-gene_id_clean), n = 40, width = 200)
