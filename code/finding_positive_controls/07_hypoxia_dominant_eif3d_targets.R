# Hypoxia-dominant eIF3d-promoted mRNAs, consistent in both cell lines.
# Relaxed version of 06_ (Kate, 2026-11-20): normoxia need not be flat, only much weaker.
# A gene qualifies if, in EACH of MDA-MB-231 and MCF7-SIX1 separately:
#   hypoxia 1hr : te_padj < PADJ and te_lfc >= MIN_LFC  (same floor as 04_/06_)
#   normoxia    : hypoxia te_lfc >= RATIO * normoxia te_lfc  (signed; no padj requirement)
# A negative normoxia LFC (eIF3d inhibits in normoxia) passes the ratio by definition; those
# genes are flagged `nor_opposite` because they are a direction switch, not a weaker effect.
# Ranked by RANK_BY: "mean" = mean hypoxia TE LFC across the two lines (default);
# "hyp_231" = MDA-MB-231 hypoxia TE LFC (Kate, 2026-11-20: large effect size in 231).
#
# Contrasts: 231 hypoxia 1hr, 231 normoxia 1+4hr time-regressed; MCF7-SIX1 hypoxia/normoxia 1hr
# (its only timepoint, files already negated - asserted). All LFCs log2(sictrl / si3d):
# positive = eIF3d promotes. TE is gene-level: hits still need the isoform / 5'UTR check (02_).
#
# Usage: Rscript code/finding_positive_controls/07_hypoxia_dominant_eif3d_targets.R [min_lfc] [ratio] [rank_by]
#   defaults 1, 2, mean
# Writes output/finding_positive_controls/07_hypoxia_dominant_eif3d_targets_lfc{min_lfc}_ratio{ratio}[_rank{rank_by}].csv (suffix only when not "mean")
suppressPackageStartupMessages({ library(tidyverse); library(here) })
args <- commandArgs(trailingOnly = TRUE)
PADJ <- 0.05
MIN_LFC <- if (length(args) >= 1) as.numeric(args[1]) else 1
RATIO   <- if (length(args) >= 2) as.numeric(args[2]) else 2
RANK_BY <- if (length(args) >= 3) args[3] else "mean"
stopifnot(RANK_BY %in% c("mean", "hyp_231"))
stopifnot(MIN_LFC > 0, RATIO >= 1)

read_te <- function(f) read_csv(
  here("output", paste0("translation_categories_si3d_vs_sictrl_", f, ".csv")),
  show_col_types = FALSE) %>% mutate(gene_id_clean = sub("\\..*", "", gene_id))
te <- map(c(h231 = "hypoxia_1hr", n231 = "normoxia_1and4hr",
            hmcf = "mcf7six1_hypoxia_1hr", nmcf = "mcf7six1_normoxia_1hr"), read_te)

# Direction checks (as in 04_/06_)
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
cat("  normoxia LFC missing in either line (excluded):",
    sum(is.na(hyp_hit$lfc_nor_231) | is.na(hyp_hit$lfc_nor_mcf)), "\n")

pass_231 <- hyp_hit$lfc_hyp_231 >= RATIO * hyp_hit$lfc_nor_231
pass_mcf <- hyp_hit$lfc_hyp_mcf >= RATIO * hyp_hit$lfc_nor_mcf
cat("  hypoxia >=", RATIO, "x normoxia in 231 only / MCF7-SIX1 only / both / neither:",
    sum(pass_231 & !pass_mcf, na.rm = TRUE), sum(!pass_231 & pass_mcf, na.rm = TRUE),
    sum(pass_231 & pass_mcf, na.rm = TRUE), sum(!pass_231 & !pass_mcf, na.rm = TRUE), "\n")

out <- hyp_hit %>%
  filter(lfc_hyp_231 >= RATIO * lfc_nor_231, lfc_hyp_mcf >= RATIO * lfc_nor_mcf) %>%
  mutate(mean_hyp_lfc = (lfc_hyp_231 + lfc_hyp_mcf) / 2,
         mean_nor_lfc = (lfc_nor_231 + lfc_nor_mcf) / 2,
         # hypoxia / normoxia per line; NA where normoxia <= 0 (ratio undefined / sign switch)
         ratio_231 = if_else(lfc_nor_231 > 0, lfc_hyp_231 / lfc_nor_231, NA_real_),
         ratio_mcf = if_else(lfc_nor_mcf > 0, lfc_hyp_mcf / lfc_nor_mcf, NA_real_),
         nor_opposite = lfc_nor_231 < 0 | lfc_nor_mcf < 0,
         nor_sig_either = (padj_nor_231 < PADJ) %in% TRUE | (padj_nor_mcf < PADJ) %in% TRUE) %>%
  arrange(desc(if (RANK_BY == "mean") mean_hyp_lfc else lfc_hyp_231)) %>% mutate(rank = row_number()) %>%
  relocate(rank, symbol, mean_hyp_lfc, mean_nor_lfc, ratio_231, ratio_mcf,
           nor_opposite, nor_sig_either)
cat("... pass in both lines:", nrow(out), "\n")
cat("    of which normoxia LFC negative in a line (direction switch):", sum(out$nor_opposite), "\n")
cat("    of which normoxia still significant (padj <", PADJ, ") in a line:", sum(out$nor_sig_either), "\n")
stopifnot(all(out$lfc_hyp_231 >= RATIO * out$lfc_nor_231),
          all(out$lfc_hyp_mcf >= RATIO * out$lfc_nor_mcf),
          all(out$lfc_hyp_231 >= MIN_LFC), all(out$lfc_hyp_mcf >= MIN_LFC))

out_dir <- here("output", "finding_positive_controls")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(out, file.path(out_dir, sprintf("07_hypoxia_dominant_eif3d_targets_lfc%s_ratio%s%s.csv",
                                          MIN_LFC, RATIO, if (RANK_BY == "mean") "" else paste0("_rank", RANK_BY))))
print(out %>% transmute(rank, symbol, hyp_231 = lfc_hyp_231, nor_231 = lfc_nor_231,
                        hyp_mcf = lfc_hyp_mcf, nor_mcf = lfc_nor_mcf, padj_nor_231, padj_nor_mcf,
                        nor_opposite) %>%
        mutate(across(where(is.double), ~ signif(.x, 2))), n = 40, width = 200)
