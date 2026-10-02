# 10t_teleman_outcome_prep.R
# Teleman / Roiuk HeLa ribosome profiling (GSE243708) as two regression outcomes for the
# HeLa-isoform matrix: ABSOLUTE translation efficiency, log2(ribo / RNA), within
#   ev - empty vector (eIF4E active)
#   bp - 4E-BP1(4A) overexpression (eIF4E inhibited; genetic, not a drug)
# Chosen by Kate on 2026-09-30 (absolute TE in both conditions, intrinsic features only).
#
# Source: output/predictive_modeling/teleman_absolute_te.rds, written by
# code/predictive_modeling/97_absolute_te_by_4e_class.Rmd (DESeq2 group contrast ribo vs RNA within
# each condition, ashr shrunk). This script does not refit anything; it re-asserts 97's own
# validation on the saved object, then writes the outcome tables.
#
# These are TE LEVELS, not a knockdown change: positive = efficiently translated in that
# condition. Run with outcome_kind = "te_level" in 10_rf_regression.Rmd.
#   Rscript code/translational_control_overview/10t_teleman_outcome_prep.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here) })
PD <- here("output", "predictive_modeling")
fit <- readRDS(file.path(PD, "teleman_absolute_te.rds"))
te_all <- readRDS(file.path(PD, "teleman_riboseq_te_fixed.rds"))$fixed   # nb76 TE change
stopifnot("absolute-TE object is not the expected list" =
            all(c("vehicle", "inhibitor", "norm", "meta") %in% names(fit)))

# --- 97's checks, re-asserted -----------------------------------------------------------------
grp_cols <- function(g) rownames(fit$meta)[fit$meta$group == g]
naive <- tibble(gene_id  = rownames(fit$norm),
  ev_naive = log2((rowMeans(fit$norm[, grp_cols("vehicle_ribo"),  drop = FALSE]) + 1) /
                  (rowMeans(fit$norm[, grp_cols("vehicle_rna"),   drop = FALSE]) + 1)),
  bp_naive = log2((rowMeans(fit$norm[, grp_cols("BP_inhib_ribo"), drop = FALSE]) + 1) /
                  (rowMeans(fit$norm[, grp_cols("BP_inhib_rna"),  drop = FALSE]) + 1)))
te <- fit$vehicle %>% dplyr::select(gene_id, ev = log2FoldChange, padj_ev = padj, baseMean) %>%
  inner_join(fit$inhibitor %>% dplyr::select(gene_id, bp = log2FoldChange, padj_bp = padj), by = "gene_id") %>%
  inner_join(naive, by = "gene_id") %>%
  inner_join(te_all %>% dplyr::select(gene_id, symbol, nb76 = log2FoldChange), by = "gene_id") %>%
  mutate(delta = bp - ev)
r_ev <- cor(te$ev, te$ev_naive); r_bp <- cor(te$bp, te$bp_naive); r_anchor <- cor(te$delta, te$nb76)
cat(sprintf("DESeq vs naive TE: E.V. r %.3f | 4E-BP1 r %.3f | TE(bp) - TE(ev) vs nb76 TE change r %.3f | n %d\n",
            r_ev, r_bp, r_anchor, nrow(te)))
s <- function(g, k) te[[k]][match(g, te$symbol)]
stopifnot("DESeq and naive TE disagree in E.V." = r_ev > 0.85,
          "DESeq and naive TE disagree in 4E-BP1(4A)" = r_bp > 0.85,
          "TE levels do not reproduce nb76's TE change - contrast may be inverted" = r_anchor > 0.9,
          "CDKN2B TE does not rise under eIF4E inhibition" = s("CDKN2B", "bp") > s("CDKN2B", "ev"),
          "GLO1 TE does not fall under eIF4E inhibition" = s("GLO1", "bp") < s("GLO1", "ev"),
          "non-finite TE" = all(is.finite(te$ev)) && all(is.finite(te$bp)),
          "duplicate gene_id" = !any(duplicated(te$gene_id)))
cat("97's validation re-asserted: all passed\n")

# --- outcomes + gene sets for 10_'s spot-check --------------------------------------------------
# Spot-check sets follow the project's 0.5 / 0.05 and 0.5 / 0.3 rules, applied to the TE level:
# efficiently translated (TE > 0.5, padj < 0.05) vs near-average (|TE| < 0.5, padj > 0.3). They
# test the join and the sign only; they are defined from the same values.
fm <- readRDS(here("output", "translational_control_overview", "feature_matrix_tco_hela.rds"))
for (k in c("ev", "bp")) {
  cond <- paste0("hela_teleman_", k)
  out <- te %>% transmute(gene_id, symbol, te_level = .data[[k]], te_level_padj = .data[[paste0("padj_", k)]])
  write_csv(out, here("output", paste0("teleman_absolute_te_", k, ".csv")))
  pos <- out %>% filter(te_level > 0.5, te_level_padj < 0.05)
  neg <- out %>% filter(abs(te_level) < 0.5, te_level_padj > 0.3)
  write_csv(pos %>% transmute(db_gene_symbol = symbol, ensembl_gene = gene_id, te_level,
                              gs_name = paste0(cond, "_high_te"),
                              gs_description = "TE level > 0.5 and padj < 0.05; spot-check set, not a factor-dependent set"),
            here("output", "genesets", paste0(cond, "_4e_promotes_TE_level_lfc0.5.csv")))
  write_csv(neg %>% dplyr::select(gene_id, symbol, te_level),
            here("output", "genesets", paste0(cond, "_4e_negative_controls_lfc0.5.csv")))
  cat(sprintf("%s: %d genes (%d in the HeLa matrix) | spot-check sets: %d high-TE, %d near-average | median TE %.3f\n",
              k, nrow(out), sum(out$gene_id %in% fm$gene_id_clean), nrow(pos), nrow(neg), median(out$te_level)))
}
