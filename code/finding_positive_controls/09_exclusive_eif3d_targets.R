# eIF3d-promoted mRNAs in hypoxia whose effect is purely or nearly purely translational in BOTH
# lines, annotated with 5'UTR features - candidates for a 5'UTR reporter (Kate, 2026-11-20:
# "focus more on exclusive or near exclusive targets", 5'UTR features noted, not filtered on).
#
# Per line, si3d hypoxia 1hr, log2(sictrl / si3d), positive = eIF3d promotes:
#   eIF3d target : te_lfc > MIN_TE and te_padj < 0.05 (the project's target rule)
#   Exclusive    : category from categorize_translation_changes() (ribo + TE sig, RNA not)
#   near-Exclusive: not Exclusive, but footprints fall significantly (ribo_lfc > 0, padj < 0.05) and the
#                  mRNA change is small: |rna_lfc - median rna_lfc of that contrast| < RNA_TOL.
#                  Centred on the median because 231 is spike-in normalized (median ~ +0.4), so
#                  a padj-only RNA test calls the typical gene "changed".
# Kept: target in both lines AND (Exclusive or near-Exclusive) in both.
# Ranked by the smaller of the two TE LFCs, so both lines must carry the effect.
#
# MCF7-SIX1 ribo/RNA/category: 00_mcf7six1_hypoxia_si3d_categories.Rmd (kd = si3d); its RNA LFC
# here is the apeglm-shrunk one, matching 231. 5'UTR features: each line's own feature matrix
# (its expressed transcript); MCF7-SIX1 transcript picks can be off where the v26 quant's top
# isoform was retired (see POSITIVE_CONTROLS.md, ACTB) - check isoforms before cloning.
#
# Usage: Rscript code/finding_positive_controls/09_exclusive_eif3d_targets.R [min_te] [rna_tol]
#   defaults 0.5 and 0.3
# Writes output/finding_positive_controls/09_exclusive_eif3d_targets_te{min_te}_rna{rna_tol}.csv
suppressPackageStartupMessages({ library(tidyverse); library(here) })
args <- commandArgs(trailingOnly = TRUE)
MIN_TE  <- if (length(args) >= 1) as.numeric(args[1]) else 0.5
RNA_TOL <- if (length(args) >= 2) as.numeric(args[2]) else 0.3
P <- 0.05
out_dir <- here("output", "finding_positive_controls")
clean <- function(x) sub("\\..*", "", x)

h231 <- read_csv(here("output", "translation_categories_si3d_vs_sictrl_hypoxia_1hr.csv"), show_col_types = FALSE) %>%
  mutate(g = clean(gene_id))
n231 <- read_csv(here("output", "translation_categories_si3d_vs_sictrl_normoxia_1and4hr.csv"), show_col_types = FALSE) %>%
  mutate(g = clean(gene_id))
hmcf <- read_csv(here("output", "translation_categories_si3d_vs_sictrl_mcf7six1_hypoxia_1hr.csv"), show_col_types = FALSE) %>%
  mutate(g = clean(gene_id))
nmcf <- read_csv(here("output", "translation_categories_si3d_vs_sictrl_mcf7six1_normoxia_1hr.csv"), show_col_types = FALSE) %>%
  mutate(g = clean(gene_id))
cmcf <- read_csv(file.path(out_dir, "mcf7six1_hypoxia_1hr_si3d_rna_ribo_te_categories.csv"), show_col_types = FALSE) %>%
  mutate(g = clean(gene_id))
stopifnot("231 direction" = h231$rna_lfc[h231$symbol %in% "EIF3D"] > 1.5,
          "MCF7-SIX1 negation" = isTRUE(all.equal(hmcf$te_lfc, -hmcf$source_log2FoldChange)),
          "MCF7-SIX1 RNA direction" = cmcf$rna_lfc_shrunk[cmcf$symbol %in% "EIF3D"] > 1.5,
          !anyDuplicated(h231$g), !anyDuplicated(hmcf$g), !anyDuplicated(cmcf$g))

med231 <- median(h231$rna_lfc, na.rm = TRUE); medmcf <- median(cmcf$rna_lfc_shrunk, na.rm = TRUE)
cat(sprintf("median rna_lfc: 231 %.3f, MCF7-SIX1 %.3f\n", med231, medmcf))

d <- h231 %>%
  transmute(g, symbol, te_231 = te_lfc, tep_231 = te_padj, ribo_231 = ribo_lfc, ribop_231 = ribo_padj,
            rna_231 = rna_lfc, rnap_231 = rna_padj, cat_231 = category) %>%
  inner_join(hmcf %>% transmute(g, te_mcf = te_lfc, tep_mcf = te_padj), by = "g") %>%
  inner_join(cmcf %>% transmute(g, ribo_mcf = ribo_lfc, ribop_mcf = ribo_padj, rna_mcf = rna_lfc_shrunk,
                                rnap_mcf = rna_padj, cat_mcf = category), by = "g") %>%
  left_join(n231 %>% transmute(g, te_nor_231 = te_lfc, tep_nor_231 = te_padj), by = "g") %>%
  left_join(nmcf %>% transmute(g, te_nor_mcf = te_lfc, tep_nor_mcf = te_padj), by = "g")

mode_call <- function(cat, ribo, ribop, rna, med) case_when(
  cat == "Exclusive" ~ "Exclusive",
  ribo > 0 & ribop < P & abs(rna - med) < RNA_TOL ~ "near-Exclusive",
  TRUE ~ "other")
d <- d %>% mutate(
  target_231 = te_231 > MIN_TE & tep_231 < P, target_mcf = te_mcf > MIN_TE & tep_mcf < P,
  mode_231 = mode_call(cat_231, ribo_231, ribop_231, rna_231, med231),
  mode_mcf = mode_call(cat_mcf, ribo_mcf, ribop_mcf, rna_mcf, medmcf))

both_t <- d %>% filter(target_231 %in% TRUE, target_mcf %in% TRUE)
cat("eIF3d targets in both lines (TE >", MIN_TE, ", padj < 0.05):", nrow(both_t), "\n")
print(table(`231` = both_t$mode_231, `MCF7-SIX1` = both_t$mode_mcf))

out <- both_t %>% filter(mode_231 != "other", mode_mcf != "other") %>%
  mutate(min_te = pmin(te_231, te_mcf),
         both_exclusive = mode_231 == "Exclusive" & mode_mcf == "Exclusive",
         # hypoxia / normoxia TE per line; NA where normoxia <= 0 (no positive normoxia effect)
         hyp_nor_ratio_231 = if_else(te_nor_231 > 0, te_231 / te_nor_231, NA_real_),
         hyp_nor_ratio_mcf = if_else(te_nor_mcf > 0, te_mcf / te_nor_mcf, NA_real_))
cat("... Exclusive or near-Exclusive in both:", nrow(out), "(both strictly Exclusive:", sum(out$both_exclusive), ")\n")

# 5'UTR features, each line on its own expressed transcript
feat <- c("transcript_id_clean", "utr5_length", "utr5_gc", "kozak_pwm_score_v2",
          "tco_struct_accessibility_cap_proximal", "tco_struct_accessibility_utr5_mean",
          "g4mer_mean_utr5", "cds_length", "cds_gc")
fm <- function(f, sfx) readRDS(here("output", "translational_control_overview", f)) %>%
  transmute(g = gene_id_clean, across(all_of(feat))) %>% distinct(g, .keep_all = TRUE) %>%
  rename_with(~ paste0(.x, "_", sfx), -g)
f231 <- fm("feature_matrix_tco.rds", "231"); fmcf <- fm("feature_matrix_tco_mcf7six1.rds", "mcf")
# percentile of each value among all genes in that line's matrix, so "long" / "GC-rich" is relative
pct <- function(x, ref) round(100 * ecdf(ref[!is.na(ref)])(x))
out <- out %>% left_join(f231, by = "g") %>% left_join(fmcf, by = "g") %>%
  mutate(utr5_length_pctile_231 = pct(utr5_length_231, f231$utr5_length_231),
         utr5_gc_pctile_231 = pct(utr5_gc_231, f231$utr5_gc_231),
         utr5_length_pctile_mcf = pct(utr5_length_mcf, fmcf$utr5_length_mcf),
         utr5_gc_pctile_mcf = pct(utr5_gc_mcf, fmcf$utr5_gc_mcf),
         same_tx = transcript_id_clean_231 == transcript_id_clean_mcf) %>%
  arrange(desc(min_te)) %>% mutate(rank = row_number()) %>%
  relocate(rank, symbol, both_exclusive, mode_231, mode_mcf, min_te, te_231, te_mcf)
stopifnot(all(out$te_231 > MIN_TE), all(out$te_mcf > MIN_TE), all(out$tep_231 < P), all(out$tep_mcf < P))

write_csv(out %>% dplyr::rename(gene_id_clean = g),
          file.path(out_dir, sprintf("09_exclusive_eif3d_targets_te%s_rna%s.csv", MIN_TE, RNA_TOL)))
print(as.data.frame(out %>% transmute(rank, symbol, m231 = substr(mode_231, 1, 4), mmcf = substr(mode_mcf, 1, 4),
  te_231 = round(te_231, 2), te_mcf = round(te_mcf, 2), rna_231 = round(rna_231 - med231, 2), rna_mcf = round(rna_mcf - medmcf, 2),
  nor_231 = round(te_nor_231, 2), nor_mcf = round(te_nor_mcf, 2),
  u5len = paste0(utr5_length_231, "/", utr5_length_mcf), u5gc = paste0(round(utr5_gc_231), "/", round(utr5_gc_mcf)),
  len_pct = utr5_length_pctile_231, same_tx)), row.names = FALSE)
