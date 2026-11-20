# One row per reporter-panel gene: our TE / ribo / RNA changes and categories, eIF3d / eIF3e
# knockdown in HeLa, eIF3 binding, and notes. Every number is read from a file, not typed.
#
# Usage: Rscript code/finding_positive_controls/08_reporter_panel_table.R
# Needs: 05_eif3_binding_lookup.R run on the same genes (writes the _<GENES> summary), and
#        00_mcf7six1_hypoxia_si3d_categories.Rmd knit for hypoxia with kd = "si3d" and kd = "si3e"
#        (MCF7-SIX1 RNA/ribo/category).
# Writes output/finding_positive_controls/08_reporter_panel_table.csv
#
# Direction: every LFC is log2(sictrl / si3d or si3e); positive = the factor promotes.
# MDA-MB-231 = spike-in normalized DESeq2 (apeglm); MCF7-SIX1 TE = saved apeglm TE, RNA/ribo =
# unshrunk (00_); HeLa = Herrmannova et al. eLife 2024 - published data, not ours (10h_herrmannova_outcome_prep.Rmd).
suppressPackageStartupMessages({ library(tidyverse); library(here) })
panel <- tribble(
  ~symbol,   ~role,
  "HERPUD1", "selected eIF3d target",
  "TBPL1",   "selected eIF3d target",
  "RPP14",   "selected eIF3d target",
  "AGPAT1",  "selected eIF3d target",
  "BCL2L12", "positive control (Roiuk et al)",
  "TPRA1",   "positive control (Roiuk et al)",
  "ACTB",    "negative control",
  "BSG",     "negative control")
P <- 0.05
out_dir <- here("output", "finding_positive_controls")
rd <- function(f) read_csv(here("output", paste0("translation_categories_", f, ".csv")),
                           show_col_types = FALSE) %>% mutate(g = sub("\\..*", "", gene_id))
fmt <- function(l, p) ifelse(is.na(l), NA, sprintf("%.2f (%s)", l, signif(p, 2)))

# --- our cell lines -------------------------------------------------------------------------
h231 <- rd("si3d_vs_sictrl_hypoxia_1hr"); n231 <- rd("si3d_vs_sictrl_normoxia_1and4hr")
hmcf <- rd("si3d_vs_sictrl_mcf7six1_hypoxia_1hr"); nmcf <- rd("si3d_vs_sictrl_mcf7six1_normoxia_1hr")
cmcf <- read_csv(file.path(out_dir, "mcf7six1_hypoxia_1hr_si3d_rna_ribo_te_categories.csv"),
                 show_col_types = FALSE) %>% mutate(g = sub("\\..*", "", gene_id))
stopifnot("231 direction" = h231$rna_lfc[h231$symbol %in% "EIF3D"] > 1.5,
          "MCF7-SIX1 negation" = isTRUE(all.equal(hmcf$te_lfc, -hmcf$source_log2FoldChange)),
          "MCF7-SIX1 RNA direction" = cmcf$rna_lfc[cmcf$symbol %in% "EIF3D"] > 1.5)

gid <- h231 %>% filter(symbol %in% panel$symbol) %>% dplyr::select(symbol, g)
stopifnot(setequal(gid$symbol, panel$symbol), !anyDuplicated(gid$symbol))

ours <- gid %>%
  left_join(h231 %>% transmute(g, te_231 = fmt(te_lfc, te_padj), ribo_231 = fmt(ribo_lfc, ribo_padj),
                               rna_231 = fmt(rna_lfc, rna_padj), category_231 = category), by = "g") %>%
  left_join(hmcf %>% transmute(g, te_mcf = fmt(te_lfc, te_padj)), by = "g") %>%
  left_join(cmcf %>% transmute(g, ribo_mcf = fmt(ribo_lfc, ribo_padj), rna_mcf = fmt(rna_lfc, rna_padj),
                               category_mcf = category), by = "g") %>%
  left_join(n231 %>% transmute(g, te_nor_231 = fmt(te_lfc, te_padj)), by = "g") %>%
  left_join(nmcf %>% transmute(g, te_nor_mcf = fmt(te_lfc, te_padj)), by = "g")

# --- eIF3e knockdown (sictrl vs si3e), hypoxia 1hr, both lines ---------------------------------
# 231: saved si3e file (TE/ribo/RNA/category). MCF7-SIX1: saved apeglm TE + 00_ knit with
# kd = "si3e" for ribo/RNA/category. eIF3e target = TE LFC > 0.5 and padj < 0.05, the same rule
# as the eIF3d target sets; TE LFC < -0.5 and padj < 0.05 is reported as eIF3e-inhibited.
e231 <- rd("si3e_vs_sictrl_hypoxia_1hr"); emcf <- rd("si3e_vs_sictrl_mcf7six1_hypoxia_1hr")
cemcf <- read_csv(file.path(out_dir, "mcf7six1_hypoxia_1hr_si3e_rna_ribo_te_categories.csv"),
                  show_col_types = FALSE) %>% mutate(g = sub("\\..*", "", gene_id))
stopifnot("231 si3e direction" = e231$rna_lfc[e231$symbol %in% "EIF3E"] > 1.5,
          "MCF7-SIX1 si3e negation" = isTRUE(all.equal(emcf$te_lfc, -emcf$source_log2FoldChange)),
          "MCF7-SIX1 si3e RNA direction" = cemcf$rna_lfc[cemcf$symbol %in% "EIF3E"] > 1.5)
# mRNA call, relative to the transcriptome median of that contrast: 231 is spike-in normalized and
# knockdown lowers total mRNA per cell, so its median rna_lfc sits near +0.4 and a call against 0
# would mark the typical gene "lower". Significant (padj < 0.05) AND >= 0.25 log2 from the median.
rna_call <- function(l, p, med) case_when(is.na(l) | is.na(p) ~ "not tested",
                                          p < P & l - med >= 0.25 ~ "lower",
                                          p < P & l - med <= -0.25 ~ "higher",
                                          TRUE ~ "unchanged")
call3e <- function(l, p) case_when(is.na(l) | is.na(p) ~ "not tested",
                                   l > 0.5 & p < P ~ "target", l < -0.5 & p < P ~ "inhibited",
                                   TRUE ~ "not a target")
si3e <- gid %>%
  left_join(e231 %>% transmute(g, si3e_te_231 = fmt(te_lfc, te_padj), si3e_ribo_231 = fmt(ribo_lfc, ribo_padj),
                               si3e_rna_231 = fmt(rna_lfc, rna_padj), si3e_category_231 = category,
                               c231 = call3e(te_lfc, te_padj),
                               r231 = rna_call(rna_lfc, rna_padj, median(e231$rna_lfc, na.rm = TRUE))), by = "g") %>%
  left_join(emcf %>% transmute(g, si3e_te_mcf = fmt(te_lfc, te_padj), cmcf = call3e(te_lfc, te_padj)), by = "g") %>%
  left_join(cemcf %>% transmute(g, si3e_ribo_mcf = fmt(ribo_lfc, ribo_padj), si3e_rna_mcf = fmt(rna_lfc, rna_padj),
                                si3e_category_mcf = category,
                                rmcf = rna_call(rna_lfc, rna_padj, median(cemcf$rna_lfc, na.rm = TRUE))), by = "g") %>%
  mutate(eif3e_target = case_when(
    c231 == "target" & cmcf == "target" ~ "target in both lines",
    c231 == "target" ~ paste0("target in MDA-MB-231 only (MCF7-SIX1: ", cmcf, ")"),
    cmcf == "target" ~ paste0("target in MCF7-SIX1 only (MDA-MB-231: ", c231, ")"),
    c231 == "inhibited" | cmcf == "inhibited" ~
      sprintf("not a target; inhibited (MDA-MB-231: %s, MCF7-SIX1: %s)", c231, cmcf),
    TRUE ~ "not a target")) %>%
  mutate(eif3e_rna_change = sprintf("mRNA after si3e: %s in MDA-MB-231, %s in MCF7-SIX1", r231, rmcf)) %>%
  dplyr::select(-c231, -cmcf, -r231, -rmcf)
ours <- ours %>% left_join(si3e %>% dplyr::select(-symbol), by = "g")

# --- published data (not ours): HeLa eIF3d / eIF3e knockdown, Herrmannova et al. eLife 2024 ----
hela3d <- rd("si3d_vs_sictrl_hela_normoxia_steadystate"); hela3e <- rd("si3e_vs_sictrl_hela_normoxia_steadystate")
one <- function(d, s) { x <- d$rna_lfc[d$symbol %in% s]; stopifnot(length(x) == 1); x }  # HeLa has NA symbols
stopifnot("HeLa si3d direction" = one(hela3d, "EIF3D") > 1.5, "HeLa si3e direction" = one(hela3e, "EIF3E") > 1.5)
other <- gid %>%
  left_join(hela3d %>% transmute(g, te_hela_si3d = fmt(te_lfc, te_padj)), by = "g") %>%
  left_join(hela3e %>% transmute(g, te_hela_si3e = fmt(te_lfc, te_padj)), by = "g")

# --- binding (05_) ----------------------------------------------------------------------------
bf <- file.path(out_dir, paste0("05_eif3_binding_lookup_summary_", paste(panel$symbol, collapse = "_"), ".csv"))
stopifnot("run 05_eif3_binding_lookup.R on the panel genes first" = file.exists(bf))
n_runs <- function(x) str_count(coalesce(x, ""), "Replicate")
bind <- read_csv(bf, show_col_types = FALSE) %>%
  transmute(symbol,
            cap_binding_subunitseq = case_when(
              !is.na(cap_lamper2020) ~ paste("Lamper 2020:", cap_lamper2020),
              !is.na(cap_thapsigargin) ~ paste("Mukhopadhyay 2023 (thapsigargin):", cap_thapsigargin),
              TRUE ~ "not called"),
            lee2015_parclip = if_else(is.na(lee_parclip), "not called", lee_parclip),
            desilva_eif3d_band = sprintf("%d of 4 runs", n_runs(desilva_eif3d_l)),
            desilva_other_bands = sprintf("%d of 16 runs", n_runs(desilva_other)),
            cate_irclip = coalesce(cate_irclip, "not called"))

# --- notes (each from a named analysis in this folder) ----------------------------------------
notes <- tribble(
  ~symbol, ~notes,
  "HERPUD1", "Largest hypoxia TE effect in 231, but nearly as large in normoxia (on the diagonal): not hypoxia-specific (06_/07_). Also strongly eIF3d- and eIF3e-promoted in published HeLa data (Herrmannova 2024; TE 1.98 / 1.05), so the regulation replicates in an independent lab and cell line. No eIF3 binding called. Mixed isoforms: main transcript 39% (231) / 52% (MCF7-SIX1).",
  "TBPL1",   "Strong in both lines but equally strong in normoxia (on the diagonal). Also eIF3d- and eIF3e-promoted in published HeLa data (Herrmannova 2024; 1.74 / 1.53). Mixed isoforms: main transcript 35% / 45%. 3'UTR only 41 nt.",
  "RPP14",   "Most hypoxia-dominant strong target: hypoxia >= 2x normoxia in both lines (07_); 231 normoxia still significant. Only panel gene cap-bound by eIF3d (thapsigargin subunit-seq). One dominant isoform (84% / 61%). Not significant in published HeLa data (0.54, padj 0.42). Buffered in MCF7-SIX1: mRNA rises on knockdown. Low expression (~3 TPM).",
  "AGPAT1",  "Hypoxia-specific in 231 (normoxia ~0), but MCF7-SIX1 hypoxia TE 0.73 is below the strong-target floor of 1. Very low expression in 231 (~0.9 TPM). Main transcript IDs differ (231 ENST00000929388 / MCF7-SIX1 ENST00000375107) but are the same mRNA: 44 nt longer 5' end in the v49 model. Borderline in published HeLa data (0.71, padj 0.053).",
  "BCL2L12", "Exclusive in both lines (purely translational). The Teleman-tested transcript ENST00000441864 (510-nt 5'UTR, 5 upstream ATGs) is only ~1-3% of BCL2L12 mRNA; our cells express ENST00000246784 (143-nt 5'UTR, the last 143 nt of the 510).",
  "TPRA1",   "eIF3d target in MCF7-SIX1 only (Buffered); just misses in 231. Teleman-tested transcript ENST00000296210 (315-nt 5'UTR) is ~30% of TPRA1 in MCF7-SIX1, ~10% in 231, where most TPRA1 uses a different first exon (82-nt 5'UTR). Low expression (~3 TPM).",
  "ACTB",    "Selected negative control (Roiuk 2024). TE flat in every si3d contrast except 231 hypoxia 24hr (TE -1.14). Bound in CLIP, as expected for a very abundant mRNA. HBB (Roiuk's reporter negative control) has 0 RNA-seq reads in both lines.",
  "BSG",     "TE flat in all six si3d contrasts in our lines, well expressed (~125-150 TPM, one dominant isoform). BUT eIF3d inhibits it in published HeLa data (Herrmannova 2024; TE -1.11), and eIF3e too (HeLa -0.86; MCF7-SIX1 normoxia si3e -0.30, padj 0.026).")
stopifnot(setequal(notes$symbol, panel$symbol))

tab <- panel %>% left_join(ours, by = "symbol") %>% left_join(other, by = c("symbol", "g")) %>%
  left_join(bind, by = "symbol") %>% left_join(notes, by = "symbol") %>%
  dplyr::rename(gene_id_clean = g)
stopifnot(nrow(tab) == nrow(panel), !anyNA(tab$te_231), !anyNA(tab$te_mcf), !anyNA(tab$category_mcf),
          !anyNA(tab$si3e_te_231), !anyNA(tab$si3e_te_mcf), !anyNA(tab$si3e_category_mcf))
write_csv(tab, file.path(out_dir, "08_reporter_panel_table.csv"))
print(as.data.frame(tab %>% dplyr::select(-notes, -gene_id_clean)), right = FALSE)
