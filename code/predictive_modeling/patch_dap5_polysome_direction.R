# patch_dap5_polysome_direction.R
#
# Patches dap5_promoted_score_polysome_lfc (and dap5_polysome_lfc) in all
# feature matrix .rds files directly, without regenerating 01c.
#
# Run with: Rscript code/predictive_modeling/patch_dap5_polysome_direction.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ── Load updated DESeq2 results (NS control vs DAP5 shKD, positive = promoted) ─
dap5_new <- read_csv(
  here("output", "dap5_deseq2_full_results.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(symbol)) %>%
  select(symbol, new_lfc = log2FoldChange) %>%
  distinct(symbol, .keep_all = TRUE)

cat("Loaded updated DAP5 results:", nrow(dap5_new), "genes\n")
cat("Median LFC (expect slightly negative):", round(median(dap5_new$new_lfc), 3), "\n\n")

# ── Files and column names to patch ───────────────────────────────────────────
rds_dir <- here("output", "predictive_modeling")

targets <- list(
  list(
    file   = "feature_matrix_1hr_combined.rds",
    col    = "dap5_promoted_score_polysome_lfc"
  ),
  list(
    file   = "feature_matrix_1hr_clip_split.rds",
    col    = "dap5_promoted_score_polysome_lfc"
  ),
  list(
    file   = "feature_matrix_te_lfc.rds",
    col    = "dap5_promoted_score_polysome_lfc"
  ),
  list(
    file   = "feature_matrix_1hr_targeted.rds",
    col    = "dap5_polysome_lfc"
  ),
  list(
    file   = "feature_matrix_1hr_structure.rds",
    col    = "dap5_promoted_score_polysome_lfc"
  ),
  list(
    file   = "feature_matrix_codon_optimality.rds",
    col    = "dap5_promoted_score_polysome_lfc"
  )
)

# ── Patch each file ────────────────────────────────────────────────────────────
for (t in targets) {
  path <- file.path(rds_dir, t$file)

  if (!file.exists(path)) {
    cat("SKIP (not found):", t$file, "\n")
    next
  }

  fm <- readRDS(path)

  if (!t$col %in% names(fm)) {
    cat("SKIP (column absent):", t$file, "/", t$col, "\n")
    next
  }

  old_median <- round(median(fm[[t$col]], na.rm = TRUE), 3)

  # Join new values by symbol
  fm <- fm %>%
    left_join(dap5_new, by = "symbol") %>%
    mutate(!!t$col := ifelse(!is.na(new_lfc), new_lfc, .data[[t$col]])) %>%
    select(-new_lfc)

  new_median <- round(median(fm[[t$col]], na.rm = TRUE), 3)

  saveRDS(fm, path)
  cat("PATCHED:", t$file, "\n")
  cat("  Column:", t$col, "\n")
  cat("  Median before:", old_median, "-> after:", new_median, "\n\n")
}

cat("Done. Verify: TP53BP1 in feature_matrix_te_lfc.rds should be ~+1.83\n")
fm_check <- readRDS(file.path(rds_dir, "feature_matrix_te_lfc.rds"))
col_check <- if ("dap5_promoted_score_polysome_lfc" %in% names(fm_check))
               "dap5_promoted_score_polysome_lfc" else "dap5_polysome_lfc"
check_row <- fm_check[!is.na(fm_check$symbol) & fm_check$symbol == "TP53BP1",
                      c("symbol", col_check)]
print(check_row)
