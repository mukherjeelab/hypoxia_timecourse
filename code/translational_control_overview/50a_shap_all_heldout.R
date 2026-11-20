# 50a_shap_all_heldout.R
# TreeSHAP on EVERY held-out gene of a saved regression model, instead of 11_shap_regression.Rmd's
# random 200 (explain_n). Feeds the all-held-out version of 50_raw_shap_heatmap_hypvsnor.Rmd
# (requested by Kate, 2026-11-19).
#
# A gene's SHAP values depend only on the trained forest (unified on its training data) and that
# gene's own features - not on which other genes are explained. So the 200 genes 11_ explained
# must come out bit-identical here; that is asserted, and it is the check that this script
# explains the same model the same way. Only the AVERAGES (mean |SHAP|, share, direction) change
# with the number of genes.
#
#   Rscript code/translational_control_overview/50a_shap_all_heldout.R <suffix> [<suffix> ...]
# Writes shapreg_per_gene_allheldout_{suffix}.csv and shapreg_feature_allheldout_{suffix}.csv.
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here); library(ranger); library(treeshap) })
O <- here("output", "translational_control_overview")
suffixes <- commandArgs(TRUE)
stopifnot("give at least one model suffix" = length(suffixes) > 0)
ann <- read_csv(file.path(O, "feature_annotation.csv"), show_col_types = FALSE)

for (s in suffixes) {
  t0 <- Sys.time()
  rf <- readRDS(file.path(O, paste0("rfreg_model_", s, ".rds")))
  md <- readRDS(file.path(O, paste0("rfreg_model_data_", s, ".rds")))
  stopifnot("saved data does not belong to this model" = identical(md$suffix, s),
            "not a regression forest" = rf$treetype == "Regression")
  fc <- md$feature_cols; te <- md$test
  # Same construction as 11_'s shap_setup chunk, so the values are the same quantity.
  unified <- ranger.unify(rf, as.data.frame(md$train[, fc]))
  X <- as.data.frame(te[, fc])
  S <- as.data.frame(treeshap(unified, X, interactions = FALSE, verbose = FALSE)$shaps[, fc])
  pred <- predict(rf, data = X)$predictions
  implied <- pred - rowSums(S)
  stopifnot("SHAP does not close on the predictions" = sd(implied) < 1e-8)

  # Reproduce 11_'s 200 genes exactly: per-gene SHAP is independent of the explain set.
  old <- read_csv(file.path(O, paste0("shapreg_per_gene_", s, ".csv")), show_col_types = FALSE)
  m <- match(old$gene_id_clean, te$gene_id_clean)
  stopifnot("an 11_ explain gene is not in the held-out set" = !anyNA(m))
  dev <- max(abs(as.matrix(S[m, fc]) - as.matrix(old[, fc])))
  stopifnot("SHAP for 11_'s 200 genes is not reproduced - not the same model or method" = dev < 1e-10)

  # Direction exactly as 11_ defines it (Spearman of value vs own SHAP; NA if either is constant).
  dir <- vapply(fc, function(f) {
    v <- X[[f]]; sv <- S[[f]]
    if (sd(v, na.rm = TRUE) == 0 || sd(sv, na.rm = TRUE) == 0) return(NA_real_)
    cor(v, sv, use = "complete.obs", method = "spearman") }, numeric(1))
  ft <- tibble(feature = fc, mean_abs_shap = colMeans(abs(S)), shap_value_cor = dir) %>%
    left_join(ann %>% dplyr::select(feature, family, block), by = "feature") %>%
    mutate(share_pct = 100 * mean_abs_shap / sum(mean_abs_shap)) %>% arrange(desc(mean_abs_shap))

  write_csv(bind_cols(tibble(gene_id_clean = te$gene_id_clean, symbol = te$symbol, y = te$y, pred = pred), S),
            file.path(O, paste0("shapreg_per_gene_allheldout_", s, ".csv")))
  write_csv(ft, file.path(O, paste0("shapreg_feature_allheldout_", s, ".csv")))
  cat(sprintf("%s\n  %d held-out genes | 200-gene reproduction max dev %.1e | baseline %.5f | %.1f min\n",
              s, nrow(X), dev, mean(implied), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
