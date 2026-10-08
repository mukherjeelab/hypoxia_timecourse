# 03g_signed_shap_lollipop.R - top-N lollipop of SHAP importance signed by direction:
#   signed share = share of total mean |SHAP| (%) x sign(shap_value_cor)
# the convention of 07c_rf_shap_eif3d_targets.Rmd ("Directional Importance") and of
# 15_shap_heatmap.Rmd's signed share. shap_value_cor is the correlation of a feature's VALUE with
# its own SHAP value over the explained genes, so a bar points right when high values push the
# prediction up. Share rather than raw mean |SHAP| so two models (e.g. hypoxia vs normoxia) sit on
# one scale. Colours and categories follow notebook 07's signed-importance lollipop.
# (Kate, 2026-10-06. No 20-draw sign filter: only the top N are shown.)
#
# Works on either model type, chosen by the suffix:
#   classifier (02_/03_, "_tco_"):     shap_feature_{suffix}.csv,    Gini importance
#   regression (10_/11_, "_tcoreg_"):  shapreg_feature_{suffix}.csv, variance-reduction importance
# Both explain a random sample of held-out genes with TreeSHAP (03_ via its default engine,
# per 03c_kernelshap_comparison; 11_ directly, as a regression forest needs no surrogate).
#   Rscript code/translational_control_overview/03g_signed_shap_lollipop.R [suffix] [top_n]
suppressPackageStartupMessages({ library(tidyverse); library(here) })

args   <- commandArgs(trailingOnly = TRUE)
suffix <- if (length(args) >= 1) args[1] else
  "eif33d_promotes_hypoxia_1hr_lfc0.5_tco_corrected_both_clipaggregate"
top_n  <- if (length(args) >= 2) as.integer(args[2]) else 30L
# "kernelshap" plots the Shapley values of the ACTUAL probability forest instead of the
# surrogate's. 03_shap.Rmd suffixes its outputs with the engine but the MODEL files are shared,
# so the engine tag is appended to the SHAP and output paths only.
engine <- if (length(args) >= 3) args[3] else "treeshap"
stopifnot("engine must be treeshap or kernelshap" = engine %in% c("treeshap", "kernelshap"))
etag <- if (engine == "treeshap") "" else paste0("_", engine)
O <- here("output", "translational_control_overview")

is_reg <- grepl("_tcoreg_", suffix)
pre <- if (is_reg) {
  c(model = "rfreg_model_", data = "rfreg_model_data_", imp = "rfreg_importance_",
    shap = "shapreg_feature_", gene = "shapreg_per_gene_")
} else {
  c(model = "rf_classifier_", data = "rf_model_data_", imp = "rf_importance_",
    shap = "shap_feature_", gene = "shap_per_gene_")
}
ext <- c(model = ".rds", data = ".rds", imp = ".csv", shap = ".csv", gene = ".csv")
# the engine tag belongs only to this notebook's SHAP outputs, never to the shared model files
tag <- c(model = "", data = "", imp = "", shap = etag, gene = etag)
f <- setNames(file.path(O, paste0(pre, suffix, tag[names(pre)], ext[names(pre)])), names(pre))
stopifnot("model or SHAP output missing - fit the model and render its SHAP notebook" = all(file.exists(f)))
# Staleness guard: SHAP written before the model was last refit explains a different model.
stopifnot("SHAP output predates the model it claims to explain - re-render its SHAP notebook" =
            file.mtime(f["shap"]) > file.mtime(f["model"]))

md  <- readRDS(f["data"])
stopifnot(identical(md$suffix, suffix))
ann <- read_csv(file.path(O, "feature_annotation.csv"), show_col_types = FALSE)
imp <- read_csv(f["imp"], show_col_types = FALSE)
n_expl <- nrow(read_csv(f["gene"], show_col_types = FALSE))

# Notebook 07's categories, applied by FAMILY. 07 matched name prefixes (^clip_, ^dap5_, ...), which
# would file is_dap5_promoted / is_dap5_repressed under "intrinsic"; the family is the definition.
# The external stability datasets (half-lives, DHX29, CNOT3) fall through to intrinsic, as in 07.
category_colors <- c("Intrinsic mRNA features" = "#3D9C60", "eIF3 binding evidence" = "#E6A817",
                     "Sensitivity to DAP5" = "#D668C6", "Sensitivity to eIF4E" = "#4A7FC1")
tbl <- read_csv(f["shap"], show_col_types = FALSE) %>%
  dplyr::select(feature, mean_abs_shap, shap_value_cor) %>%
  left_join(ann %>% dplyr::select(feature, block, family), by = "feature") %>%
  left_join(imp %>% dplyr::select(feature, impurity = importance), by = "feature") %>%
  mutate(share_pct = 100 * mean_abs_shap / sum(mean_abs_shap),
         signed_share = share_pct * sign(shap_value_cor),
         category = case_when(family == "clip" ~ "eIF3 binding evidence",
                              family == "dap5" ~ "Sensitivity to DAP5",
                              block  == "eif4e" ~ "Sensitivity to eIF4E",
                              TRUE ~ "Intrinsic mRNA features"))
stopifnot("SHAP table does not cover exactly the model's features" = setequal(tbl$feature, md$feature_cols),
          "feature missing from annotation or impurity importance" = !anyNA(tbl$block), !anyNA(tbl$impurity))

# Impurity is unsigned, so it is compared with mean |SHAP|, never with the signed value.
rho <- cor(tbl$mean_abs_shap, tbl$impurity, method = "spearman")
imp_name <- if (is_reg) "impurity (variance)" else "Gini"
cat(sprintf("%s vs mean |SHAP|: Spearman %.3f (n = %d features) | %d genes explained\n",
            imp_name, rho, nrow(tbl), n_expl))

top <- tbl %>% slice_max(mean_abs_shap, n = top_n, with_ties = FALSE)
stopifnot("a top feature has no direction (constant value or SHAP) - cannot be signed" =
            !any(is.na(top$shap_value_cor) | top$shap_value_cor == 0))
cat(sprintf("top %d: %d point right, %d left | |direction| < 0.3: %s\n", nrow(top),
            sum(top$signed_share > 0), sum(top$signed_share < 0),
            paste(top$feature[abs(top$shap_value_cor) < 0.3], collapse = ", ")))

what <- sub("_lfc.*", "", sub("^eif33d_promotes_", "si3d ", suffix))
p <- top %>%
  mutate(feature = fct_reorder(feature, signed_share),
         category = factor(category, levels = names(category_colors))) %>%
  ggplot(aes(signed_share, feature, colour = category)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = signed_share, yend = feature), colour = "grey80", linewidth = 0.5) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = category_colors, name = NULL, drop = TRUE) +
  labs(x = if (is_reg) "signed SHAP share (%)\n(right = higher values predict a larger TE drop on si3d)"
           else        "signed SHAP share (%)\n(right = higher values push toward eIF3d-promoted)",
       y = NULL,
       title = sprintf("Top %d features: SHAP importance signed by direction", nrow(top)),
       subtitle = sprintf("%s | %s | n = %d held-out genes explained\n%s vs mean |SHAP| Spearman = %.2f, n = %d features",
                          paste0(if (is_reg) "Regression (TE LFC, all genes)" else "Classifier (promoted vs negative controls)",
                                 if (engine == "kernelshap") ", Kernel SHAP on the probability forest" else ""),
                          what, n_expl, imp_name, rho, nrow(tbl)),
       caption = "Bar = share of total mean |SHAP| x sign of cor(feature value, its SHAP value) over the explained genes.") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        plot.caption.position = "plot", plot.title.position = "plot")

out_pdf <- here("plots", paste0("03g_signed_shap_lollipop_", suffix, etag, ".pdf"))
ggsave(out_pdf, p, width = 8.5, height = 0.25 * nrow(top) + 2.8)
write_csv(tbl %>% arrange(desc(mean_abs_shap)), file.path(O, paste0("shap_signed_importance_", suffix, etag, ".csv")))
cat("wrote", out_pdf, "\n")
