# 03k_shap_bar_gene_bootstrap.R - the beeswarm's dots summarised as bars with a bootstrap CI.
#   bar      = mean |SHAP| over the explained genes x sign(shap_value_cor) - 03g's signed
#              convention, in raw SHAP units so it sits on the same axis as 03j's beeswarm
#   interval = 95% percentile bootstrap over GENES: resample the explained genes with
#              replacement and recompute the bar (sign included), B times
# What the interval does and does not cover: it is how precisely 200 genes pin down the summary
# for ONE fitted model. It does not include model refit variability (the negative draw), which
# 03i's 20-draw whiskers measure - so these intervals are narrower than a refit-based one and
# must not be read as the uncertainty of the feature's importance in general.
# A whisker crossing zero means the direction is not settled even within this one model.
#   Rscript code/translational_control_overview/03k_shap_bar_gene_bootstrap.R [suffix] [top_n] [engine]
# Interactively: edit the settings block below and run top to bottom.
suppressPackageStartupMessages({ library(tidyverse); library(here) })

## ---- settings: edit these when running interactively (Positron / RStudio) ----
# Command-line arguments, when given, override them. Defaults match 03j's beeswarm.
suffix <- "eif33d_promotes_hypoxia_1hr_lfc0.5_tco4e_corrected_both_clipaggregate"
top_n  <- 20L
engine <- "kernelshap"   # or "treeshap"
n_boot <- 2000L
## ---------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) suffix <- args[1]
if (length(args) >= 2) top_n  <- as.integer(args[2])
if (length(args) >= 3) engine <- args[3]
stopifnot("engine must be treeshap or kernelshap" = engine %in% c("treeshap", "kernelshap"))
# 03_shap.Rmd tags its SHAP outputs with the engine; the MODEL files are shared between engines
etag <- if (engine == "treeshap") "" else paste0("_", engine)
O <- here("output", "translational_control_overview")

f <- c(model = file.path(O, paste0("rf_classifier_", suffix, ".rds")),
       data  = file.path(O, paste0("rf_model_data_", suffix, ".rds")),
       gene  = file.path(O, paste0("shap_per_gene_", suffix, etag, ".csv")),
       feat  = file.path(O, paste0("shap_feature_", suffix, etag, ".csv")))
stopifnot("model or SHAP output missing" = all(file.exists(f)),
          "SHAP output predates the model it claims to explain - re-render 03_shap.Rmd" =
            file.mtime(f["gene"]) > file.mtime(f["model"]))

md   <- readRDS(f["data"])
shap <- read_csv(f["gene"], show_col_types = FALSE)
ann  <- read_csv(file.path(O, "feature_annotation.csv"), show_col_types = FALSE)
ref  <- read_csv(f["feat"], show_col_types = FALSE)
stopifnot(identical(md$suffix, suffix),
          setequal(setdiff(names(shap), c("gene_id_clean", "symbol")), md$feature_cols),
          !anyDuplicated(md$test$gene_id_clean), all(shap$gene_id_clean %in% md$test$gene_id_clean))
# feature values = the exact imputed inputs 03_ explained, matched by gene (as in 03j)
vals <- md$test[match(shap$gene_id_clean, md$test$gene_id_clean), md$feature_cols] %>%
  mutate(across(everything(), as.numeric))
S <- as.matrix(shap[, md$feature_cols])
V <- as.matrix(vals)
n_genes <- nrow(S)

# signed summary for a set of gene rows: mean |SHAP| x sign of Spearman(value, SHAP).
# A feature constant over the resample has no direction and gets sign 0.
signed_summary <- function(rows, feats) {
  vapply(feats, function(j) {
    s <- S[rows, j]; v <- V[rows, j]
    r <- if (sd(s) == 0 || sd(v) == 0) 0 else cor(v, s, method = "spearman")
    mean(abs(s)) * sign(r)
  }, numeric(1))
}

point <- tibble(feature = md$feature_cols,
                mean_abs_shap = unname(colMeans(abs(S))),
                signed = unname(signed_summary(seq_len(n_genes), md$feature_cols)))
# must reproduce 03_'s feature table, or this is summarising something else
stopifnot("mean |SHAP| does not reproduce 03_'s feature table" =
            isTRUE(all.equal(point$mean_abs_shap, ref$mean_abs_shap[match(point$feature, ref$feature)])),
          "direction does not reproduce 03_'s shap_value_cor sign" =
            all(sign(point$signed) == sign(ref$shap_value_cor[match(point$feature, ref$feature)]), na.rm = TRUE))

top <- point %>% slice_max(mean_abs_shap, n = top_n, with_ties = FALSE)

set.seed(9)
boot <- replicate(n_boot, signed_summary(sample(n_genes, replace = TRUE), top$feature))
ci <- tibble(feature = top$feature,
             lo = apply(boot, 1, quantile, 0.025), hi = apply(boot, 1, quantile, 0.975),
             frac_same_sign = rowMeans(sign(boot) == sign(top$signed)))

category_colors <- c("Intrinsic mRNA features" = "#3D9C60", "eIF3 binding evidence" = "#E6A817",
                     "Sensitivity to DAP5" = "#D668C6", "Sensitivity to eIF4E" = "#4A7FC1")
plot_df <- top %>% left_join(ci, by = "feature") %>%
  left_join(ann %>% dplyr::select(feature, block, family), by = "feature") %>%
  mutate(category = case_when(family == "clip" ~ "eIF3 binding evidence",
                              family == "dap5" ~ "Sensitivity to DAP5",
                              block  == "eif4e" ~ "Sensitivity to eIF4E",
                              TRUE ~ "Intrinsic mRNA features"),
         category = factor(category, levels = names(category_colors)),
         feature  = fct_reorder(feature, signed))
stopifnot("feature missing from annotation" = !anyNA(plot_df$block))
crosses <- plot_df$feature[plot_df$lo < 0 & plot_df$hi > 0]
cat(sprintf("top %d of %d features | %d genes | %d bootstrap resamples\n",
            nrow(plot_df), length(md$feature_cols), n_genes, n_boot))
cat("95% CI crosses zero:", if (length(crosses)) paste(crosses, collapse = ", ") else "none", "\n")

what <- sub("_lfc.*", "", sub("^eif33d_promotes_", "si3d ", suffix))
p <- ggplot(plot_df, aes(signed, feature, fill = category)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.35, linewidth = 0.4, colour = "grey20") +
  scale_fill_manual(values = category_colors, name = NULL, drop = TRUE) +
  labs(x = "mean |SHAP| signed by direction\n(right = higher values push toward eIF3d-promoted)",
       y = NULL,
       title = sprintf("Top %d features: SHAP importance with 95%% CI", nrow(plot_df)),
       subtitle = sprintf("Classifier, %s | %s | n = %d held-out genes explained, n = %d features",
                          if (engine == "kernelshap") "Kernel SHAP on the probability forest" else "TreeSHAP",
                          what, n_genes, length(md$feature_cols)),
       caption = sprintf(paste0("Bar = mean |SHAP| x sign of Spearman(feature value, SHAP). Whiskers = 95%% bootstrap over ",
                                "the %d explained genes (%d resamples, seed 9).\n",
                                "Covers which genes were explained, NOT model refit variability (see 03i)."),
                         n_genes, n_boot)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        plot.caption.position = "plot", plot.title.position = "plot")
print(p)   # shows in the Plots pane when run interactively

out_pdf <- here("plots", paste0("03k_shap_bar_gene_bootstrap_", suffix, etag, "_top", top_n, ".pdf"))
ggsave(out_pdf, p, width = 9, height = 0.3 * nrow(plot_df) + 2.8)
write_csv(plot_df %>% arrange(desc(mean_abs_shap)),
          file.path(O, paste0("shap_bar_gene_bootstrap_", suffix, etag, "_top", top_n, ".csv")))
cat("wrote", out_pdf, "\n")
