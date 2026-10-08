# 03j_shap_beeswarm.R - standard SHAP summary (beeswarm) plot for the top-N features.
# One point per explained gene per feature: x = that gene's SHAP value, colour = the gene's value
# of the feature, as a within-feature percentile so every row uses the full colour range
# (the shap package's convention). Direction is read directly: high-value (vermillion) points on
# the right mean high values push toward eIF3d-promoted. Rows ordered by mean |SHAP|, so the
# ranking matches 03g / 03i.
# Feature values come from the model's saved, imputed test split (02_'s rf_model_data), i.e. the
# exact inputs 03_ explained, matched by gene.
#   Rscript code/translational_control_overview/03j_shap_beeswarm.R [suffix] [top_n]
suppressPackageStartupMessages({ library(tidyverse); library(here); library(ggbeeswarm) })

args   <- commandArgs(trailingOnly = TRUE)
suffix <- if (length(args) >= 1) args[1] else
  "eif33d_promotes_hypoxia_1hr_lfc0.5_tco_corrected_both_clipaggregate"
top_n  <- if (length(args) >= 2) as.integer(args[2]) else 30L
O <- here("output", "translational_control_overview")

f <- c(model = file.path(O, paste0("rf_classifier_", suffix, ".rds")),
       data  = file.path(O, paste0("rf_model_data_", suffix, ".rds")),
       gene  = file.path(O, paste0("shap_per_gene_", suffix, ".csv")))
stopifnot("model or SHAP output missing" = all(file.exists(f)),
          "SHAP output predates the model it claims to explain - re-render 03_shap.Rmd" =
            file.mtime(f["gene"]) > file.mtime(f["model"]))

md   <- readRDS(f["data"])
shap <- read_csv(f["gene"], show_col_types = FALSE)
stopifnot(identical(md$suffix, suffix), setequal(setdiff(names(shap), c("gene_id_clean", "symbol")), md$feature_cols),
          !anyDuplicated(md$test$gene_id_clean), all(shap$gene_id_clean %in% md$test$gene_id_clean))
vals <- md$test[match(shap$gene_id_clean, md$test$gene_id_clean), md$feature_cols]

long <- shap %>% dplyr::select(gene_id_clean, all_of(md$feature_cols)) %>%
  pivot_longer(-gene_id_clean, names_to = "feature", values_to = "shap") %>%
  left_join(bind_cols(gene_id_clean = shap$gene_id_clean, vals) %>%
              mutate(across(-gene_id_clean, as.numeric)) %>%
              pivot_longer(-gene_id_clean, names_to = "feature", values_to = "value"),
            by = c("gene_id_clean", "feature")) %>%
  group_by(feature) %>%
  mutate(value_pct = percent_rank(value)) %>%
  ungroup()
stopifnot("feature value missing for an explained gene" = !anyNA(long$value), nrow(long) == nrow(shap) * length(md$feature_cols))

ranking <- long %>% group_by(feature) %>% summarise(mean_abs_shap = mean(abs(shap))) %>%
  slice_max(mean_abs_shap, n = top_n, with_ties = FALSE)
# must agree with the importance table 03_ wrote, or this is explaining something else
ref <- read_csv(file.path(O, paste0("shap_feature_", suffix, ".csv")), show_col_types = FALSE)
stopifnot("mean |SHAP| does not reproduce 03_'s feature table" =
            isTRUE(all.equal(ranking$mean_abs_shap, ref$mean_abs_shap[match(ranking$feature, ref$feature)])))

what <- sub("_lfc.*", "", sub("^eif33d_promotes_", "si3d ", suffix))
p <- long %>% filter(feature %in% ranking$feature) %>%
  mutate(feature = factor(feature, levels = rev(ranking$feature))) %>%
  ggplot(aes(shap, feature, colour = value_pct)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_quasirandom(orientation = "y", size = 0.8, width = 0.35, alpha = 0.8) +
  scale_colour_gradient2(low = "#0072B2", mid = "grey85", high = "#D55E00", midpoint = 0.5,
                         breaks = c(0, 1), labels = c("low", "high"), name = "feature value") +
  guides(colour = guide_colourbar(barheight = unit(8, "lines"), barwidth = unit(0.6, "lines"))) +
  labs(x = "SHAP value (impact on predicted probability of eIF3d-promoted)", y = NULL,
       title = sprintf("Top %d features: SHAP summary", nrow(ranking)),
       subtitle = sprintf("Classifier, TreeSHAP | %s | n = %d held-out genes explained, n = %d features\nRows ordered by mean |SHAP|; one point per gene",
                          what, nrow(shap), length(md$feature_cols)),
       caption = "Colour = the gene's value of that feature, as a percentile within the feature.") +
  theme_classic(base_size = 12) +
  theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        plot.caption.position = "plot", plot.title.position = "plot")

out_pdf <- here("plots", paste0("03j_shap_beeswarm_", suffix, ".pdf"))
ggsave(out_pdf, p, width = 9, height = 0.3 * nrow(ranking) + 2.5)
cat("wrote", out_pdf, "\n")
