# 03j_shap_beeswarm.R - standard SHAP summary (beeswarm) plot for the top-N features.
# One point per explained gene per feature: x = that gene's SHAP value, colour = the gene's value
# of the feature, as a within-feature percentile so every row uses the full colour range
# (the shap package's convention). Direction is read directly: high-value (vermillion) points on
# the right mean high values push toward eIF3d-promoted. Rows ordered by mean |SHAP|, so the
# ranking matches 03g / 03i.
# Feature values come from the model's saved, imputed test split (02_'s rf_model_data), i.e. the
# exact inputs 03_ explained, matched by gene.
#   Rscript code/translational_control_overview/03j_shap_beeswarm.R [suffix] [top_n] [engine]
# Interactively: edit the settings block below and run top to bottom.
# Colour (2026-10-09): each row runs from grey (low feature value) to its category colour (high),
# categories from feature_colours.R, shared with 03i_ so a feature is one colour everywhere.
# Dots are outlined so the yellow (3'UTR) rows stay visible on white.
suppressPackageStartupMessages({ library(tidyverse); library(here); library(ggbeeswarm) })
source(here("code", "translational_control_overview", "feature_colours.R"))

## ---- settings: edit these when running interactively (Positron / RStudio) ----
# Command-line arguments, when given, override them. The defaults match 03g's lollipop: the
# eIF4E-inclusive hypoxia 1hr classifier explained with Kernel SHAP.
suffix <- "eif33d_promotes_hypoxia_1hr_lfc0.5_tco4e_corrected_both_clipaggregate"
top_n  <- 12L
# "kernelshap" = Shapley values of the actual probability forest; "treeshap" = the surrogate's
engine <- "kernelshap"
# "category" = each row grey -> its feature_colours.R category colour (2026-10-09);
# "value"    = the original palette: every row dark blue (low) - grey - vermillion (high)
palette <- "category"
# split_n > 0: two panels by feature DIRECTION (sign of shap_value_cor - do higher values push
# toward or away from eIF3d-promoted), each showing its own top split_n by mean |SHAP|;
# top_n is then ignored. 0 = one panel, top_n overall.
split_n <- 0L
# fixed x-axis range, e.g. c(-0.15, 0.10), so hypoxia and normoxia figures share one scale;
# NULL = each figure's own range. Stops if a point would fall outside.
x_limits <- NULL
## ---------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) suffix <- args[1]
if (length(args) >= 2) top_n  <- as.integer(args[2])
if (length(args) >= 3) engine <- args[3]
if (length(args) >= 4) palette <- args[4]
if (length(args) >= 5) split_n <- as.integer(args[5])
if (length(args) >= 6) x_limits <- as.numeric(strsplit(args[6], ",")[[1]])
stopifnot("engine must be treeshap or kernelshap" = engine %in% c("treeshap", "kernelshap"),
          "palette must be category or value" = palette %in% c("category", "value"))
# 03_shap.Rmd tags its SHAP outputs with the engine; the MODEL files are shared between engines
etag <- if (engine == "treeshap") "" else paste0("_", engine)
O <- here("output", "translational_control_overview")

f <- c(model = file.path(O, paste0("rf_classifier_", suffix, ".rds")),
       data  = file.path(O, paste0("rf_model_data_", suffix, ".rds")),
       gene  = file.path(O, paste0("shap_per_gene_", suffix, etag, ".csv")))
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

ref <- read_csv(file.path(O, paste0("shap_feature_", suffix, etag, ".csv")), show_col_types = FALSE)
PANELS <- c("Higher values push toward eIF3d-promoted", "Higher values push away from eIF3d-promoted")
ranking <- long %>% group_by(feature) %>% summarise(mean_abs_shap = mean(abs(shap))) %>%
  left_join(ref %>% dplyr::select(feature, shap_value_cor), by = "feature") %>%
  mutate(panel = factor(if_else(shap_value_cor > 0, PANELS[1], PANELS[2]), levels = PANELS))
ranking <- if (split_n > 0) {
  # a feature constant over the explained genes (e.g. a binary flag) has no direction and fits
  # neither panel; it is dropped, reported, and asserted not to be a top feature overall
  no_dir <- ranking %>% filter(is.na(shap_value_cor))
  if (nrow(no_dir)) cat("no direction, left out of the panels:", paste(no_dir$feature, collapse = ", "), "\n")
  stopifnot("a feature with no direction ranks among the top - it would be silently missing" =
              all(rank(-ranking$mean_abs_shap)[is.na(ranking$shap_value_cor)] > 2 * split_n))
  ranking %>% filter(!is.na(shap_value_cor)) %>% group_by(panel) %>% slice_max(mean_abs_shap, n = split_n, with_ties = FALSE) %>%
    ungroup() %>% arrange(panel, desc(mean_abs_shap))
} else {
  ranking %>% slice_max(mean_abs_shap, n = top_n, with_ties = FALSE)
}
# must agree with the importance table 03_ wrote, or this is explaining something else
stopifnot("mean |SHAP| does not reproduce 03_'s feature table" =
            isTRUE(all.equal(ranking$mean_abs_shap, ref$mean_abs_shap[match(ranking$feature, ref$feature)])))

ann <- read_csv(file.path(O, "feature_annotation.csv"), show_col_types = FALSE)
ranking <- ranking %>% mutate(category = feature_category(feature, ann))
# per-point fill: interpolate grey -> category colour by the gene's within-feature percentile
ramp_fill <- function(pct, cat) {
  out <- character(length(pct))
  for (k in unique(as.character(cat))) {
    i <- which(as.character(cat) == k)
    m <- grDevices::colorRamp(c(FEATURE_GRADIENT_LOW[[k]], FEATURE_COLOURS[[k]]))(pct[i])
    out[i] <- grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
  }
  out
}
plot_df <- long %>% filter(feature %in% ranking$feature) %>%
  left_join(ranking %>% dplyr::select(feature, category, panel), by = "feature") %>%
  mutate(feature = factor(feature, levels = rev(ranking$feature)),
         fill_hex = ramp_fill(value_pct, category))
stopifnot("a point has no fill colour" = !anyNA(plot_df$fill_hex))

what <- sub("_lfc.*", "", sub("^eif33d_promotes_", "si3d ", suffix))
colour_layers <- if (palette == "category") list(
  geom_quasirandom(aes(fill = fill_hex), orientation = "y", shape = 21, size = 1.3,
                   stroke = 0.15, colour = "grey20", width = 0.35),
  scale_fill_identity(),
  # invisible layer that carries the category legend
  geom_point(aes(colour = category), data = distinct(plot_df, feature, category, panel) %>% mutate(shap = 0),
             alpha = 0),
  scale_colour_manual(values = FEATURE_COLOURS, name = NULL, drop = TRUE),
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3.5, shape = 21,
                                                   fill = FEATURE_COLOURS[levels(droplevels(plot_df$category))],
                                                   colour = "grey20", stroke = 0.3)))
) else list(
  geom_quasirandom(aes(colour = value_pct), orientation = "y", size = 0.8, width = 0.35, alpha = 0.8),
  scale_colour_gradient2(low = "#0072B2", mid = "grey85", high = "#D55E00", midpoint = 0.5,
                         breaks = c(0, 1), labels = c("low", "high"), name = "feature value"),
  guides(colour = guide_colourbar(barheight = unit(8, "lines"), barwidth = unit(0.6, "lines")))
)
p <- ggplot(plot_df, aes(shap, feature)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  colour_layers +
  labs(x = "SHAP value (impact on predicted probability of eIF3d-promoted)", y = NULL,
       title = if (split_n > 0) sprintf("Top %d features in each direction: SHAP summary", split_n)
               else sprintf("Top %d features: SHAP summary", nrow(ranking)),
       subtitle = sprintf("Classifier, %s | %s\nn = %d held-out genes explained, n = %d features | rows ordered by mean |SHAP|; one point per gene",
                          if (engine == "kernelshap") "Kernel SHAP on the probability forest" else "TreeSHAP",
                          what, nrow(shap), length(md$feature_cols)),
       caption = if (palette == "category")
         "Each row runs from grey (low value of that feature) to its category colour (high value); value is the gene's percentile within the feature."
       else "Colour = the gene's value of that feature, as a percentile within the feature.") +
  theme_classic(base_size = 12) +
  theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        plot.caption.position = "plot", plot.title.position = "plot")

if (!is.null(x_limits)) {
  stopifnot("x_limits must be c(low, high)" = length(x_limits) == 2 && x_limits[1] < x_limits[2],
            "a point falls outside x_limits - widen them" =
              all(plot_df$shap >= x_limits[1] & plot_df$shap <= x_limits[2]))
  p <- p + scale_x_continuous(limits = x_limits)
}
if (split_n > 0)
  # side by side (Kate 2026-10-09); both panels hold split_n rows, so row spacing matches
  p <- p + facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    theme(strip.text = element_text(hjust = 0, size = 11, face = "bold"),
          strip.background = element_blank())
print(p)   # shows in the Plots pane when run interactively

# top_n is in the name so a top-12 figure does not overwrite the top-30 one
out_pdf <- here("plots", paste0("03j_shap_beeswarm_", suffix, etag,
                                  if (split_n > 0) paste0("_split", split_n) else paste0("_top", top_n),
                                  if (palette == "value") "_valuecolour" else "",
                                  if (!is.null(x_limits)) "_fixedx" else "", ".pdf"))
ggsave(out_pdf, p, width = if (split_n > 0) 15 else 9,
       height = if (split_n > 0) 0.3 * split_n + 3 else 0.3 * nrow(ranking) + 2.5)
cat("wrote", out_pdf, "\n")
