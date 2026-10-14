# 03i_signed_shap_bar_ci.R - bar versions of 03g's SHAP lollipop, with an interval per bar.
#   bar      = mean |SHAP| on the reference model (neg_seed 9, the model 03g plots)
#   interval = empirical 5-95% of the same quantity across the 20 negative draws of
#              03d_shap_seed_sweep.R. Each draw keeps the same positives, split, forest seed and
#              CV and changes only which negatives are sampled, so the interval is how far the
#              value moves for no biological reason - the refit-based spread, which is the usual
#              way uncertainty is put on feature importance. It is NOT an interval over genes:
#              TreeSHAP is exact and the 200 explained genes are fixed within a draw.
# Two figures:
#   unsigned - the standard mean |SHAP| bar (direction is read off 03j's beeswarm)
#   signed   - 03g's convention, x sign(shap_value_cor), with the sign taken PER DRAW, so a
#              feature whose direction flips between draws gets whiskers that cross zero
# 5-95% rather than a normal-theory interval: values are bounded and skewed for small features.
# Sweeps: 03d_ (tco, TreeSHAP) and 03l_ (tco4e, Kernel SHAP). Under kernelshap the SHAP files
# carry an _kernelshap tag after the seed; the bar is still the neg_seed 9 reference model.
#   Rscript code/translational_control_overview/03i_signed_shap_bar_ci.R [suffix] [top_n] [engine]
# Interactively: edit the settings block below and run top to bottom.
# Colour (2026-10-09): bar fill = feature category from feature_colours.R (shared with 03j_).
# In the signed plot a FADED bar means its direction is not settled: the sign flips in at least
# one of the 20 draws. Bars are outlined so yellow (3'UTR) stays visible, faded or not.
suppressPackageStartupMessages({ library(tidyverse); library(here) })
source(here("code", "translational_control_overview", "feature_colours.R"))

## ---- settings: edit these when running interactively (Positron / RStudio) ----
# Command-line arguments, when given, override them. Defaults match 03g's lollipop and 03j's
# beeswarm: the eIF4E-inclusive hypoxia 1hr classifier explained with Kernel SHAP.
suffix <- "eif33d_promotes_hypoxia_1hr_lfc0.5_tco4e_corrected_both_clipaggregate"
top_n  <- 20L
engine <- "kernelshap"   # or "treeshap"
# what the bar shows: "mean" = mean over the 20 draws (default, Kate 2026-10-09), "median", or
# "reference" = the neg_seed 9 model alone (the one 03j_'s beeswarm explains; extreme for some
# features - 6 of the top 20 on the tco model sit outside their own 5-95%).
bar_stat <- "mean"
# what the whisker shows: "spread" = 5-95% of the 20 draws (how far one run can land), or
# "ci_mean" = 95% t-interval of the 20-draw mean. ci_mean shrinks with more draws and treats the
# draws as independent although they share ~half their negatives, so it is optimistic.
whisker <- "spread"
## ---------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) suffix <- args[1]
if (length(args) >= 2) top_n  <- as.integer(args[2])
if (length(args) >= 3) engine <- args[3]
if (length(args) >= 4) bar_stat <- args[4]
if (length(args) >= 5) whisker  <- args[5]
stopifnot("engine must be treeshap or kernelshap" = engine %in% c("treeshap", "kernelshap"),
          "bar_stat must be mean, median or reference" = bar_stat %in% c("mean", "median", "reference"),
          "whisker must be spread or ci_mean" = whisker %in% c("spread", "ci_mean"))
etag <- if (engine == "treeshap") "" else paste0("_", engine)
O <- here("output", "translational_control_overview")

# the sweep writes neg_seed 9 under the bare suffix, every other seed under _negseed{s};
# 03_shap.Rmd appends the engine tag after that
seed_file <- function(s) file.path(O, paste0("shap_feature_", suffix,
                                             if (s == 9) "" else paste0("_negseed", s), etag, ".csv"))
files <- setNames(map_chr(1:20, seed_file), 1:20)
stopifnot("SHAP seed sweep incomplete - run 03d_ (tco) or 03l_ (tco4e kernelshap)" = all(file.exists(files)))
stopifnot("SHAP output predates the model it claims to explain - re-render 03_shap.Rmd" =
            file.mtime(files[["9"]]) > file.mtime(file.path(O, paste0("rf_classifier_", suffix, ".rds"))))

# shap_value_cor is NA when a feature is constant over a draw's explained genes (only the binary
# is_dap5_repressed / kozak_optimal, bottom of the ranking): no direction, so signed = 0
sweep <- imap_dfr(files, ~ read_csv(.x, show_col_types = FALSE) %>%
                    dplyr::select(feature, mean_abs_shap, shap_value_cor, block, family) %>%
                    mutate(seed = as.integer(.y))) %>%
  mutate(signed = mean_abs_shap * coalesce(sign(shap_value_cor), 0))
stopifnot("draws do not share one feature set" =
            all(table(sweep$feature) == 20), n_distinct(table(sweep$seed)) == 1)
n_feat <- n_distinct(sweep$feature)
n_expl <- nrow(read_csv(file.path(O, paste0("shap_per_gene_", suffix, etag, ".csv")), show_col_types = FALSE))

summ <- sweep %>%
  group_by(feature, block, family) %>%
  summarise(abs_lo = quantile(mean_abs_shap, 0.05), abs_hi = quantile(mean_abs_shap, 0.95),
            signed_lo = quantile(signed, 0.05), signed_hi = quantile(signed, 0.95),
            n_right = sum(signed > 0), n_left = sum(signed < 0),
            med_abs = median(mean_abs_shap), med_signed = median(signed),
            mean_abs = mean(mean_abs_shap), mean_signed = mean(signed),
            # 95% t-interval of the 20-draw mean, used when whisker = "ci_mean"
            abs_se = sd(mean_abs_shap) / sqrt(n()), signed_se = sd(signed) / sqrt(n()),
            .groups = "drop") %>%
  inner_join(sweep %>% filter(seed == 9) %>% dplyr::select(feature, ref_abs = mean_abs_shap, ref_signed = signed),
             by = "feature") %>%
  # the plotted bar; ranking and the top-N cut follow it too
  mutate(mean_abs_shap = switch(bar_stat, mean = mean_abs, median = med_abs, reference = ref_abs),
         signed        = switch(bar_stat, mean = mean_signed, median = med_signed, reference = ref_signed)) %>%
  { if (whisker == "ci_mean")
      mutate(., abs_lo = mean_abs - qt(0.975, 19) * abs_se, abs_hi = mean_abs + qt(0.975, 19) * abs_se,
             signed_lo = mean_signed - qt(0.975, 19) * signed_se, signed_hi = mean_signed + qt(0.975, 19) * signed_se)
    else . } %>%
  mutate(category = feature_category(feature, read_csv(file.path(O, "feature_annotation.csv"),
                                                        show_col_types = FALSE)),
         sign_stable = n_right == 20 | n_left == 20)

top <- summ %>% slice_max(mean_abs_shap, n = top_n, with_ties = FALSE)
stopifnot("a top feature has no direction on the plotted bar" = all(top$signed != 0))
cat(sprintf("top %d of %d features | %d genes explained | 20 draws\n", nrow(top), n_feat, n_expl))
cat(sprintf("top features whose sign flips across draws: %s\n",
            paste(top$feature[!top$sign_stable], collapse = ", ")))

top <- top %>% mutate(direction = factor(if_else(sign_stable, "same sign in all 20 draws",
                                                  "sign flips between draws"),
                                          levels = c("same sign in all 20 draws", "sign flips between draws")))
what <- sub("_lfc.*", "", sub("^eif33d_promotes_", "si3d ", suffix))
subtitle <- sprintf("Classifier, %s | %s\nn = %d held-out genes explained, n = %d features\nBar = %s; whiskers = %s",
                    if (engine == "kernelshap") "Kernel SHAP on the probability forest" else "TreeSHAP",
                    what, n_expl, n_feat,
                    switch(bar_stat, mean = "mean of the 20 negative draws", median = "median of the 20 negative draws",
                           reference = "reference model (neg_seed 9)"),
                    if (whisker == "ci_mean") "95% CI of that mean (t, 19 df)" else "5-95% of the 20 draws")
base_theme <- theme_classic(base_size = 12) +
  theme(legend.position = "right", plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        plot.caption.position = "plot", plot.title.position = "plot")

p_abs <- top %>%
  mutate(feature = fct_reorder(feature, mean_abs_shap)) %>%
  ggplot(aes(mean_abs_shap, feature, fill = category)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.25) +
  geom_errorbarh(aes(xmin = abs_lo, xmax = abs_hi), height = 0.35, linewidth = 0.4, colour = "grey20") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = FEATURE_COLOURS, name = NULL, drop = TRUE) +
  labs(x = "mean |SHAP|", y = NULL,
       title = sprintf("Top %d features: SHAP importance", nrow(top)), subtitle = subtitle,
       caption = if (whisker == "ci_mean") "Whiskers: precision of the 20-draw mean; draws share ~half their negatives, so this is optimistic."
                 else "Whiskers: spread when only the sampled negative genes change (same positives, split and forest seed).") +
  base_theme

p_signed <- top %>%
  mutate(feature = fct_reorder(feature, signed)) %>%
  ggplot(aes(signed, feature, fill = category, alpha = direction)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.25) +
  # whiskers stay full strength on faded bars, so the interval is always readable
  geom_errorbarh(aes(xmin = signed_lo, xmax = signed_hi), height = 0.35, linewidth = 0.4,
                 colour = "grey20", alpha = 1) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = FEATURE_COLOURS, name = NULL, drop = TRUE) +
  scale_alpha_manual(values = c("same sign in all 20 draws" = 1, "sign flips between draws" = 0.35),
                     name = NULL, drop = TRUE) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
         alpha = guide_legend(order = 2, override.aes = list(fill = "grey40"))) +
  labs(x = "mean |SHAP| signed by direction\n(right = higher values push toward eIF3d-promoted)",
       y = NULL,
       title = sprintf("Top %d features: SHAP importance signed by direction", nrow(top)),
       subtitle = subtitle,
       caption = "Direction = sign of cor(feature value, its SHAP value) over the explained genes, taken per draw. Faded bar = the sign flips in at least one of the 20 draws.") +
  base_theme

print(p_abs); print(p_signed)   # show in the Plots pane when run interactively

h <- 0.25 * nrow(top) + 2.8
tagn <- paste0(etag, "_top", top_n, "_", bar_stat, if (whisker == "ci_mean") "_cimean" else "")   # engine and N in the name so variants never overwrite
out <- c(unsigned = here("plots", paste0("03i_shap_bar_ci_", suffix, tagn, ".pdf")),
         signed   = here("plots", paste0("03i_signed_shap_bar_ci_", suffix, tagn, ".pdf")))
ggsave(out[["unsigned"]], p_abs, width = 8.5, height = h)
ggsave(out[["signed"]], p_signed, width = 8.5, height = h)
write_csv(summ %>% arrange(desc(mean_abs_shap)),
          file.path(O, paste0("shap_importance_ci_", suffix, tagn, ".csv")))
cat("wrote", out, sep = "\n")
