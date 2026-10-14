# 03h_condition_specific_shap.R - "hypoxia-specific" SHAP importance: % hypoxia - % normoxia.
#
# A SHAP translation of notebook 09's `hypoxia_score_barplot_top30`, which was built on Gini
# importance. Importance share and direction come from SHAP instead of from impurity:
#
#   share_c  = 100 * mean|SHAP|_c / sum(mean|SHAP|_c)      (per condition, so the two are
#                                                           comparable despite different scales)
#   delta    = share_hypoxia - share_normoxia              <- the bar
#   signed_c = share_c * sign(shap_value_cor_c)            <- direction, for the fading only
#
# Shares are per-condition normalised because raw mean|SHAP| depends on how separable each
# model's classes are; a difference of raw values would mix that in with the thing being asked.
#
# 2026-10-09 (Kate): colour = the feature's category from feature_colours.R, shared with 03i_ /
# 03j_, replacing the four-colour concordance scheme. Concordance is now the FADE: a faded bar
# pushes opposite ways in the two conditions. Whiskers = 5-95% of delta over all 20 x 20 pairings
# of a hypoxia draw with a normoxia draw (03l_ sweeps). Draws are not paired by seed number - the
# two models have different positives and draw their negatives independently, so seed s in one
# has nothing to do with seed s in the other; every pairing is an equally valid sample.
#
#   Rscript code/translational_control_overview/03h_condition_specific_shap.R [top_n] [engine] [label] [bar_stat]
# Interactively: edit the settings block below and run top to bottom.
suppressPackageStartupMessages({ library(tidyverse); library(here) })
source(here("code", "translational_control_overview", "feature_colours.R"))

## ---- settings: edit these when running interactively (Positron / RStudio) ----
top_n    <- 30L
engine   <- "kernelshap"   # or "treeshap"
flabel   <- "tco4e"
# bar: "mean" (default) = mean over the 400 pairings = mean hypoxia share - mean normoxia share;
#      "median" = median over the pairings; "reference" = delta between the two neg_seed 9 models
bar_stat <- "mean"
# whisker: "spread" = 5-95% of the 400 pairings; "ci_mean" = 95% Welch t-interval of the
#          difference of the two 20-draw means (optimistic - draws share ~half their negatives)
whisker  <- "spread"
## ---------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) top_n    <- as.integer(args[1])
if (length(args) >= 2) engine   <- args[2]
if (length(args) >= 3) flabel   <- args[3]
if (length(args) >= 4) bar_stat <- args[4]
if (length(args) >= 5) whisker  <- args[5]
stopifnot("engine must be treeshap or kernelshap" = engine %in% c("treeshap", "kernelshap"),
          "bar_stat must be mean, median or reference" = bar_stat %in% c("mean", "median", "reference"),
          "whisker must be spread or ci_mean" = whisker %in% c("spread", "ci_mean"))
etag <- if (engine == "treeshap") "" else paste0("_", engine)
O <- here("output", "translational_control_overview")

suffix_of <- function(cond, tp) sprintf("eif33d_promotes_%s_%s_lfc0.5_%s_corrected_both_clipaggregate",
                                        cond, tp, flabel)
SFX <- c(hypoxia = suffix_of("hypoxia", "1hr"), normoxia = suffix_of("normoxia", "1and4hr"))

# all 20 draws of one model, share computed WITHIN each draw
read_sweep <- function(sfx) {
  files <- setNames(file.path(O, paste0("shap_feature_", sfx, ifelse(1:20 == 9, "", paste0("_negseed", 1:20)),
                                        etag, ".csv")), 1:20)
  stopifnot("SHAP seed sweep incomplete - run 03l_shap_seed_sweep_tco4e_kernel.R" = all(file.exists(files)))
  m <- file.path(O, paste0("rf_classifier_", sfx, ".rds"))
  stopifnot("reference SHAP output predates its model - re-render 03_shap.Rmd" =
              file.mtime(files[["9"]]) > file.mtime(m))
  imap_dfr(files, ~ read_csv(.x, show_col_types = FALSE) %>%
             dplyr::select(feature, mean_abs_shap, shap_value_cor) %>%
             mutate(seed = as.integer(.y), share = 100 * mean_abs_shap / sum(mean_abs_shap),
                    signed = share * coalesce(sign(shap_value_cor), 0)))
}
h <- read_sweep(SFX[["hypoxia"]]); n <- read_sweep(SFX[["normoxia"]])
stopifnot("the two models do not share a feature set - a delta would be meaningless" =
            setequal(h$feature, n$feature),
          "a draw is missing features" = all(table(h$feature) == 20), all(table(n$feature) == 20))

# all 400 pairings per feature
pairs <- inner_join(h %>% dplyr::select(feature, seed_h = seed, share_h = share),
                    n %>% dplyr::select(feature, seed_n = seed, share_n = share),
                    by = "feature", relationship = "many-to-many") %>%
  mutate(delta = share_h - share_n)
stopifnot("every feature should have 20 x 20 pairings" = all(dplyr::count(pairs, feature)$n == 400))
ci <- pairs %>% group_by(feature) %>%
  summarise(lo = quantile(delta, 0.05), hi = quantile(delta, 0.95), med_delta = median(delta),
            mean_delta = mean(delta), .groups = "drop")
if (whisker == "ci_mean") {
  welch <- inner_join(h %>% group_by(feature) %>% summarise(mh = mean(share), vh = var(share) / n(), nh = n()),
                      n %>% group_by(feature) %>% summarise(mn = mean(share), vn = var(share) / n(), nn = n()),
                      by = "feature") %>%
    mutate(se = sqrt(vh + vn), df = (vh + vn)^2 / (vh^2 / (nh - 1) + vn^2 / (nn - 1)),
           lo = (mh - mn) - qt(0.975, df) * se, hi = (mh - mn) + qt(0.975, df) * se)
  ci <- ci %>% dplyr::select(-lo, -hi) %>% inner_join(welch %>% dplyr::select(feature, lo, hi), by = "feature")
}

pick <- function(df, cond) {
  ref <- df %>% filter(seed == 9) %>% dplyr::select(feature, share, signed)
  agg <- df %>% group_by(feature) %>%
    summarise(share_med = median(share), signed_med = median(signed),
              share_mean = mean(share), signed_mean = mean(signed))
  ref %>% inner_join(agg, by = "feature") %>%
    transmute(feature,
              !!paste0("share_", cond)  := switch(bar_stat, mean = share_mean, median = share_med, reference = share),
              !!paste0("signed_", cond) := switch(bar_stat, mean = signed_mean, median = signed_med, reference = signed))
}
ann <- read_csv(file.path(O, "feature_annotation.csv"), show_col_types = FALSE)
d <- pick(h, "hyp") %>% inner_join(pick(n, "nor"), by = "feature") %>%
  inner_join(ci, by = "feature") %>%
  mutate(delta = switch(bar_stat, mean = mean_delta, median = med_delta, reference = share_hyp - share_nor),
         category = feature_category(feature, ann),
         direction = factor(if_else(sign(signed_hyp) == sign(signed_nor),
                                    "same direction in both conditions", "opposite directions"),
                            levels = c("same direction in both conditions", "opposite directions")),
         ci_crosses_zero = lo < 0 & hi > 0)

top <- d %>% slice_max(abs(delta), n = top_n, with_ties = FALSE) %>% arrange(desc(delta))
cat(sprintf("%d features | engine %s | bar = %s | %d more important in hypoxia, %d in normoxia\n",
            nrow(d), engine, bar_stat, sum(d$delta > 0), sum(d$delta < 0)))
cat(sprintf("top %d: %d whiskers cross zero, %d push opposite directions\n", nrow(top),
            sum(top$ci_crosses_zero), sum(top$direction == "opposite directions")))
print(as.data.frame(top %>% transmute(feature, category, share_hyp = round(share_hyp, 2),
                                      share_nor = round(share_nor, 2), delta = round(delta, 2),
                                      lo = round(lo, 2), hi = round(hi, 2), direction)), row.names = FALSE)

p <- top %>% mutate(feature = fct_reorder(feature, -delta)) %>%
  ggplot(aes(feature, delta, fill = category, alpha = direction)) +
  geom_col(colour = "grey20", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.35, linewidth = 0.4, colour = "grey20", alpha = 1) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = FEATURE_COLOURS, name = NULL, drop = TRUE) +
  scale_alpha_manual(values = c("same direction in both conditions" = 1, "opposite directions" = 0.35),
                     name = NULL, drop = TRUE) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
         # the fade key only means something when a top feature actually is faded
         alpha = if (n_distinct(top$direction) > 1)
                   guide_legend(order = 2, override.aes = list(fill = "grey40")) else "none") +
  labs(x = NULL, y = "Hypoxia-specific SHAP importance\n(% hypoxia - % normoxia)",
       title = sprintf("Feature importance: hypoxia-specific score (top %d by |score|)", nrow(top)),
       subtitle = sprintf("si3d | hypoxia 1hr vs normoxia 1+4hr | %s on the probability forest | %d features\nBar = %s; whiskers = %s",
                          if (engine == "kernelshap") "Kernel SHAP" else "TreeSHAP (surrogate)", nrow(d),
                          switch(bar_stat, mean = "mean hypoxia - mean normoxia over 20 draws each",
                                 median = "median over the 20 x 20 draw pairings", reference = "the two reference models (neg_seed 9)"),
                          if (whisker == "ci_mean") "95% Welch CI of the difference of the 20-draw means"
                          else "5-95% over all 20 x 20 draw pairings"),
       caption = paste("Bar = difference in share of total mean |SHAP| between the two models, each normalised within its own condition (and within each draw for the whiskers).",
                       if (any(top$direction == "opposite directions"))
                         "Faded = the feature pushes opposite ways in the two conditions, from sign(cor(value, its own SHAP)). Colour = feature category."
                       else "Every feature shown pushes the same way in both conditions (sign of cor(value, its own SHAP), mean over draws). Colour = feature category.",
                       sep = "\n")) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 13), legend.text = element_text(size = 13),
        # title and subtitle at the original base_size-12 sizes; only the data text is enlarged
        plot.title = element_text(size = 14.4), plot.subtitle = element_text(size = 12),
        legend.position = "right", plot.caption = element_text(hjust = 0, size = 10, colour = "grey30"),
        plot.caption.position = "plot")
print(p)   # shows in the Plots pane when run interactively

tagn <- paste0(etag, "_top", top_n, "_", bar_stat, if (whisker == "ci_mean") "_cimean" else "")
out <- here("plots", paste0("03h_hypoxia_specific_shap_ci_", flabel, tagn, ".pdf"))
ggsave(out, p, width = 13, height = 6.5)   # the original long, short shape
write_csv(d %>% arrange(desc(delta)), file.path(O, paste0("shap_hypoxia_specific_ci_", flabel, tagn, ".csv")))
cat("\nwrote", out, "\n")
