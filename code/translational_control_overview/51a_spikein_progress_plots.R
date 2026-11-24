# 51a_spikein_progress_plots.R - interim spike-in vs no-spike figures from whatever splits of the
# `spikein` sweep exist so far (11d_shap_regression_seed_sweep.R spikein). Safe to re-run at any
# point; every panel prints how many splits each arm contributes. The final, paired analysis is
# 51_spikein_regression_comparison.Rmd, which refuses to run until all 20 splits exist.
#
#   Rscript code/translational_control_overview/51a_spikein_progress_plots.R
#
# Writes plots/51a_spikein_progress_rho2.pdf and plots/51a_spikein_progress_top15_shap.pdf.
suppressPackageStartupMessages({ library(tidyverse); library(here) })
src <- here("code", "translational_control_overview")
source(file.path(src, "02_sweep_presets.R")); source(file.path(src, "10_regression_presets.R"))
O <- here("output", "translational_control_overview")
arms <- REG_PRESETS$spikein
fpath <- function(stem, arm, s) file.path(O, paste0(stem, reg_suffix(arms[[arm]], s), ".csv"))
arm_lab  <- c(spike = "Spike-in", nospike = "No spike-in")
arm_fill <- c("Spike-in" = "#E69F00", "No spike-in" = "#87610e")

# --- which splits exist ------------------------------------------------------------------------
have <- expand_grid(arm = names(arms), seed = REG_SPLIT_SEEDS) %>%
  mutate(perf = map2_lgl(arm, seed, ~ file.exists(fpath("rfreg_performance_", .x, .y))),
         shap = map2_lgl(arm, seed, ~ file.exists(fpath("shapreg_feature_", .x, .y))))

# A split counts as paired only if both arms put the SAME genes in train and test.
paired_seeds <- have %>% filter(perf) %>% dplyr::count(seed) %>% filter(n == 2) %>% pull(seed)
for (s in paired_seeds) {
  a <- read_csv(fpath("rfreg_predictions_", "spike", s), show_col_types = FALSE)
  b <- read_csv(fpath("rfreg_predictions_", "nospike", s), show_col_types = FALSE)
  j <- inner_join(a, b, by = "gene_id_clean")
  stopifnot("a 'paired' split has different train/test genes in the two arms" =
              nrow(j) == nrow(a) && all(j$split.x == j$split.y))
}
cat("splits with performance:", paste(sprintf("%s %d", arm_lab[names(arms)],
    map_int(names(arms), ~ sum(have$perf[have$arm == .x]))), collapse = ", "),
    "| with SHAP:", paste(sprintf("%s %d", arm_lab[names(arms)],
    map_int(names(arms), ~ sum(have$shap[have$arm == .x]))), collapse = ", "),
    "| paired:", paste(paired_seeds, collapse = ","), "\n")

# --- 1. rho2 -----------------------------------------------------------------------------------
perf <- have %>% filter(perf) %>%
  mutate(rho = map2_dbl(arm, seed, ~ read_csv(fpath("rfreg_performance_", .x, .y),
                                              show_col_types = FALSE) %>%
                          filter(metric == "spearman_raw_y") %>% pull(value)),
         rho2 = rho^2, arm = factor(arm_lab[arm], levels = arm_lab))
n_by <- perf %>% dplyr::count(arm)
lab <- perf %>% group_by(arm) %>%
  summarise(txt = sprintf("n = %d split%s\nmedian rho2 = %.3f", n(), if (n() == 1) "" else "s",
                          median(rho2)), y = max(rho2), .groups = "drop")
p1 <- ggplot(perf, aes(arm, rho2)) +
  geom_line(data = perf %>% filter(seed %in% paired_seeds), aes(group = seed),
            colour = "grey55", linewidth = 0.5) +
  geom_point(aes(fill = arm), shape = 21, size = 2.6, colour = "black",
             position = position_jitter(width = 0.04, height = 0, seed = 9)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.35, linewidth = 0.4) +
  geom_text(data = lab, aes(x = arm, y = y, label = txt), vjust = -0.4, size = 3) +
  scale_fill_manual(values = arm_fill, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "held-out rho2 (Spearman squared)",
       title = "si3d hypoxia 1hr intrinsic regression",
       subtitle = paste0("PRELIMINARY: sweep in progress\n",
                         "point = one train/test split; lines join the\n",
                         "same split in both arms (", length(paired_seeds), " paired so far)")) +
  theme_classic(base_size = 12) +
  theme(plot.title.position = "plot", plot.subtitle = element_text(size = 9.5))
ggsave(here("plots", "51a_spikein_progress_rho2.pdf"), p1, width = 4.6, height = 4.4)

# --- 2. top-15 SHAP ----------------------------------------------------------------------------
feat <- have %>% filter(shap) %>%
  mutate(d = map2(arm, seed, ~ read_csv(fpath("shapreg_feature_", .x, .y), show_col_types = FALSE) %>%
                    dplyr::select(feature, share_pct, shap_value_cor))) %>%
  unnest(d) %>% mutate(arm = factor(arm_lab[arm], levels = arm_lab))
summ <- feat %>% group_by(arm, feature) %>%
  summarise(n = n(), med = median(share_pct),
            lo = if (n() > 2) quantile(share_pct, 0.05) else NA_real_,
            hi = if (n() > 2) quantile(share_pct, 0.95) else NA_real_,
            dir = if_else(median(shap_value_cor, na.rm = TRUE) > 0, "+", "-"), .groups = "drop")
top <- summ %>% group_by(arm) %>% slice_max(med, n = 15) %>% pull(feature) %>% unique()
ord <- summ %>% filter(arm == "Spike-in", feature %in% top) %>% arrange(med) %>% pull(feature)
ord <- c(setdiff(top, ord), ord)        # any no-spike-only top feature goes at the bottom
dirs <- summ %>% filter(feature %in% top) %>%
  dplyr::select(arm, feature, dir) %>% pivot_wider(names_from = arm, values_from = dir)
dir_lab <- setNames(sprintf("%s  (%s / %s)", dirs$feature, coalesce(dirs$`Spike-in`, "?"),
                            coalesce(dirs$`No spike-in`, "?")), dirs$feature)
t15 <- summ %>% filter(feature %in% top) %>%
  mutate(feature = factor(feature, levels = ord),
         arm = factor(arm, levels = rev(arm_lab)))
n_shap <- feat %>% distinct(arm, seed) %>% dplyr::count(arm)
p2 <- ggplot(t15, aes(med, feature, fill = arm)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.8),
                 height = 0.25, linewidth = 0.35, na.rm = TRUE) +
  scale_y_discrete(labels = dir_lab) +
  scale_fill_manual(values = arm_fill, breaks = arm_lab,
                    labels = sprintf("%s (n = %d split%s)", n_shap$arm, n_shap$n,
                                     if_else(n_shap$n == 1, "", "s"))[match(arm_lab, n_shap$arm)]) +
  labs(x = "share of total |SHAP| (%)", y = NULL, fill = NULL,
       title = "Top-15 features by SHAP: spike-in vs no spike-in",
       subtitle = paste0("PRELIMINARY: sweep in progress. Bar = median over splits; whiskers 5-95% (> 2 splits only)\n",
                         "200 held-out genes per split. (+/-) = direction of effect, spike-in / no spike-in")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", plot.title.position = "plot",
        plot.subtitle = element_text(size = 9))
ggsave(here("plots", "51a_spikein_progress_top15_shap.pdf"), p2, width = 8, height = 6.5)

print(perf %>% dplyr::select(arm, seed, rho, rho2) %>% arrange(arm, seed), n = Inf)
print(summ %>% filter(feature %in% top) %>% dplyr::select(arm, feature, n, med, dir) %>%
        pivot_wider(names_from = arm, values_from = c(n, med, dir)) %>%
        arrange(desc(`med_Spike-in`)), n = Inf, width = Inf)
