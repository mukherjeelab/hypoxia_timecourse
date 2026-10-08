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
#   Rscript code/translational_control_overview/03i_signed_shap_bar_ci.R [suffix] [top_n]
suppressPackageStartupMessages({ library(tidyverse); library(here) })

args   <- commandArgs(trailingOnly = TRUE)
suffix <- if (length(args) >= 1) args[1] else
  "eif33d_promotes_hypoxia_1hr_lfc0.5_tco_corrected_both_clipaggregate"
top_n  <- if (length(args) >= 2) as.integer(args[2]) else 30L
O <- here("output", "translational_control_overview")

# the sweep writes neg_seed 9 under the bare suffix, every other seed under _negseed{s}
seed_file <- function(s) file.path(O, paste0("shap_feature_", suffix,
                                             if (s == 9) "" else paste0("_negseed", s), ".csv"))
files <- setNames(map_chr(1:20, seed_file), 1:20)
stopifnot("SHAP seed sweep incomplete - run 03d_shap_seed_sweep.R" = all(file.exists(files)))
stopifnot("SHAP output predates the model it claims to explain - re-render 03_shap.Rmd" =
            file.mtime(files[["9"]]) > file.mtime(file.path(O, paste0("rf_classifier_", suffix, ".rds"))))

# shap_value_cor is NA when a feature is constant over a draw's explained genes (only the binary
# is_dap5_repressed / kozak_optimal, bottom of the ranking): no direction, so signed = 0
sweep <- imap_dfr(files, ~ read_csv(.x, show_col_types = FALSE) %>%
                    dplyr::select(feature, mean_abs_shap, shap_value_cor, block) %>%
                    mutate(seed = as.integer(.y))) %>%
  mutate(signed = mean_abs_shap * coalesce(sign(shap_value_cor), 0))
stopifnot("draws do not share one feature set" =
            all(table(sweep$feature) == 20), n_distinct(table(sweep$seed)) == 1)
n_feat <- n_distinct(sweep$feature)
n_expl <- nrow(read_csv(file.path(O, paste0("shap_per_gene_", suffix, ".csv")), show_col_types = FALSE))

summ <- sweep %>%
  group_by(feature, block) %>%
  summarise(abs_lo = quantile(mean_abs_shap, 0.05), abs_hi = quantile(mean_abs_shap, 0.95),
            signed_lo = quantile(signed, 0.05), signed_hi = quantile(signed, 0.95),
            n_right = sum(signed > 0), n_left = sum(signed < 0), .groups = "drop") %>%
  inner_join(sweep %>% filter(seed == 9) %>% dplyr::select(feature, mean_abs_shap, signed),
             by = "feature") %>%
  mutate(category = if_else(block == "intrinsic", "intrinsic mRNA feature", "external dataset"),
         sign_stable = n_right == 20 | n_left == 20)

top <- summ %>% slice_max(mean_abs_shap, n = top_n, with_ties = FALSE)
stopifnot("a top feature has no direction on the reference model" = all(top$signed != 0))
cat(sprintf("top %d of %d features | %d genes explained | 20 draws\n", nrow(top), n_feat, n_expl))
cat(sprintf("top features whose sign flips across draws: %s\n",
            paste(top$feature[!top$sign_stable], collapse = ", ")))

category_colors <- c("external dataset" = "#CC79A7", "intrinsic mRNA feature" = "#3D9C60")
what <- sub("_lfc.*", "", sub("^eif33d_promotes_", "si3d ", suffix))
subtitle <- sprintf("Classifier, TreeSHAP | %s | n = %d held-out genes explained, n = %d features\nBar = reference model (neg_seed 9); whiskers = 5-95%% over 20 negative draws",
                    what, n_expl, n_feat)
base_theme <- theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        plot.caption.position = "plot", plot.title.position = "plot")

p_abs <- top %>%
  mutate(feature = fct_reorder(feature, mean_abs_shap)) %>%
  ggplot(aes(mean_abs_shap, feature, fill = category)) +
  geom_col(width = 0.7) +
  geom_errorbarh(aes(xmin = abs_lo, xmax = abs_hi), height = 0.35, linewidth = 0.4, colour = "grey20") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = category_colors, name = NULL) +
  labs(x = "mean |SHAP|", y = NULL,
       title = sprintf("Top %d features: SHAP importance", nrow(top)), subtitle = subtitle,
       caption = "Whiskers: spread when only the sampled negative genes change (same positives, split and forest seed).") +
  base_theme

p_signed <- top %>%
  mutate(feature = fct_reorder(feature, signed)) %>%
  ggplot(aes(signed, feature, fill = category)) +
  geom_col(width = 0.7) +
  geom_errorbarh(aes(xmin = signed_lo, xmax = signed_hi), height = 0.35, linewidth = 0.4, colour = "grey20") +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = category_colors, name = NULL) +
  labs(x = "mean |SHAP| signed by direction\n(right = higher values push toward eIF3d-promoted)",
       y = NULL,
       title = sprintf("Top %d features: SHAP importance signed by direction", nrow(top)),
       subtitle = subtitle,
       caption = "Direction = sign of cor(feature value, its SHAP value) over the explained genes, taken per draw; whiskers crossing 0 mean the sign flips between draws.") +
  base_theme

h <- 0.25 * nrow(top) + 2.8
out <- c(unsigned = here("plots", paste0("03i_shap_bar_ci_", suffix, ".pdf")),
         signed   = here("plots", paste0("03i_signed_shap_bar_ci_", suffix, ".pdf")))
ggsave(out[["unsigned"]], p_abs, width = 8.5, height = h)
ggsave(out[["signed"]], p_signed, width = 8.5, height = h)
write_csv(summ %>% arrange(desc(mean_abs_shap)),
          file.path(O, paste0("shap_importance_ci_", suffix, ".csv")))
cat("wrote", out, sep = "\n")
