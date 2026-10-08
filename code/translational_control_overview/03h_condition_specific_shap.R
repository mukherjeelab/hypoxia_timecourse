# 03h_condition_specific_shap.R - "hypoxia-specific" SHAP importance: % hypoxia - % normoxia.
#
# A SHAP translation of notebook 09's `hypoxia_score_barplot_top30`, which was built on Gini
# importance. Same quantity and same four-colour concordance scheme; the only change is that
# importance share and direction come from SHAP instead of from impurity:
#
#   share_c  = 100 * mean|SHAP|_c / sum(mean|SHAP|_c)      (per condition, so the two are
#                                                           comparable despite different scales)
#   delta    = share_hypoxia - share_normoxia              <- the bar
#   signed_c = share_c * sign(shap_value_cor_c)            <- direction, for the colour only
#
# Shares are per-condition normalised because raw mean|SHAP| depends on how separable each
# model's classes are; a difference of raw values would mix that in with the thing being asked.
#
# The colour is the CONCORDANCE QUADRANT, not the sign of the bar: a feature can matter more in
# hypoxia while pushing the same way in both conditions.
#
#   Rscript code/translational_control_overview/03h_condition_specific_shap.R [top_n] [engine] [label]
suppressPackageStartupMessages({ library(tidyverse); library(here) })

args   <- commandArgs(trailingOnly = TRUE)
top_n  <- if (length(args) >= 1) as.integer(args[1]) else 30L
engine <- if (length(args) >= 2) args[2] else "kernelshap"
flabel <- if (length(args) >= 3) args[3] else "tco4e"
stopifnot("engine must be treeshap or kernelshap" = engine %in% c("treeshap", "kernelshap"))
etag <- if (engine == "treeshap") "" else paste0("_", engine)
O <- here("output", "translational_control_overview")

suffix_of <- function(cond, tp) sprintf("eif33d_promotes_%s_%s_lfc0.5_%s_corrected_both_clipaggregate",
                                        cond, tp, flabel)
SFX <- c(hypoxia = suffix_of("hypoxia", "1hr"), normoxia = suffix_of("normoxia", "1and4hr"))

read_shap <- function(sfx) {
  f <- file.path(O, paste0("shap_feature_", sfx, etag, ".csv"))
  stopifnot("SHAP output missing - render 03_shap.Rmd for this model and engine first" = file.exists(f))
  m <- file.path(O, paste0("rf_classifier_", sfx, ".rds"))
  stopifnot("SHAP output predates its model - re-render 03_shap.Rmd" = file.mtime(f) > file.mtime(m))
  read_csv(f, show_col_types = FALSE) %>%
    dplyr::select(feature, mean_abs_shap, shap_value_cor) %>%
    mutate(share = 100 * mean_abs_shap / sum(mean_abs_shap),
           signed = share * sign(shap_value_cor))
}
h <- read_shap(SFX[["hypoxia"]]); n <- read_shap(SFX[["normoxia"]])
stopifnot("the two models do not share a feature set - a delta would be meaningless" =
            setequal(h$feature, n$feature))

QUAD <- c("Up in both\n(same direction)"   = "#009E73",
          "Down in both\n(same direction)" = "#D55E00",
          "Discordant\n(hyp+, nor-)"       = "#CC79A7",
          "Discordant\n(hyp-, nor+)"       = "#0072B2")
ann <- read_csv(file.path(O, "feature_annotation.csv"), show_col_types = FALSE)
d <- h %>% dplyr::select(feature, share_hyp = share, signed_hyp = signed, dir_hyp = shap_value_cor) %>%
  inner_join(n %>% dplyr::select(feature, share_nor = share, signed_nor = signed, dir_nor = shap_value_cor),
             by = "feature") %>%
  left_join(ann %>% dplyr::select(feature, block, family), by = "feature") %>%
  mutate(delta = share_hyp - share_nor,
         quadrant = factor(case_when(
           signed_nor >= 0 & signed_hyp >= 0 ~ names(QUAD)[1],
           signed_nor <  0 & signed_hyp <  0 ~ names(QUAD)[2],
           signed_nor <  0 & signed_hyp >= 0 ~ names(QUAD)[3],
           TRUE                              ~ names(QUAD)[4]), levels = names(QUAD)))
stopifnot("a feature has no quadrant" = !anyNA(d$quadrant))

cat(sprintf("%d features | engine %s | %d more important in hypoxia, %d in normoxia\n",
            nrow(d), engine, sum(d$delta > 0), sum(d$delta < 0)))
print(as.data.frame(dplyr::count(d, quadrant)))
cat("\ntop", top_n, "by |hypoxia - normoxia| share:\n")
top <- d %>% slice_max(abs(delta), n = top_n, with_ties = FALSE) %>% arrange(desc(delta))
print(as.data.frame(top %>% transmute(feature, block, share_hyp = round(share_hyp, 2),
                                      share_nor = round(share_nor, 2), delta = round(delta, 2),
                                      quadrant = sub("\n", " ", quadrant))), row.names = FALSE)

p <- top %>% mutate(feature = fct_reorder(feature, -delta)) %>%
  ggplot(aes(feature, delta, fill = quadrant)) +
  geom_col(colour = "grey20", linewidth = 0.25) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = QUAD, name = NULL, drop = FALSE) +
  labs(x = NULL, y = "Hypoxia-specific SHAP importance\n(% hypoxia - % normoxia)",
       title = sprintf("Feature importance: hypoxia-specific score (top %d by |score|)", nrow(top)),
       subtitle = sprintf("si3d | hypoxia 1hr vs normoxia 1+4hr | %s on the probability forest | %d features",
                          if (engine == "kernelshap") "Kernel SHAP" else "TreeSHAP (surrogate)", nrow(d)),
       caption = paste("Bar = difference in share of total mean |SHAP| between the two models, each normalised within its own condition.",
                       "Colour = whether the feature pushes the SAME way in both conditions, from sign(cor(value, its own SHAP)); it is not the sign of the bar.",
                       sep = "\n")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        legend.position = "right", plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30"),
        plot.caption.position = "plot")
out <- here("plots", paste0("03h_hypoxia_specific_shap_", flabel, etag, ".pdf"))
ggsave(out, p, width = 13, height = 6)
write_csv(d %>% arrange(desc(delta)), file.path(O, paste0("shap_hypoxia_specific_", flabel, etag, ".csv")))
cat("\nwrote", out, "\n")
