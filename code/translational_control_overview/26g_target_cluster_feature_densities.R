# Density of one feature for the si3d hypoxia 1hr targets vs every other modelled gene, then
# per target cluster. Answers "is cluster X actually long / GC-rich, or just less short?" -
# the 26d_ heatmap shows z-scores against ALL modelled genes, which says the direction but not
# where a cluster sits on the absolute scale.
# Clusters: 26b_'s "features, raw" profile at k = 3 (sizes 801 / 410 / 225), the grouping the
# 26d_ heatmap uses. Features: feature_matrix_tco.rds, unimputed, as in 26d_.
# Run: Rscript code/translational_control_overview/26g_target_cluster_feature_densities.R
suppressPackageStartupMessages({library(tidyverse); library(here); library(patchwork)})

O <- here("output", "translational_control_overview")
suffix <- "eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw"

cl <- read_csv(file.path(O, paste0("target_shap_feature_clusters_", suffix, ".csv")), show_col_types = FALSE) %>%
  filter(profile == "features, raw", k == 3) %>% dplyr::select(gene_id_clean, cluster)
stopifnot("cluster sizes are not 801/410/225" = identical(as.integer(table(cl$cluster)), c(801L, 410L, 225L)),
          "duplicate targets" = !anyDuplicated(cl$gene_id_clean))

oof <- readRDS(file.path(O, paste0("shap_oof_", suffix, ".rds")))
modelled <- oof$genes$gene_id_clean
stopifnot("targets not all modelled" = all(cl$gene_id_clean %in% modelled))

fm <- readRDS(file.path(O, "feature_matrix_tco.rds"))
stopifnot("matrix has duplicate genes" = !anyDuplicated(fm$gene_id_clean))
d0 <- fm %>% filter(gene_id_clean %in% modelled) %>%
  transmute(gene_id_clean,
            cds = as.numeric(unlist(cds_length)), utr5 = as.numeric(unlist(utr5_length)),
            gc3 = as.numeric(unlist(gc3_internal))) %>%
  left_join(cl, by = "gene_id_clean") %>%
  mutate(target = !is.na(cluster))
stopifnot("modelled genes lost" = nrow(d0) == length(modelled),
          "GC3 out of [0,100]" = all(d0$gc3 >= 0 & d0$gc3 <= 100, na.rm = TRUE),
          "lengths not positive" = all(d0$cds > 0, na.rm = TRUE) && all(d0$utr5 >= 0, na.rm = TRUE))

cl_cols <- c("Cluster 1" = "#0072B2", "Cluster 2" = "#CC79A7", "Cluster 3" = "#009E73")  # as in 26d_

features <- list(
  cds  = list(col = "cds",  x = "log2 CDS length (nt)",    log = TRUE,  unit = "nt"),
  gc3  = list(col = "gc3",  x = "GC3 content (%)",         log = FALSE, unit = "%"),
  utr5 = list(col = "utr5", x = "log2 5' UTR length (nt)", log = TRUE,  unit = "nt")
)

summ_all <- list()
for (nm in names(features)) {
  f <- features[[nm]]
  d <- d0 %>% mutate(v = .data[[f$col]]) %>% filter(!is.na(v))
  if (f$log) d <- d %>% mutate(xv = log2(v + 1)) else d <- d %>% mutate(xv = v)
  bg <- d %>% filter(!target)
  # Where each group's median falls in the background distribution
  pct <- function(m) 100 * mean(bg$v <= m)
  grp <- bind_rows(
    d %>% filter(target) %>% mutate(group = "All si3d targets"),
    d %>% filter(target) %>% mutate(group = paste("Cluster", cluster))
  )
  summ <- grp %>% group_by(group) %>%
    summarise(n = n(), median = median(v),
              p = wilcox.test(v, bg$v)$p.value, .groups = "drop") %>%
    mutate(bg_percentile = map_dbl(median, pct)) %>%
    add_row(group = "Everything else", n = nrow(bg), median = median(bg$v), p = NA, bg_percentile = 50) %>%
    mutate(feature = nm, .before = 1)
  summ_all[[nm]] <- summ
  fmt <- function(m) if (f$unit == "nt") sprintf("%s nt", scales::comma(round(m))) else sprintf("%.1f%%", m)
  lab <- function(g) { s <- summ[summ$group == g, ]
    sprintf("%s (n = %s, median %s)", g, scales::comma(s$n), fmt(s$median)) }
  bg_lab <- lab("Everything else")
  med_x <- function(m) if (f$log) log2(m + 1) else m
  # Median value printed at the top of each dashed line; neighbours (sorted by x) alternate between
  # rows so close medians do not overprint. Headroom added on y so labels clear the curves.
  med_text <- function(s) {
    s <- s %>% arrange(median) %>% mutate(row = (row_number() - 1) %% nrow(.))
    # White-backed label so the dashed line does not strike through the number; the background
    # grey is darkened for the text only, since grey70 text is unreadable on white.
    s <- s %>% mutate(tcol = ifelse(group == "Everything else", "grey40", NA))
    list(geom_label(data = s, aes(x = med_x(median), y = Inf, label = map_chr(median, fmt), colour = g),
                    vjust = 1.2 + 1.5 * s$row, size = 3.3, fill = "white", label.size = 0,
                    label.padding = unit(0.08, "lines"), show.legend = FALSE),
         geom_label(data = s %>% filter(!is.na(tcol)), aes(x = med_x(median), y = Inf, label = map_chr(median, fmt)),
                    colour = "grey40", vjust = 1.2 + 1.5 * s$row[!is.na(s$tcol)], size = 3.3, fill = "white",
                    label.size = 0, label.padding = unit(0.08, "lines"), show.legend = FALSE),
         scale_y_continuous(expand = expansion(mult = c(0, 0.1 + 0.11 * nrow(s)))))
  }

  # Panel A: all targets vs everything else
  da <- bind_rows(bg %>% mutate(g = bg_lab), d %>% filter(target) %>% mutate(g = lab("All si3d targets")))
  cols_a <- setNames(c("grey70", "#E69F00"), c(bg_lab, lab("All si3d targets")))
  sa <- summ %>% filter(group %in% c("Everything else", "All si3d targets"))
  pa <- ggplot(da, aes(xv, colour = g)) +
    geom_density(linewidth = 0.8) +
    geom_vline(data = sa %>% mutate(g = map_chr(group, lab)), aes(xintercept = med_x(median), colour = g),
               linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    scale_colour_manual(values = cols_a, name = NULL, breaks = rev(names(cols_a))) +
    med_text(sa %>% mutate(g = map_chr(group, lab))) +
    labs(x = f$x, y = "Density", title = sprintf("si3d targets vs everything else (Wilcoxon p = %.2g)",
                                                 summ$p[summ$group == "All si3d targets"])) +
    theme_classic(base_size = 12) + theme(legend.position = "bottom", legend.direction = "vertical")

  # Panel B: each cluster vs everything else
  cl_labs <- setNames(map_chr(names(cl_cols), lab), names(cl_cols))
  db <- bind_rows(bg %>% mutate(g = bg_lab),
                  d %>% filter(target) %>% mutate(g = cl_labs[paste("Cluster", cluster)]))
  cols_b <- c(setNames("grey70", bg_lab), setNames(cl_cols, cl_labs))
  sb <- summ %>% filter(group %in% c("Everything else", names(cl_cols))) %>% mutate(g = map_chr(group, lab))
  pb <- ggplot(db, aes(xv, colour = g)) +
    geom_density(linewidth = 0.8) +
    geom_vline(data = sb, aes(xintercept = med_x(median), colour = g),
               linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    scale_colour_manual(values = cols_b, name = NULL, breaks = names(cols_b)) +
    med_text(sb) +
    labs(x = f$x, y = "Density", title = "Each target cluster vs everything else") +
    theme_classic(base_size = 12) + theme(legend.position = "bottom", legend.direction = "vertical")

  p <- pa + pb + plot_annotation(
    title = sprintf("MDA-MB-231 si3d hypoxia 1hr targets: %s", sub(" \\(.*", "", f$x)),
    subtitle = "Everything else = all other modelled genes; dashed lines = medians; clusters = raw SHAP, k = 3 (as in the 26d heatmap)")
  ggsave(here("plots", sprintf("26g_density_%s_targets_and_clusters_si3d_hypoxia_1hr.pdf", nm)), p, width = 12, height = 5.5)
}
out <- bind_rows(summ_all)
print(out %>% mutate(median = round(median, 1), bg_percentile = round(bg_percentile), p = signif(p, 2)), n = Inf)
write_csv(out, file.path(O, "26g_target_cluster_feature_medians_si3d_hypoxia_1hr.csv"))
