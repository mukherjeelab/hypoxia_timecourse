# 01n_value_parity_check.R
# Every modelled feature must be BIT-IDENTICAL between a cell set's own-isoform matrix and the
# MDA-MB-231 matrix on the transcripts both contain (Kate's rule; see the memory note
# feedback-shared-definitions-across-cell-lines). A difference means a definition was rebuilt per
# matrix, and the same column would mean different things in different heatmap columns.
#   Rscript 01n_value_parity_check.R feature_matrix_tco_hela.rds
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here) })
o <- here("output", "translational_control_overview")
f <- commandArgs(trailingOnly = TRUE)[1]
a <- readRDS(file.path(o, "feature_matrix_tco.rds")); m <- readRDS(file.path(o, f))
full <- read_csv(file.path(o, "rfreg_importance_eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_corrected_both_clipaggregate_yraw.csv"), show_col_types = FALSE)$feature
intr <- read_csv(file.path(o, "rfreg_importance_eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw.csv"), show_col_types = FALSE)$feature
sh <- intersect(a$transcript_id_clean, m$transcript_id_clean)
x <- a[match(sh, a$transcript_id_clean), ]; y <- m[match(sh, m$transcript_id_clean), ]
res <- bind_rows(lapply(full, function(cc) {
  u <- as.numeric(unlist(x[[cc]])); v <- as.numeric(unlist(y[[cc]]))
  k <- !is.na(u) & !is.na(v)
  tibble(feature = cc, intrinsic = cc %in% intr, n = sum(k), na_mismatch = sum(is.na(u) != is.na(v)),
         max_abs_dev = if (any(k)) max(abs(u[k] - v[k])) else NA_real_)
})) %>% mutate(identical = na_mismatch == 0 & (is.na(max_abs_dev) | max_abs_dev < 1e-12))
cat(sprintf("%s vs feature_matrix_tco.rds: %d shared transcripts\n", f, length(sh)))
cat(sprintf("modelled features identical: %d of %d | intrinsic: %d of %d\n",
            sum(res$identical), nrow(res), sum(res$identical & res$intrinsic), sum(res$intrinsic)))
bad <- res %>% filter(!identical)
if (nrow(bad)) print(as.data.frame(bad), row.names = FALSE)
stopifnot("an INTRINSIC feature differs from MDA-MB-231 on shared transcripts" =
            all(res$identical[res$intrinsic]))
cat("VALUE PARITY PASSED for all intrinsic features\n")
