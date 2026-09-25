# Which experiments can support a meaningful negative-draw sweep, and which cannot.
#
# The noise floor measured on one experiment does NOT transfer to another, because it is
# driven by the POOL RATIO - how large the negative pool is relative to the positive set.
# Notebook 02 keeps every positive and draws a size-matched negative class, so:
#
#   ratio ~ 1      every seed draws almost the same genes. The floor collapses toward zero
#                  and a sweep produces N copies of one model. A sweep here is meaningless,
#                  and worse, it LOOKS reassuringly stable.
#   ratio ~ 2      draws overlap about half (the reference experiment: 2,789 pool vs 1,436
#                  positives, 52% measured overlap)
#   ratio >> 2     draws barely overlap, the floor is wide, and single-run rankings are
#                  correspondingly less trustworthy
#
# Expected overlap between two independent size-matched draws from a pool of size P taking
# n each is n/P, i.e. 1/ratio. So the ratio predicts the floor before any model is fit.
#
# This costs no model fitting: it reads the gene-set CSVs only.
#
# Usage: Rscript code/translational_control_overview/03e_experiment_pool_ratios.R
suppressMessages({library(tidyverse); library(here)})

gs_dir <- here("output", "genesets")
fm <- readRDS(here("output", "translational_control_overview", "feature_matrix_tco.rds"))
in_matrix <- unique(fm$gene_id_clean)

read_set <- function(f) {
  d <- suppressWarnings(read_csv(f, show_col_types = FALSE, progress = FALSE))
  col <- intersect(c("ensembl_gene", "gene_id_clean", "gene_id"), colnames(d))[1]
  if (is.na(col)) return(character(0))
  unique(sub("\\..*", "", d[[col]]))
}

# Positive sets follow {condition}_{gene}_{direction}_{readout}_{timepoint}_lfc{lfc}.csv
pos_files <- list.files(gs_dir, pattern = "_(TE|RNA)_.*lfc[0-9.]+\\.csv$", full.names = TRUE)
pos_files <- pos_files[!str_detect(basename(pos_files), "negative")]

neg_candidates <- c(
  here("output", "predictive_modeling", "negative_control_genes.csv"),
  list.files(gs_dir, pattern = "negative", full.names = TRUE)
)
negs <- set_names(lapply(neg_candidates, read_set), basename(neg_candidates))
negs <- negs[lengths(negs) > 0]
cat("Negative pools available:\n")
for (n in names(negs))
  cat(sprintf("  %-52s %5d genes (%d in matrix)\n", n, length(negs[[n]]),
              length(intersect(negs[[n]], in_matrix))))

default_neg <- intersect(negs[[1]], in_matrix)

res <- map_dfr(pos_files, function(f) {
  p <- intersect(read_set(f), in_matrix)
  if (!length(p)) return(tibble())
  pool <- setdiff(default_neg, p)
  tibble(geneset = basename(f), n_pos = length(p), n_pool = length(pool),
         ratio = length(pool) / length(p))
}) %>%
  filter(n_pos >= 50) %>%
  mutate(expected_overlap_pct = round(100 * pmin(1, 1 / ratio), 0),
         sweep_verdict = case_when(
           ratio < 1.1 ~ "POOL TOO SMALL - every seed draws ~the same genes",
           ratio < 1.5 ~ "marginal - draws overlap heavily, floor understated",
           ratio < 4   ~ "usable - comparable to the reference experiment",
           TRUE        ~ "wide pool - floor will be larger than the reference")) %>%
  arrange(ratio)

cat("\n=== experiments by pool ratio (against the default si3d negative pool) ===\n")
print(as.data.frame(res %>% mutate(ratio = round(ratio, 2))), right = FALSE)

cat("\n=== summary ===\n")
print(as.data.frame(res %>% count(sweep_verdict, name = "n_experiments")))

f <- here("output", "translational_control_overview", "experiment_pool_ratios.csv")
write_csv(res, f)
cat("\nWrote", f, "-", nrow(res), "experiments\n")
cat("\nNOTE: ratios use the DEFAULT si3d negative pool. Experiments that should use a\n")
cat("different pool (si3e, RNA readout, MCF7-SIX1) need their own neg_geneset_csv, and\n")
cat("their ratio recomputed against it.\n")
