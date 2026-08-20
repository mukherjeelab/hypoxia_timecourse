#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 90: replicate negative-subsample driver
#
# Renders notebook 07 once per negative-draw seed, holding EVERYTHING else fixed.
# Only `neg_seed` varies: the train/test split, the forest and the CV are pinned at
# seed 9 inside 07, so any spread across replicates is attributable to which background
# genes were drawn and to nothing else.
#
# Pool:   output/genesets/si3e_expressed_background_pool.csv  (notebook 02c)
#         = expressed (DESeq2-tested) genes minus ALL 3e-promoted genes, 10,778 genes,
#           9,183 of them in the model matrix. NOT screened on effect size.
# Draw:   399 per seed, size-matched to the eIF3e-specific positive set.
#
# neg_tag = "exprbg" is REQUIRED: neg_geneset_csv is not part of 07's output suffix, so
# without a tag these runs would overwrite the runs built on the screened
# si3e_negative_controls_padj0.05.csv pool.
#
# Usage:  Rscript code/predictive_modeling/90_replicate_neg_subsample_driver.R [n_seeds]
# ---------------------------------------------------------------------------

suppressMessages(library(here))

if (!nzchar(Sys.getenv("RSTUDIO_PANDOC"))) {
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
}

args    <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1) as.integer(args[1]) else 20L
SEEDS   <- seq_len(n_seeds)

NEG_POOL <- "si3e_expressed_background_pool.csv"
NEG_TAG  <- "exprbg"

BASE_PARAMS <- list(
  condition               = "hypoxia",
  gene                    = "3e",
  direction               = "promotes",
  timepoint               = "1hr_not3d",
  lfc                     = "0.5",
  feature_set             = "dhx29_kd_feature",
  feature_matrix_path     = "feature_matrix_dhx29_kd_feature.rds",
  exclude_4e_features     = TRUE,
  exclude_own_te_features = TRUE,
  neg_geneset_csv         = NEG_POOL,
  neg_tag                 = NEG_TAG,
  save_plots              = FALSE   # 8 near-identical PDFs per replicate otherwise
)

stopifnot("negative pool not found - knit notebook 02c first" =
  file.exists(here("output", "genesets", NEG_POOL)))

rmd     <- here("code", "predictive_modeling",
                "07_rf_classification_eif3d_targets_no_boruta.Rmd")
html_dir <- here("output", "predictive_modeling", "replicates")
dir.create(html_dir, showWarnings = FALSE, recursive = TRUE)

# Rebuild 07's suffix here so the driver can tell whether a seed has already been run.
# Must track 07's construction exactly; the assertion after the first render catches drift.
suffix_for <- function(seed) {
  paste0("eif3", BASE_PARAMS$gene, "_", BASE_PARAMS$direction, "_",
         BASE_PARAMS$condition, "_", BASE_PARAMS$timepoint, "_lfc", BASE_PARAMS$lfc,
         "_", BASE_PARAMS$feature_set, "_clipaggregate_no_boruta",
         "_no4e", "_noOwnTE",
         "_", NEG_TAG,
         if (seed != 9L) paste0("_negseed", seed) else "")
}
perf_path <- function(seed) {
  here("output", "predictive_modeling",
       paste0("rf_classification_performance_", suffix_for(seed), ".csv"))
}

cat("Replicate sweep:", length(SEEDS), "seeds\n")
cat("Pool:", NEG_POOL, " tag:", NEG_TAG, "\n\n")

t0 <- Sys.time()
for (seed in SEEDS) {
  if (file.exists(perf_path(seed))) {
    cat(sprintf("[seed %2d] already run, skipping\n", seed)); next
  }
  cat(sprintf("[seed %2d] rendering ... ", seed)); flush.console()
  ts <- Sys.time()
  rmarkdown::render(
    rmd,
    params      = c(BASE_PARAMS, list(neg_seed = seed)),
    output_file = file.path(html_dir, paste0("07_", suffix_for(seed), ".html")),
    envir       = new.env(),   # each replicate gets a clean environment
    quiet       = TRUE
  )
  cat(sprintf("%.1fs\n", as.numeric(difftime(Sys.time(), ts, units = "secs"))))

  # The suffix is reconstructed above rather than read from 07. If 07's construction ever
  # changes, the file lands somewhere this driver cannot see and the sweep would silently
  # re-render the same seed forever.
  stopifnot("07's output suffix no longer matches this driver's reconstruction" =
    file.exists(perf_path(seed)))
}

cat(sprintf("\nDone in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("Analyze with notebook 91.\n")
