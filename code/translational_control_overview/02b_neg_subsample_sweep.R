# Negative-subsample sweep for 02_rf_model.Rmd
#
# Notebook 02 keeps every positive and draws a size-matched negative class from the pool
# (1,436 of 2,789 here), so a single run reports the AUC and feature ranking of ONE arbitrary
# draw. Any current-vs-corrected difference has to be read against that spread, not against
# zero. This renders 02 once per neg_seed for both variant sets; 02c_ reports the result.
#
# neg_seed seeds ONLY the negative draw - the split, the forest and the CV stay pinned at 9 -
# so every difference across seeds is negative-class composition and nothing else.
#
# Usage:  Rscript code/translational_control_overview/02b_neg_subsample_sweep.R
suppressMessages({library(here); library(rmarkdown)})

if (Sys.getenv("RSTUDIO_PANDOC") == "")
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")

SEEDS <- 1:20

# Each config is one arm of a paired comparison. Because neg_seed seeds only the draw, two
# arms at the same seed see identical genes, so differences are the config alone.
CONFIGS <- list(
  no_peptide  = list(feature_set_label = "tco_pepnone",
                     exclude_families  = "nascent_peptide,polya_track"),
  with_peptide = list(feature_set_label = "tco_pepreal"),
  shuf_peptide = list(feature_set_label = "tco_pepshuf",
                      shuffle_families  = "nascent_peptide,polya_track")
)
rmd      <- here("code", "translational_control_overview", "02_rf_model.Rmd")
scratch  <- file.path(tempdir(), "tco_sweep"); dir.create(scratch, showWarnings = FALSE)

grid <- expand.grid(seed = SEEDS, config = names(CONFIGS), stringsAsFactors = FALSE)
cat("Rendering", nrow(grid), "models (", length(SEEDS), "seeds x",
    length(CONFIGS), "configs )\n\n")

t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  s <- grid$seed[i]; nm <- grid$config[i]; cfg <- CONFIGS[[nm]]
  cat(sprintf("[%2d/%d] config=%-17s neg_seed=%2d ... ", i, nrow(grid), nm, s))
  ok <- tryCatch({
    render(rmd,
           params = c(cfg, list(variant_set = "corrected", neg_seed = s,
                                gate_against_baseline = FALSE, save_plots = FALSE)),
           # keep the 40 throwaway HTML files out of the repo
           output_file = sprintf("sweep_%s_seed%02d.html", nm, s),
           output_dir  = scratch,
           quiet = TRUE)
    TRUE
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); FALSE })
  if (ok) cat("ok\n")
}
cat("\nElapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("Outputs in output/translational_control_overview/ with suffix _negseed<N>\n")
cat("(seed 9 writes no suffix - it is the default draw)\n")
