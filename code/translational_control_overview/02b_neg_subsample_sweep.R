# Negative-subsample sweep for 02_rf_model.Rmd
#
# Notebook 02 keeps every positive and draws a size-matched negative class from the pool, so a
# single run reports the AUC and feature ranking of ONE arbitrary draw. On this model the draw
# alone moves test AUC by ~0.07. No claim is readable against zero; it has to be read against
# an arm that differs only in the thing being tested.
#
# neg_seed seeds ONLY the negative draw - the split, the forest and the CV stay pinned at 9 -
# so within a seed every arm sees identical genes in an identical split. The sweep is paired
# by construction, and the paired-difference sd runs ~5-10x smaller than the draw-to-draw sd.
#
# Arms come from 02_sweep_presets.R, shared with 02c_sweep_stability.Rmd which reports them.
#
# Usage:
#   Rscript code/translational_control_overview/02b_neg_subsample_sweep.R g4
#   Rscript code/translational_control_overview/02b_neg_subsample_sweep.R      # lists presets
suppressMessages({library(here); library(rmarkdown)})
source(here("code", "translational_control_overview", "02_sweep_presets.R"))

if (Sys.getenv("RSTUDIO_PANDOC") == "")
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")

args   <- commandArgs(trailingOnly = TRUE)
preset <- if (length(args)) args[1] else NA_character_

if (is.na(preset) || !preset %in% names(SWEEP_PRESETS)) {
  cat("Usage: Rscript 02b_neg_subsample_sweep.R <preset>\n\nAvailable presets:\n")
  for (nm in names(SWEEP_PRESETS))
    cat(sprintf("  %-12s %d arms: %s\n", nm, length(SWEEP_PRESETS[[nm]]),
                paste(names(SWEEP_PRESETS[[nm]]), collapse = ", ")))
  if (!is.na(preset)) cat("\nUnknown preset:", preset, "\n")
  quit(status = if (is.na(preset)) 0 else 1)
}

CONFIGS <- SWEEP_PRESETS[[preset]]
rmd     <- here("code", "translational_control_overview", "02_rf_model.Rmd")
scratch <- file.path(tempdir(), "tco_sweep"); dir.create(scratch, showWarnings = FALSE)

grid <- expand.grid(seed = SWEEP_SEEDS, config = names(CONFIGS), stringsAsFactors = FALSE)
cat("Preset:", preset, "-", nrow(grid), "models (", length(SWEEP_SEEDS), "seeds x",
    length(CONFIGS), "arms )\n\n")

t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  s <- grid$seed[i]; nm <- grid$config[i]; cfg <- CONFIGS[[nm]]
  cat(sprintf("[%2d/%d] arm=%-12s neg_seed=%2d ... ", i, nrow(grid), nm, s))
  ok <- tryCatch({
    render(rmd,
           params = c(cfg, list(neg_seed = s, gate_against_baseline = FALSE,
                                save_plots = FALSE, save_model_data = FALSE)),
           # keep the throwaway HTML out of the repo
           output_file = sprintf("sweep_%s_%s_seed%02d.html", preset, nm, s),
           output_dir  = scratch,
           quiet = TRUE)
    TRUE
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); FALSE })
  if (ok) cat("ok\n")
}
cat("\nElapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("Report it with 02c_sweep_stability.Rmd, params: preset =", shQuote(preset), "\n")
cat("(seed 9 writes no suffix - it is the default draw)\n")
