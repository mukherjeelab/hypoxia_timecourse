# SHAP noise floor: how much does a family's attribution move when NOTHING biological
# changes - only the negative draw?
#
# This is the yardstick for the cross-experiment comparison. One run per experiment is the
# plan, which is defensible only if a between-experiment difference can be judged against
# the spread a single experiment shows on its own. Without that floor, "this feature matters
# more in hypoxia" cannot be distinguished from the draw.
#
# It is therefore run ONCE on a reference experiment, not per experiment.
#
# neg_seed seeds only the negative draw; the split, the forest and the CV stay pinned at 9,
# so every difference across seeds is negative-class composition and nothing else.
#
# 02_ must save its model data here (the sweep driver 02b_ turns that off to avoid 60 copies),
# because 03_ explains that run's saved split rather than re-deriving it.
#
# Usage: Rscript code/translational_control_overview/03d_shap_seed_sweep.R [n_seeds]
suppressMessages({library(here); library(rmarkdown)})

if (Sys.getenv("RSTUDIO_PANDOC") == "")
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")

args  <- commandArgs(trailingOnly = TRUE)
SEEDS <- seq_len(if (length(args)) as.integer(args[1]) else 20)

d       <- here("code", "translational_control_overview")
scratch <- file.path(tempdir(), "shap_sweep"); dir.create(scratch, showWarnings = FALSE)

cat("SHAP seed sweep:", length(SEEDS), "seeds\n\n")
t0 <- Sys.time()
for (s in SEEDS) {
  cat(sprintf("[%2d/%d] seed %2d ... ", match(s, SEEDS), length(SEEDS), s))
  ok <- tryCatch({
    render(file.path(d, "02_rf_model.Rmd"),
           params = list(neg_seed = s, gate_against_baseline = FALSE,
                         save_plots = FALSE, save_model_data = TRUE),
           output_file = sprintf("m_%02d.html", s), output_dir = scratch, quiet = TRUE)
    render(file.path(d, "03_shap.Rmd"),
           params = list(neg_seed = s, sweep_mode = TRUE),
           output_file = sprintf("s_%02d.html", s), output_dir = scratch, quiet = TRUE)
    TRUE
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); FALSE })
  if (ok) cat("ok\n")
}
cat("\nElapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("Report with 03e_shap_stability.Rmd\n")
