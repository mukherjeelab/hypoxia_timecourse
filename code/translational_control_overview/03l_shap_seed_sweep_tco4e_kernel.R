# 03l_shap_seed_sweep_tco4e_kernel.R - the 20-negative-draw SHAP sweep (03d_'s design) on the
# eIF4E-inclusive classifier (all 75 features), explained with Kernel SHAP. Feeds 03i_'s interval
# bars for the model 03g_ / 03j_ plot. Experiments: hypoxia 1hr (default) and normoxia 1and4hr
# (the time-regressed classifier), chosen by `experiment` below or the 2nd argument.
#
# neg_seed seeds only the negative draw; split, forest and CV stay pinned at 9 (see 03d_).
# Seed 9 is the REFERENCE model and is deliberately NOT re-rendered: its files are what 03g_ / 03j_
# read, and 02_ / 03_ write seed 9 under the bare suffix. Every other seed writes _negseed{s}.
# Each draw is asserted to have exactly the reference model's feature set, which pins it to
# the same blocks / variant set as the model it is the interval for.
#
# Kernel SHAP is ~minutes per draw (sweep_mode skips the extra explain sets), so run it in the
# background. A seed whose SHAP output already exists and is newer than its model is skipped,
# so an interrupted sweep resumes.
# Each notebook is rendered in a FRESH environment. render() defaults to the caller's, and 02_
# assigns its own `include_blocks` there - which overwrote this driver's setting after the first
# render and silently fitted every later draw on the intrinsic block only. The feature-set
# assertion below is what caught it.
#   Rscript code/translational_control_overview/03l_shap_seed_sweep_tco4e_kernel.R [seeds] [experiment]
#   e.g. ... 1:20 normoxia_1and4hr
suppressMessages({ library(here); library(rmarkdown) })

if (Sys.getenv("RSTUDIO_PANDOC") == "")
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")

## ---- settings ----
experiment <- "hypoxia_1hr"   # or "normoxia_1and4hr"
## ------------------

args  <- commandArgs(trailingOnly = TRUE)
SEEDS <- if (length(args) >= 1) eval(parse(text = args[1])) else 1:20
if (length(args) >= 2) experiment <- args[2]
cond_tp <- switch(experiment,
                  hypoxia_1hr      = c(condition = "hypoxia",  timepoint = "1hr"),
                  normoxia_1and4hr = c(condition = "normoxia", timepoint = "1and4hr"),
                  stop("experiment must be hypoxia_1hr or normoxia_1and4hr"))
model_params <- list(condition = cond_tp[["condition"]], gene = "3d", direction = "promotes", readout = "TE",
                     timepoint = cond_tp[["timepoint"]], lfc = "0.5", clip_mode = "aggregate",
                     feature_set_label = "tco4e", variant_set = "corrected", initiation_score = "both")
include_blocks <- "intrinsic,external,eif4e"
engine         <- "kernelshap"
SEEDS <- setdiff(SEEDS, 9L)   # the reference model, never overwritten here

d <- here("code", "translational_control_overview")
O <- here("output", "translational_control_overview")
scratch <- file.path(tempdir(), "shap_sweep_tco4e"); dir.create(scratch, showWarnings = FALSE)
base <- paste0("eif33d_promotes_", experiment, "_lfc0.5_tco4e_corrected_both_clipaggregate")
ref_cols <- readRDS(file.path(O, paste0("rf_model_data_", base, ".rds")))$feature_cols

cat("SHAP seed sweep (", engine, ", ", base, "): seeds ", paste(SEEDS, collapse = ","), "\n\n", sep = "")
t0 <- Sys.time()
for (s in SEEDS) {
  sfx  <- paste0(base, "_negseed", s)
  shap_out  <- file.path(O, paste0("shap_feature_", sfx, "_", engine, ".csv"))
  model_out <- file.path(O, paste0("rf_classifier_", sfx, ".rds"))
  if (file.exists(shap_out) && file.exists(model_out) && file.mtime(shap_out) > file.mtime(model_out)) {
    cat(sprintf("seed %2d: already done, skipping\n", s)); next
  }
  cat(sprintf("seed %2d ... ", s)); t1 <- Sys.time()
  ok <- tryCatch({
    render(file.path(d, "02_rf_model.Rmd"),
           params = c(model_params, list(include_blocks = include_blocks, neg_seed = s,
                                         gate_against_baseline = FALSE, save_plots = FALSE,
                                         save_model_data = TRUE)),
           output_file = sprintf("m_%02d.html", s), output_dir = scratch, quiet = TRUE,
           envir = new.env())
    md <- readRDS(file.path(O, paste0("rf_model_data_", sfx, ".rds")))
    stopifnot("draw has a different feature set from the reference model" =
                identical(md$feature_cols, ref_cols))
    render(file.path(d, "03_shap.Rmd"),
           params = c(model_params, list(neg_seed = s, shap_engine = engine, sweep_mode = TRUE)),
           output_file = sprintf("s_%02d.html", s), output_dir = scratch, quiet = TRUE,
           envir = new.env())
    TRUE
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); FALSE })
  if (ok) cat(sprintf("ok (%.1f min)\n", as.numeric(difftime(Sys.time(), t1, units = "mins"))))
}
cat("\nElapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("Plot with 03i_signed_shap_bar_ci.R (engine = kernelshap)\n")
