#!/usr/bin/env Rscript
# SHAP noise floor for the regression, over the split_seed axis.
#
#   Rscript 11d_shap_regression_seed_sweep.R <preset> [arm] [n_seeds]
#   Rscript 11d_shap_regression_seed_sweep.R transform raw 20
#
# One run per experiment is defensible only against a MEASURED floor. For the classifier that
# floor came from the negative draw (03d_/03f_); there is no draw here, so it comes from the
# train/test split. Nothing biological varies: the same 10,172 genes every time, a different 80%
# trained on, forest pinned at 9.
#
# The floor is NOT interchangeable with the classifier's. Two 80/20 splits share ~81% of their
# training rows, against the ~52% gene overlap a negative draw gave, so this floor should come
# out tighter - and quantifying that is the concrete payoff of dropping the draw.
#
# Each seed needs BOTH renders: 10_ to fit and save the model plus its imputed split, then 11_ to
# explain it. 11_ runs with sweep_mode = TRUE, which computes ONE SHAP matrix instead of three.

suppressMessages({ library(here) })
if (Sys.getenv("RSTUDIO_PANDOC") == "")
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")

src <- here("code", "translational_control_overview")
source(file.path(src, "02_sweep_presets.R"))
source(file.path(src, "10_regression_presets.R"))

args    <- commandArgs(trailingOnly = TRUE)
preset  <- if (length(args) >= 1) args[1] else "transform"
arm_nm  <- if (length(args) >= 2) args[2] else NA_character_
n_seeds <- if (length(args) >= 3) as.integer(args[3]) else length(REG_SPLIT_SEEDS)

if (!preset %in% names(REG_PRESETS))
  stop("unknown preset '", preset, "'. Options: ", paste(names(REG_PRESETS), collapse = ", "))
arms <- REG_PRESETS[[preset]]
if (!is.na(arm_nm)) {
  if (!arm_nm %in% names(arms))
    stop("arm '", arm_nm, "' not in preset '", preset, "'. Arms: ",
         paste(names(arms), collapse = ", "))
  arms <- arms[arm_nm]
}
seeds <- REG_SPLIT_SEEDS[seq_len(min(n_seeds, length(REG_SPLIT_SEEDS)))]

# An arm carries whatever params its question needs, and 10_ declares all of them. 11_ declares a
# subset - it never selects features itself, so it has no exclude_families or shuffle_families -
# and rmarkdown hard-errors on a param its YAML does not declare. Forwarding the whole arm
# therefore crashed the `none` and `shuf` arms of every family preset before they reached SHAP.
# The declared names are read from 11_'s own front matter rather than listed here, so this cannot
# drift out of date when a param is added on either side.
shap_params <- names(rmarkdown::yaml_front_matter(
  file.path(src, "11_shap_regression.Rmd"))$params)
stopifnot("could not read 11_shap_regression.Rmd's declared params" = length(shap_params) > 0)
arm_for_shap <- function(arm) arm[intersect(names(arm), shap_params)]

total <- length(arms) * length(seeds)
cat(sprintf("SHAP floor sweep: preset '%s', arms [%s], %d seeds = %d x (10_ + 11_)\n",
            preset, paste(names(arms), collapse = ", "), length(seeds), total))

t0 <- Sys.time(); i <- 0L; failures <- character(0)
for (arm_name in names(arms)) {
  arm <- arms[[arm_name]]
  for (s in seeds) {
    i <- i + 1L
    suf <- reg_suffix(arm, s)
    ok <- tryCatch({
      # save_model_data = TRUE is required: 11_ explains THAT model's saved split rather than
      # re-deriving it from params, which is the drift 07c_ has to assert its way back from.
      rmarkdown::render(
        file.path(src, "10_rf_regression.Rmd"),
        params = utils::modifyList(arm, list(
          split_seed = s, run_cv = FALSE, save_plots = FALSE,
          save_model_data = TRUE, auc_bridge = TRUE, gate_feature_parity = FALSE)),
        output_file = file.path(tempdir(), paste0("10_", suf, ".html")), quiet = TRUE)
      rmarkdown::render(
        file.path(src, "11_shap_regression.Rmd"),
        params = utils::modifyList(arm_for_shap(arm), list(
          split_seed = s, sweep_mode = TRUE, save_shap_plots = FALSE)),
        output_file = file.path(tempdir(), paste0("11_", suf, ".html")), quiet = TRUE)
      TRUE
    }, error = function(e) { message("FAILED ", arm_name, " seed ", s, ": ",
                                     conditionMessage(e)); FALSE })
    if (!ok) failures <- c(failures, paste0(arm_name, "/", s))
    # The model RDS is ~63 MB, so 20 seeds would leave over a gigabyte behind for no reason:
    # the floor is computed from the shapreg_family_ CSVs, which are kept.
    for (f in c(paste0("rfreg_model_", suf, ".rds"), paste0("rfreg_model_data_", suf, ".rds")))
      unlink(file.path(here("output", "translational_control_overview"), f))
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%3d/%3d] %-8s seed %2d  %6.1f min elapsed, ~%6.1f min left\n",
                i, total, arm_name, s, el, el / i * (total - i)))
  }
}
cat(sprintf("\nDone in %.1f min. %d/%d succeeded.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), total - length(failures), total))
if (length(failures)) { cat("FAILURES:", paste(failures, collapse = ", "), "\n"); quit(status = 1) }
cat("Now run 11f_shap_regression_stability.Rmd, then 10c_ and 12_.\n")
