#!/usr/bin/env Rscript
# Renders 10_rf_regression.Rmd once per (arm, split_seed) for one preset.
#
#   Rscript 10b_regression_split_sweep.R <preset>
#   Rscript 10b_regression_split_sweep.R            # lists the presets and exits
#
# Counterpart to 02b_neg_subsample_sweep.R, with the axis changed. Notebook 02_ keeps every
# positive and draws a size-matched negative class, so neg_seed drives its noise. A regression on
# all 10,172 genes has NO draw, so split_seed is the axis: it seeds ONLY createDataPartition,
# while the forest stays pinned at forest_seed. Within a split_seed every arm therefore sees
# identical rows, and paired differences isolate the arm alone - the paired-difference sd runs
# several-fold below the seed-to-seed sd, which is the whole reason to pair.
#
# CV, plots and model data are off: CV is five extra forests per render and dominates the cost,
# and a sweep needs neither figures nor 60 copies of an imputed split. ~32 s per render.

suppressMessages({ library(here); library(purrr) })

if (Sys.getenv("RSTUDIO_PANDOC") == "")
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")

src <- here("code", "translational_control_overview")
source(file.path(src, "02_sweep_presets.R"))       # family_arms(), sweep_suffix()
source(file.path(src, "10_regression_presets.R"))  # REG_PRESETS, REG_SPLIT_SEEDS, reg_suffix()

args   <- commandArgs(trailingOnly = TRUE)
preset <- if (length(args) >= 1) args[1] else ""

if (!nzchar(preset) || !preset %in% names(REG_PRESETS)) {
  cat("Usage: Rscript 10b_regression_split_sweep.R <preset>\n\nPresets:\n")
  for (nm in names(REG_PRESETS))
    cat(sprintf("  %-12s %d arms: %s\n", nm, length(REG_PRESETS[[nm]]),
                paste(names(REG_PRESETS[[nm]]), collapse = ", ")))
  cat(sprintf("\n%d split seeds each. ~32 s per render.\n", length(REG_SPLIT_SEEDS)))
  quit(status = if (nzchar(preset)) 1 else 0)
}

# A second argument overrides the seed list, for a smoke test: ... <preset> 1,2,3
seeds <- if (length(args) >= 2)
  as.integer(strsplit(args[2], ",")[[1]]) else REG_SPLIT_SEEDS

arms <- REG_PRESETS[[preset]]
# Two arms that resolve to the SAME suffix would overwrite each other's outputs - the bug
# 02_sweep_presets.R's header records, where a sweep silently kept only its last arm. The
# invariant is distinct suffixes, NOT distinct feature_set_labels: the transform arms
# deliberately share a label because target_transform and winsorize_y already appear in the
# suffix, whereas the family arms (none/real/shuf) are distinguished by nothing else and so
# must differ by label. Asserting on the suffix covers both cases.
suffixes <- vapply(arms, function(a) reg_suffix(a, seeds[1]), character(1))
if (anyDuplicated(suffixes))
  stop("arms resolve to the same output suffix and would overwrite each other: ",
       paste(unique(suffixes[duplicated(suffixes)]), collapse = ", "))

total <- length(arms) * length(seeds)
cat(sprintf("Preset '%s': %d arms x %d seeds = %d renders (~%.0f min)\n",
            preset, length(arms), length(seeds), total, total * 32 / 60))

t0 <- Sys.time(); i <- 0L; failures <- character(0)
for (arm_name in names(arms)) {
  arm <- arms[[arm_name]]
  for (s in seeds) {
    i <- i + 1L
    suf <- reg_suffix(arm, s)
    ok <- tryCatch({
      rmarkdown::render(
        file.path(src, "10_rf_regression.Rmd"),
        params = utils::modifyList(arm, list(
          split_seed = s, run_cv = FALSE, save_plots = FALSE,
          save_model_data = FALSE, gate_feature_parity = FALSE)),
        # Every render must go to its own temp file: a shared output_file makes concurrent or
        # interrupted runs clobber each other's HTML while the CSVs look fine.
        output_file = file.path(tempdir(), paste0("10_", suf, ".html")),
        quiet = TRUE)
      TRUE
    }, error = function(e) { message("FAILED ", arm_name, " seed ", s, ": ",
                                     conditionMessage(e)); FALSE })
    if (!ok) failures <- c(failures, paste0(arm_name, "/", s))
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%3d/%3d] %-7s seed %2d  %5.1f min elapsed, ~%4.1f min left\n",
                i, total, arm_name, s, el, el / i * (total - i)))
  }
}

cat(sprintf("\nDone in %.1f min. %d/%d succeeded.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), total - length(failures), total))
if (length(failures)) {
  cat("FAILURES:", paste(failures, collapse = ", "), "\n")
  quit(status = 1)
}
cat("Report with: rmarkdown::render('10c_regression_sweep_stability.Rmd', params=list(preset='",
    preset, "'))\n", sep = "")
