#!/usr/bin/env Rscript
# 33h_construct_score_sweep.R - score the five reporter constructs across the replication axis of
# each intrinsic model, so every number is read against its own noise rather than taken from one
# arbitrary fit.
#
# WHY A SWEEP. Notebook 02_ keeps every positive and draws a size-matched negative class, so one
# run reports the probability of ONE draw; on these models the draw alone moves test AUC by ~0.07.
# A construct's P(target) inherits that. The two model types have DIFFERENT noise sources and so
# different axes:
#
#   classifier   neg_seed    1-20   seeds only the negative draw; split, forest and CV stay at 9
#   regression   split_seed  1-20   all genes are used, so there is NO draw - the split is the axis
#
# Both are therefore paired across seeds: within a seed every construct is scored by a model that
# saw identical genes.
#
# PER SEED: render the notebook with save_model_data = TRUE, load the fitted forest, predict the
# five constructs, compute SHAP for the regression, then DELETE that seed's model files. Disk stays
# flat instead of holding 80 forests.
#
# SHAP IS REGRESSION-ONLY, deliberately. The classifier is a ranger PROBABILITY forest, which
# treeshap cannot unify; 03_shap.Rmd handles that with a surrogate regression forest and its own
# fidelity bound, which is not worth importing for five sequences. So the classifier answers "is
# it called a target" and the regression answers "what drives it", on the predicted-TE-log2FC
# scale. This mirrors what 33_plasmid_predictions.Rmd already does.
#
# Usage:
#   Rscript code/translational_control_overview/33h_construct_score_sweep.R
#   Rscript code/translational_control_overview/33h_construct_score_sweep.R 1,2    # smoke test

suppressPackageStartupMessages({
  library(tidyverse); library(here); library(rmarkdown); library(ranger); library(treeshap)
})

if (Sys.getenv("RSTUDIO_PANDOC") == "")
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")

SRC <- here("code", "translational_control_overview")
O   <- here("output", "translational_control_overview")
source(file.path(SRC, "02_sweep_presets.R"))       # sweep_suffix()
source(file.path(SRC, "10_regression_presets.R"))  # reg_suffix()

args  <- commandArgs(trailingOnly = TRUE)
SEEDS <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else 1:20
stopifnot("no seeds to run" = length(SEEDS) > 0, "seeds must be positive" = all(SEEDS > 0))

# --- the two conditions, and the shared model configuration --------------------------------
CONDS <- tribble(
  ~label,            ~condition,  ~timepoint,
  "hypoxia 1hr",     "hypoxia",   "1hr",
  "normoxia 1+4hr",  "normoxia",  "1and4hr")

BASE <- list(gene = "3d", direction = "promotes", readout = "TE", lfc = "0.5",
             include_blocks = "intrinsic", variant_set = "corrected",
             initiation_score = "both", clip_mode = "aggregate",
             matrix_rds = "feature_matrix_tco.rds")
CLF_LABEL <- "tco_intrinsic"
REG_LABEL <- "tcoreg_intrinsic"

clf_sfx <- function(cond, tp, seed)
  sweep_suffix(c(BASE[c("gene","direction","readout","lfc","variant_set",
                        "initiation_score","clip_mode")],
                 list(condition = cond, timepoint = tp, feature_set_label = CLF_LABEL)), seed)

reg_sfx <- function(cond, tp, seed)
  reg_suffix(c(BASE[c("gene","direction","readout","lfc","variant_set",
                      "initiation_score","clip_mode")],
               list(condition = cond, timepoint = tp, feature_set_label = REG_LABEL,
                    target_transform = "raw")), seed)

# --- the constructs, and the feature list the SAVED models use -----------------------------
f <- readRDS(file.path(O, "plasmid_features.rds"))
stopifnot("expected five constructs - run 33a_ then recompute the features" = nrow(f) == 5)

ref_sfx <- clf_sfx(CONDS$condition[1], CONDS$timepoint[1], 9L)
ref_md  <- file.path(O, paste0("rf_model_data_", ref_sfx, ".rds"))
stopifnot("the saved seed-9 intrinsic classifier is missing - the params no longer reproduce it" =
            file.exists(ref_md))
FC <- readRDS(ref_md)$feature_cols
cat(sprintf("Reference model: %s\n  %d features\n", ref_sfx, length(FC)))
stopifnot("a model feature is missing for the constructs" = all(FC %in% names(f)),
          "a construct feature is NA - the models would see an imputed value" = !anyNA(f[, FC]))

X <- as.data.frame(f[, FC]); rownames(X) <- f$id

ref_matrix <- readRDS(file.path(O, "feature_matrix_tco.rds"))
REF <- ref_matrix %>% dplyr::select(transcript_id_clean, all_of(FC))

#' Predict the five constructs, the held-out split and the whole transcriptome from one forest.
#' The reference matrix is imputed with the TRAINING medians, exactly as the model was fitted;
#' imputing with a median leaves it unchanged, so the saved (already imputed) training split
#' carries the same medians.
score_one <- function(mdl, md, kind) {
  meds <- md$train %>% dplyr::select(all_of(FC)) %>%
    summarise(across(everything(), ~ median(.x, na.rm = TRUE)))
  ref <- REF
  for (cl in FC) ref[[cl]][is.na(ref[[cl]])] <- meds[[cl]]

  pr <- function(nd) {
    p <- predict(mdl, data = as.data.frame(nd))$predictions
    if (is.matrix(p)) p[, md$pos_label] else p
  }
  ref_pred  <- pr(ref[, FC])
  test_pred <- pr(md$test[, FC])
  is_pos <- if (kind == "classifier") md$test$eif3d_target == md$pos_label else md$test$y > 0.5
  cpred  <- pr(X[, FC])

  tibble(construct = rownames(X), pred = cpred,
         pct_transcriptome = map_dbl(cpred, ~ 100 * mean(ref_pred  < .x)),
         pct_heldout_pos   = map_dbl(cpred, ~ 100 * mean(test_pred[is_pos] < .x)),
         pct_heldout_neg   = map_dbl(cpred, ~ 100 * mean(test_pred[!is_pos] < .x)))
}

#' Exact per-feature SHAP for the five constructs. Regression forests only - treeshap cannot
#' unify a probability forest.
shap_one <- function(mdl, md) {
  uni <- ranger.unify(mdl, as.data.frame(md$train[, FC]))
  s   <- treeshap(uni, X[, FC], verbose = FALSE)$shaps
  as_tibble(s) %>% mutate(construct = rownames(X), .before = 1) %>%
    pivot_longer(-construct, names_to = "feature", values_to = "shap")
}

scratch <- file.path(tempdir(), "tco_constructs"); dir.create(scratch, showWarnings = FALSE)
preds <- list(); shaps <- list(); failed <- character(0)
t0 <- Sys.time()
n_tot <- nrow(CONDS) * length(SEEDS) * 2L
i <- 0L

for (k in seq_len(nrow(CONDS))) {
  cond <- CONDS$condition[k]; tp <- CONDS$timepoint[k]; lab <- CONDS$label[k]

  # The regression gates its feature list against the classifier's. The default builds that
  # suffix with feature_set_label = "tco", which is the full model; an intrinsic run has to be
  # pointed at the intrinsic classifier or the gate compares the wrong two lists.
  clf_ref <- clf_sfx(cond, tp, 9L)

  for (s in SEEDS) {
    for (kind in c("classifier", "regression")) {
      i <- i + 1L
      sfx <- if (kind == "classifier") clf_sfx(cond, tp, s) else reg_sfx(cond, tp, s)
      cat(sprintf("[%2d/%d] %-15s %-11s seed %2d ... ", i, n_tot, lab, kind, s))

      rmd <- file.path(SRC, if (kind == "classifier") "02_rf_model.Rmd" else "10_rf_regression.Rmd")
      p <- c(BASE, list(condition = cond, timepoint = tp,
                        save_plots = FALSE, save_model_data = TRUE))
      p <- if (kind == "classifier")
        c(p, list(feature_set_label = CLF_LABEL, neg_seed = s, gate_against_baseline = FALSE))
      else
        c(p, list(feature_set_label = REG_LABEL, target_transform = "raw", split_seed = s,
                  run_cv = FALSE, permutation_check = FALSE, clf_suffix_override = clf_ref))

      ok <- tryCatch({
        render(rmd, params = p, quiet = TRUE,
               output_file = sprintf("cs_%s_%s_%02d.html", kind, cond, s),
               output_dir  = scratch)
        TRUE
      }, error = function(e) { cat("RENDER FAILED:", conditionMessage(e), "\n"); FALSE })
      if (!ok) { failed <- c(failed, sfx); next }

      stem <- if (kind == "classifier") c("rf_classifier_", "rf_model_data_")
              else                      c("rfreg_model_",   "rfreg_model_data_")
      fm <- file.path(O, paste0(stem[1], sfx, ".rds"))
      fd <- file.path(O, paste0(stem[2], sfx, ".rds"))
      if (!file.exists(fm) || !file.exists(fd)) {
        cat("MODEL NOT FOUND:", basename(fm), "\n"); failed <- c(failed, sfx); next
      }

      mdl <- readRDS(fm); md <- readRDS(fd)
      stopifnot("this seed's model uses a different feature list" =
                  identical(md$feature_cols, FC))

      preds[[length(preds) + 1L]] <- score_one(mdl, md, kind) %>%
        mutate(label = lab, kind = kind, seed = s, .before = 1)
      if (kind == "regression")
        shaps[[length(shaps) + 1L]] <- shap_one(mdl, md) %>%
          mutate(label = lab, seed = s, .before = 1)

      # Delete this seed's forests. Seed 9 is the SAVED default run that 33_ reads, so it is
      # never removed - only the sweep's own throwaway fits are.
      if (s != 9L) unlink(c(fm, fd))
      cat("ok\n")
    }
  }
}

preds <- bind_rows(preds); shaps <- bind_rows(shaps)
cat(sprintf("\nElapsed: %.1f min | %d predictions | %d SHAP rows | %d failures\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            nrow(preds), nrow(shaps), length(failed)))
if (length(failed)) { cat("FAILED:\n"); cat(paste0("  ", failed, collapse = "\n"), "\n") }
stopifnot("no predictions were produced" = nrow(preds) > 0)

write_csv(preds, file.path(O, "construct_sweep_predictions.csv"))
write_csv(shaps, file.path(O, "construct_sweep_shap.csv"))
cat("wrote construct_sweep_predictions.csv, construct_sweep_shap.csv\n")

# --- what the sweep says, in one table -----------------------------------------------------
cat("\nP(target) / predicted TE log2FC across seeds:\n")
preds %>%
  group_by(kind, label, construct) %>%
  summarise(n = n(), median = median(pred), sd = sd(pred),
            p05 = quantile(pred, 0.05), p95 = quantile(pred, 0.95),
            pct_tx = median(pct_transcriptome), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print()

# For the classifier, how often does the CALL flip across draws? 0.5 is the nominal cut on a
# training set balanced 50/50, so it is the boundary the draw can move a construct across.
cat("\nClassifier call stability (P > 0.5):\n")
preds %>% filter(kind == "classifier") %>%
  group_by(label, construct) %>%
  summarise(n = n(), n_called_target = sum(pred > 0.5),
            frac = round(mean(pred > 0.5), 2), .groups = "drop") %>%
  as.data.frame() %>% print()
