#!/usr/bin/env Rscript
# Mutation test for the validation checks in 10_rf_regression.Rmd.
#
# A check that has never been shown to FAIL is uninformative. In this repo that is not a
# hypothetical: the 01h_ Kozak PWM passed its own tier-monotonicity check for years while being
# wrong at position -3, and g4_region_swap passed every check in 01_ because the motif spot-check
# asserted only an ordering WITHIN a region. Range checks pass on a scrambled join.
#
# So this renders 10_ once per deliberate corruption and asserts the corresponding check fires.
# An uncaught mutation is the useful output: it names a blind spot and is reported, not hidden.
# Nothing here writes a real output - 10_'s export chunk is gated on mutate_target == "".
#
# Usage: Rscript code/translational_control_overview/10z_regression_mutation_tests.R
suppressMessages({ library(here); library(rmarkdown) })

if (Sys.getenv("RSTUDIO_PANDOC") == "")
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")

src <- here("code", "translational_control_overview")
source(file.path(src, "02_sweep_presets.R"))
source(file.path(src, "10_regression_presets.R"))

# `extra` carries params the mutation needs in order to REACH the check it targets. Two examples
# matter: ecdf_on_all_rows only has an ECDF to corrupt under rank_ecdf, and
# drop_variant_pair_guard only trips the two-start-context assert when initiation_score is
# exclusive (under "both" that assert is deliberately skipped).
MUTATIONS <- list(
  list(name = "target_is_matrix_te_lfc",
       breaks = "F7 class: the matrix's own te_lfc (siCTRL hypoxia-vs-normoxia) used as the target",
       expect = "target provenance, max|dev| vs matrix te_lfc",
       extra  = list()),
  list(name = "target_is_indiv_te",
       breaks = "indiv_te_hypoxia_1hr (an absolute TE level) used as the target",
       expect = "target provenance, correlation with an absolute TE level",
       extra  = list()),
  list(name = "flip_target_sign",
       breaks = "the target negated, so eIF3d-promoted genes look inhibited",
       expect = "geneset separation spot-check (the knockdown gene is held fixed so the contrast anchor passes)",
       extra  = list()),
  list(name = "permute_target",
       breaks = "the target permuted across genes, destroying every real association",
       expect = "geneset separation spot-check (knockdown gene held fixed)",
       extra  = list()),
  list(name = "leak_target_as_feature",
       breaks = "the target written INTO an existing modelled column, leaving names, order and count intact",
       expect = "value-based collinearity guard (name regex AND column-order cannot see this)",
       extra  = list()),
  list(name = "impute_from_all_rows",
       breaks = "imputation medians computed on train AND test",
       expect = "train-vs-all median divergence (the discriminative form)",
       extra  = list()),
  list(name = "ecdf_on_all_rows",
       breaks = "the rank map built on train AND test rows",
       expect = "discriminative ECDF check",
       extra  = list(target_transform = "rank_ecdf")),
  list(name = "drop_variant_pair_guard",
       breaks = "both members of every variant pair admitted to one forest",
       expect = "two start-context scores / feature parity",
       extra  = list(initiation_score = "pwm")),
  # Not a mutate_target hook: the corruption is the PARAM, pointing at a real file for the wrong
  # contrast. The filename sanity check is skipped when outcome_csv is given explicitly, which is
  # exactly the gap the golden EIF3D snapshot is there to close.
  list(name = "wrong_condition_file",
       breaks = "a real outcome file for the WRONG contrast (normoxia 1hr) passed explicitly",
       expect = "contrast anchor, parsed from the params-implied file",
       extra  = list(outcome_csv = "translation_categories_si3d_vs_sictrl_normoxia_1hr.csv"))
)

rmd     <- file.path(src, "10_rf_regression.Rmd")
scratch <- file.path(tempdir(), "tcoreg_mutation"); dir.create(scratch, showWarnings = FALSE)

run_one <- function(mut, extra) {
  err <- NULL
  p <- utils::modifyList(
    list(mutate_target = mut, run_cv = FALSE, save_plots = FALSE,
         save_model_data = FALSE, auc_bridge = FALSE, gate_feature_parity = TRUE,
         feature_set_label = "tcoreg_mut"),
    extra)
  # wrong_condition_file is a param corruption, not a hook, so it must not also set the hook.
  if (mut == "wrong_condition_file") p$mutate_target <- ""
  # Each render gets its OWN intermediates_dir. knitr writes <name>.knit.md next to the .Rmd by
  # default, so a concurrent render of the same notebook - a sweep driver running in another
  # shell - deletes this one's intermediate and the render dies in pandoc. That failure would be
  # recorded as CAUGHT while the targeted check never ran, which is the same way the ungated
  # auc_bridge chunks fooled this suite.
  inter <- file.path(scratch, paste0("inter_", mut)); dir.create(inter, showWarnings = FALSE)
  tryCatch(
    render(rmd, params = p, output_file = paste0("mut_", mut, ".html"),
           output_dir = scratch, intermediates_dir = inter, quiet = TRUE),
    error = function(e) err <<- conditionMessage(e))
  err
}

# `caught <- !is.null(err)` alone cannot tell a check from a crash: every mutation renders with
# auc_bridge = FALSE, so any chunk that reads an auc_bridge-only object would abort the render and
# be scored CAUGHT while the targeted check never ran. That is exactly what happened before
# fig_fit and fig_roc_and_calibration were gated on the param. Fixing the cause is not enough -
# the suite must also be unable to be fooled the same way again, so it now records WHICH error
# fired and flags anything that reads as infrastructure rather than a deliberate assert.
INFRA_PATTERNS <- c("object 'br' not found", "object 'bridge' not found",
                    # A clobbered knitr intermediate, i.e. a concurrent render of this notebook.
                    "knit.md", "pandoc document conversion failed")
looks_infrastructural <- function(msg) {
  # Guarded against length != 1 rather than assuming conditionMessage() returns a scalar: `if`
  # on a length-2 logical is an error in R >= 4.2, and that would fail the whole suite.
  if (!length(msg) || anyNA(msg)) return(FALSE)
  any(vapply(INFRA_PATTERNS,
             function(pat) any(grepl(pat, msg, fixed = TRUE)), logical(1)))
}

cat("Mutation test:", length(MUTATIONS), "corruptions of 10_rf_regression.Rmd\n")
cat("A caught mutation means the check works. An uncaught one is a blind spot.\n\n")

results <- vector("list", length(MUTATIONS)); t0 <- Sys.time()
for (i in seq_along(MUTATIONS)) {
  m <- MUTATIONS[[i]]
  cat(sprintf("[%2d/%2d] %-26s ... ", i, length(MUTATIONS), m$name))
  err <- run_one(m$name, m$extra)
  caught <- !is.null(err)
  infra  <- looks_infrastructural(err)
  fired_expected <- caught && !infra
  cat(if (!caught) "*** NOT CAUGHT ***\n"
      else if (infra) "*** CRASH, NOT A CHECK ***\n" else "CAUGHT\n")
  if (caught) cat("         ", sub("^.*Error *: *", "", trimws(err)), "\n")
  results[[i]] <- data.frame(mutation = m$name, breaks = m$breaks, expected_check = m$expect,
                             caught = caught, fired_expected = fired_expected,
                             message = if (caught) trimws(err) else NA_character_,
                             stringsAsFactors = FALSE)
}
res <- do.call(rbind, results)

# Unit checks for the namespace guard. Deliberately NOT render-based: 10_'s export chunk is
# gated off during a mutation run, so a render could never exercise it. Testing the guard
# directly is both faster and the only way to actually reach it.
cat("\nNamespace guard (direct, not via a render):\n")
guard_cases <- list(
  list(path = "output/translational_control_overview/rf_importance_x.csv",
       must_fail = TRUE,  what = "a classification importance filename"),
  list(path = "output/translational_control_overview/shap_family_x.csv",
       must_fail = TRUE,  what = "a classification SHAP filename"),
  list(path = "output/predictive_modeling/rfreg_importance_x.csv",
       must_fail = TRUE,  what = "the frozen pipeline's directory"),
  list(path = "output/translational_control_overview/rfreg_importance_x.csv",
       must_fail = FALSE, what = "a legitimate regression filename"),
  list(path = "output/translational_control_overview/shapreg_family_x.csv",
       must_fail = FALSE, what = "a legitimate regression SHAP filename")
)
guard_rows <- lapply(guard_cases, function(g) {
  blocked <- inherits(try(assert_no_collision(g$path), silent = TRUE), "try-error")
  ok <- blocked == g$must_fail
  # Print the DIRECTORY too: two of these cases share a basename and differ only in the
  # directory, so printing basename alone made the table look like it repeated a row.
  cat(sprintf("  %-58s %-9s %s\n", g$path,
              if (blocked) "BLOCKED" else "allowed", if (ok) "ok" else "*** WRONG ***"))
  # These call the guard directly rather than through a render, so there is no render error to
  # misattribute: fired_expected is the same fact as caught here.
  data.frame(mutation = paste0("namespace:", basename(g$path)), breaks = g$what,
             expected_check = "assert_no_collision", caught = ok, fired_expected = ok,
             message = NA_character_, stringsAsFactors = FALSE)
})
res <- rbind(res, do.call(rbind, guard_rows))

cat("\n", strrep("-", 74), "\n", sep = "")
cat(sprintf("Caught %d of %d in %.1f min (%d by the expected check)\n",
            sum(res$caught), nrow(res),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            sum(res$fired_expected)))
infra_rows <- which(res$caught & !res$fired_expected)
if (length(infra_rows)) {
  cat("\n*** THE SUITE IS NOT MEASURING WHAT IT CLAIMS ***\n")
  cat("These mutations aborted the render before their check could run, so CAUGHT here says\n")
  cat("nothing about the check. Gate the offending chunk on its param, then re-run:\n")
  for (j in infra_rows)
    cat(sprintf("  %-26s %s\n", res$mutation[j], res$message[j]))
}
if (any(!res$caught)) {
  cat("\nBLIND SPOTS - these corruptions passed every check:\n")
  for (j in which(!res$caught))
    cat(sprintf("  %-26s %s\n", res$mutation[j], res$breaks[j]))
  cat("\nThat is a finding about the checks, not a failure of this script.\n")
}

# A mutation run must never leave an output behind.
stray <- list.files(here("output", "translational_control_overview"), pattern = "tcoreg_mut")
if (length(stray)) {
  cat("\n*** A mutation run wrote output, which it must never do:\n")
  cat(paste("   ", stray, collapse = "\n"), "\n")
} else cat("\nNo mutation run wrote an output file, as intended.\n")

out <- here("output", "translational_control_overview", "regression_mutation_test_results.csv")
write.csv(res, out, row.names = FALSE)
cat("Wrote", out, "\n")
if (any(!res$fired_expected)) quit(status = 1)
