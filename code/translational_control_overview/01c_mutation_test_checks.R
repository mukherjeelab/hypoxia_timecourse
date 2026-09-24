# Mutation test for the validation checks in 01_build_feature_matrix.Rmd.
#
# A check that has never been shown to FAIL is uninformative. The 01h_ Kozak PWM passed its
# own tier-monotonicity check for years while being wrong at position -3, because a
# within-tier error cancels in the tier means. Range checks likewise pass on a scrambled
# join. So this renders 01_ once per deliberate corruption and asserts that the corresponding
# check fires.
#
# A mutation that is NOT caught is the useful output: it names a blind spot, and is reported
# rather than hidden. Nothing here writes the real matrix - 01_'s write chunk is gated on
# mutate_input == "".
#
# Usage: Rscript code/translational_control_overview/01c_mutation_test_checks.R
suppressMessages({library(here); library(rmarkdown)})

if (Sys.getenv("RSTUDIO_PANDOC") == "")
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")

# Each mutation names the check it is meant to trip. "expect" is documentation of intent; the
# script reports what actually happened either way.
MUTATIONS <- list(
  list(name = "kozak_minus3_bug",
       breaks = "the historical 01h_ PWM error: G at -3 set to 0.13",
       expect = "within-tier purine check"),
  list(name = "noderer_window_shift",
       breaks = "TIS context read at -5..+6 instead of -6..+5",
       expect = "Noderer coverage / -3 purine / +4 G on our transcripts"),
  list(name = "cridge_permute",
       breaks = "readthrough values permuted across the 192 contexts",
       expect = "TGA > TAG > TAA readthrough ordering"),
  list(name = "term_stop_shuffle",
       breaks = "stop-codon identity permuted across transcripts",
       expect = "TGA > TAG > TAA readthrough ordering"),
  list(name = "peptide_frame_shift",
       breaks = "CDS translated in the +1 frame",
       expect = "collagens are proline-rich"),
  list(name = "g4_region_swap",
       breaks = "5'UTR and 3'UTR G4mer summaries exchanged",
       expect = "G4 motif spot-check / residual orthogonality"),
  # One per defect class the independent audit found. Each of these WAS a real bug that
  # every check of the day passed, so each now has a mutation proving the new check fires.
  list(name = "struct_zero_fill",
       breaks = "F1 class: structure NA refilled with the 0 sentinel",
       expect = "invariant D, unexplained sentinel block"),
  list(name = "gene_level_join",
       breaks = "F3 class: transcript-keyed source joined on gene",
       expect = "invariant E, occupancy vs a fresh transcript-level join"),
  list(name = "duplicate_feature",
       breaks = "F4 class: two modelled columns made bit-identical",
       expect = "invariant A, identical modelled columns"),
  list(name = "constant_feature",
       breaks = "F5 class: a modelled feature made constant",
       expect = "invariant B, constant modelled columns"),
  list(name = "length_proxy_feature",
       breaks = "F2 class: a non-length feature set to a region length",
       expect = "invariant C, length proxy outside the length family"),
  list(name = "golden_drift",
       breaks = "every cds_gc shifted by 1e-6 - too small for any range or family check",
       expect = "golden-transcript snapshot")
)

rmd     <- here("code", "translational_control_overview", "01_build_feature_matrix.Rmd")
scratch <- file.path(tempdir(), "tco_mutation"); dir.create(scratch, showWarnings = FALSE)

run_one <- function(mut) {
  err <- NULL
  tryCatch(
    render(rmd, params = list(mutate_input = mut),
           output_file = paste0("mut_", mut, ".html"), output_dir = scratch, quiet = TRUE),
    error = function(e) err <<- conditionMessage(e)
  )
  err
}

cat("Mutation test:", length(MUTATIONS), "corruptions\n")
cat("A caught mutation means the check works. An uncaught one is a blind spot.\n\n")

results <- vector("list", length(MUTATIONS))
t0 <- Sys.time()
for (i in seq_along(MUTATIONS)) {
  m <- MUTATIONS[[i]]
  cat(sprintf("[%d/%d] %-22s ... ", i, length(MUTATIONS), m$name))
  err <- run_one(m$name)
  caught <- !is.null(err)
  cat(if (caught) "CAUGHT\n" else "*** NOT CAUGHT ***\n")
  if (caught) cat("        ", sub("^.*Error *: *", "", trimws(err)), "\n")
  results[[i]] <- data.frame(mutation = m$name, breaks = m$breaks,
                             expected_check = m$expect, caught = caught,
                             message = if (caught) trimws(err) else NA_character_,
                             stringsAsFactors = FALSE)
}
res <- do.call(rbind, results)

cat("\n", strrep("-", 70), "\n", sep = "")
cat(sprintf("Caught %d of %d in %.1f min\n", sum(res$caught), nrow(res),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (any(!res$caught)) {
  cat("\nBLIND SPOTS - these corruptions produced a matrix that passed every check:\n")
  for (j in which(!res$caught))
    cat(sprintf("  %-22s %s\n", res$mutation[j], res$breaks[j]))
  cat("\nThat is a finding about the checks, not a failure of this script.\n")
}

out <- here("output", "translational_control_overview", "mutation_test_results.csv")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
write.csv(res, out, row.names = FALSE)
cat("\nWrote", out, "\n")
