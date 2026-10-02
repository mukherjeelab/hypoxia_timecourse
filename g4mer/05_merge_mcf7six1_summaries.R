# 05_merge_mcf7six1_summaries.R
# Build g4mer_summary_{region}_mcf7six1.tsv for the MCF7-SIX1 transcript set, by combining
# rows REUSED from the existing (MDA-MB-231) summaries with rows from the new inference run.
#
# Reuse is exact, not an approximation: regenerating the 5'UTR summary reproduced the stored
# values bit-identically, so the model is deterministic and a shared transcript's score does not
# depend on which run produced it.
#
# FALLBACK, deliberate: if a region's new-inference file is absent or empty, the merged summary is
# written from the reused rows ALONE. The transcripts that were never scored are then simply
# missing from the summary, so their g4mer_* columns arrive as NA in the matrix and are imputed
# like any other missing feature. That keeps the feature SET identical across cell lines, which is
# the property the heatmap depends on, at the cost of more imputation in the g4 family. The
# coverage is printed and asserted against a floor so the cost is never silent.
#
# Usage: Rscript g4mer/05_merge_mcf7six1_summaries.R [new_dir]
#   new_dir defaults to g4mer/out_mcf7six1

suppressPackageStartupMessages({ library(dplyr); library(readr); library(here) })
args    <- commandArgs(trailingOnly = TRUE)
# Optional second arg: cell = mcf7six1 (default) | hela.
CELL    <- if (length(args) >= 2) args[2] else "mcf7six1"
TX_RDS  <- c(mcf7six1 = "precomputed_mcf7six1_tx.rds", hela = "precomputed_teleman_tx.rds",
            dhx29 = "precomputed_dhx29_tx.rds", mane = "precomputed_mane_tx.rds",
           mdamb231 = "precomputed_mdamb231_tx.rds")
stopifnot("unknown cell set" = CELL %in% names(TX_RDS))
new_dir <- if (length(args) >= 1) args[1] else here("g4mer", paste0("out_", CELL))
g4_dir  <- here("output", "g4mer")
want <- unique(readRDS(here("output", "predictive_modeling", TX_RDS[[CELL]]))$transcript_id_clean)

# Below this, the modelling code drops the column outright (features under 50% non-NA), which
# would silently break feature parity with the 231 runs - the one thing that must not happen.
COVERAGE_FLOOR <- 0.50
rows <- list()

for (rg in c("utr5", "cds", "utr3")) {
  old_f <- file.path(g4_dir, paste0("g4mer_summary_", rg, ".tsv"))
  stopifnot("the MDA-MB-231 summary for this region is missing" = file.exists(old_f))
  old <- read_tsv(old_f, show_col_types = FALSE)
  # Reuse from the 231 base AND from every cell set's merged summary, INCLUDING this cell's own
  # previous output: it is read in full here, before being overwritten below. Excluding it (the
  # original rule) silently dropped everything a first run had scored whenever 01k_ ran a second
  # time for the same cell - 1,100 MANE transcripts on 2026-10-01, because 01b_ counts this file
  # as "already scored" and so never re-sends them. Scores are sequence properties, bit-identical.
  others <- list.files(g4_dir, pattern = paste0("^g4mer_summary_", rg, "_.*\\.tsv$"), full.names = TRUE)
  pool <- bind_rows(old, lapply(others, function(f) read_tsv(f, show_col_types = FALSE) %>%
                                  dplyr::select(all_of(colnames(old))))) %>%
    distinct(transcript_id_clean, .keep_all = TRUE)
  reused <- pool %>% filter(transcript_id_clean %in% want)

  new_f <- file.path(new_dir, paste0("g4mer_summary_", rg, ".tsv"))
  new <- if (file.exists(new_f)) read_tsv(new_f, show_col_types = FALSE) else old[0, ]
  new <- new %>% filter(transcript_id_clean %in% want,
                        !transcript_id_clean %in% reused$transcript_id_clean)

  stopifnot("the new inference output has different columns from the 231 summary" =
              nrow(new) == 0 || setequal(colnames(new), colnames(old)))
  merged <- bind_rows(reused, new %>% dplyr::select(all_of(colnames(old))))
  stopifnot("a transcript appears twice in the merged summary" =
              !any(duplicated(merged$transcript_id_clean)),
            "the merged summary contains a transcript outside this cell line's set" =
              all(merged$transcript_id_clean %in% want),
            "region column disagrees with the file it was written to" =
              all(merged$region == rg))

  out_f <- file.path(g4_dir, paste0("g4mer_summary_", rg, "_", CELL, ".tsv"))
  write_tsv(merged, out_f)
  cov <- nrow(merged) / length(want)
  cat(sprintf("%-5s reused %6d | new %6d | merged %6d of %6d %s transcripts (%.1f%%)%s\n",
              rg, nrow(reused), nrow(new), nrow(merged), length(want), CELL, 100 * cov,
              if (nrow(new) == 0) "   [FALLBACK: no new inference, reused only]" else ""))
  rows[[rg]] <- tibble(region = rg, reused = nrow(reused), new = nrow(new),
                       merged = nrow(merged), coverage = cov)
}

r <- bind_rows(rows); cat("\n"); print(as.data.frame(r), row.names = FALSE)
low <- r %>% filter(coverage < COVERAGE_FLOOR)
if (nrow(low)) {
  cat("\nWARNING - these regions fall below the ", COVERAGE_FLOOR * 100,
      "% floor at which the modelling code DROPS a feature:\n", sep = "")
  print(as.data.frame(low), row.names = FALSE)
  cat("Dropping would break feature parity with the MDA-MB-231 runs. Either finish the G4mer\n",
      "inference for this region, or exclude the g4 family from BOTH cell lines.\n", sep = "")
}
stopifnot("a merged G4mer summary is below the coverage floor - see the warning above" =
            all(r$coverage >= COVERAGE_FLOOR))
cat("\nwrote g4mer_summary_{utr5,cds,utr3}_", CELL, ".tsv - use g4mer_tag = '_", CELL, "'\n", sep = "")
