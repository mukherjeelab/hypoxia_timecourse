# 33d_g4_features.R - the 9 G4mer rG4 features, computed from sequence.
#
# G4mer takes at most 70 nt, so a region is scanned with a 70 nt sliding window and summarised.
# STRIDE DIFFERS BY REGION - 5'UTR 10, CDS and 3'UTR 20 - which is how the stored summaries were
# produced (verified here against their own n_windows column, 100% of rows in all three regions).
# The three regions' scores are therefore not comparable to each other and are never differenced.
#
# g4mer_max rises with region length by construction (it is a maximum over windows) and tracks GC
# because rG4 needs G-runs, so the MODELLED column is the residual of g4mer_max on
# log2(region_length + 1) + region GC. Those coefficients are LOADED from the MDA-MB-231 fit,
# never refitted: refitting makes the column a within-population quantity, so the same sequence
# would score differently alongside a different set of transcripts.
# g4mer_mean and g4mer_frac_above have no length bias and enter raw.
#
# Requires the g4mer conda env (torch + transformers). Inference is CPU-fine at this scale.
# Validated by 33g_.
suppressPackageStartupMessages({ library(tidyverse); library(here); library(Biostrings) })

G4_STRIDE <- c(utr5 = 10L, cds = 20L, utr3 = 20L)
G4_WINDOW <- 70L
G4_PY     <- "/opt/homebrew/Caskroom/mambaforge/base/envs/g4mer/bin/python"
G4_COLS   <- as.vector(outer(c("g4mer_max_resid_", "g4mer_mean_", "g4mer_frac_above_"),
                             c("utr5", "cds", "utr3"), paste0))

#' Score one region's sequences with G4mer and return its summary table.
run_g4mer <- function(ids, seqs, region, workdir = NULL) {
  stopifnot("unknown region" = region %in% names(G4_STRIDE))
  d <- workdir %||% file.path(tempdir(), paste0("g4_", region, "_", as.integer(Sys.time())))
  dir.create(file.path(d, "fasta"), showWarnings = FALSE, recursive = TRUE)
  fa <- file.path(d, "fasta", paste0(region, "_sequences.fa"))
  writeXStringSet(DNAStringSet(setNames(toupper(seqs), ids)), fa)
  stopifnot("the g4mer python environment was not found" = file.exists(G4_PY))
  cmd <- paste(shQuote(G4_PY), shQuote(here("g4mer", "02_run_g4mer.py")),
               "--region", region, "--fasta", shQuote(fa), "--outdir", shQuote(d),
               "--basedir", shQuote(d), "--window", G4_WINDOW, "--stride", G4_STRIDE[[region]],
               "--device cpu")
  st <- system(cmd, ignore.stdout = TRUE)
  stopifnot("G4mer inference failed" = st == 0)
  f <- file.path(d, paste0("g4mer_summary_", region, ".tsv"))
  stopifnot("G4mer wrote no summary" = file.exists(f))
  read_tsv(f, show_col_types = FALSE)
}

#' Residualise g4mer_max on log2(region length + 1) + region GC, using the saved 231 coefficients.
g4_residual <- function(y, region_length, region_gc, region, coefs) {
  cf <- coefs[[region]]
  stopifnot("no saved g4 residual coefficients for this region" = !is.null(cf), length(cf) == 3)
  l2 <- log2(region_length + 1)
  r <- y - (cf[1] + cf[2] * l2 + cf[3] * region_gc)
  r[!is.finite(y) | !is.finite(l2) | !is.finite(region_gc)] <- NA_real_
  r
}

#' Compute all 9 G4 features.
#' @param d tibble with columns id, utr5, cds, utr3
#' @param gc named list/tibble of per-region GC percentages (utr5_gc, cds_gc, utr3_gc), which
#'   must be the SAME values the model's gc features carry, not recomputed differently.
g4_features <- function(d, gc) {
  stopifnot(all(c("id", "utr5", "cds", "utr3") %in% names(d)),
            all(c("utr5_gc", "cds_gc", "utr3_gc") %in% names(gc)))
  coefs <- readRDS(file.path(here("output", "translational_control_overview"),
                             "g4_resid_coefs_mdamb231.rds"))
  out <- tibble(id = d$id)
  for (rg in names(G4_STRIDE)) {
    s <- run_g4mer(d$id, d[[rg]], rg)
    s <- s[match(d$id, s$transcript_id_clean), ]
    stopifnot("G4mer did not score every sequence" = !anyNA(s$transcript_id_clean))
    out[[paste0("g4mer_mean_", rg)]]       <- s$g4mer_mean
    out[[paste0("g4mer_frac_above_", rg)]] <- s$g4mer_frac_above
    out[[paste0("g4mer_max_resid_", rg)]]  <-
      g4_residual(s$g4mer_max, s$region_length, gc[[paste0(rg, "_gc")]], rg, coefs)
  }
  stopifnot("a G4 feature is missing" = all(G4_COLS %in% names(out)),
            "g4mer_mean outside [0,1]" =
              all(map_lgl(out[grep("^g4mer_mean_", names(out))], ~ all(.x >= 0 & .x <= 1, na.rm = TRUE))),
            "g4mer_frac_above outside [0,1]" =
              all(map_lgl(out[grep("^g4mer_frac_above_", names(out))], ~ all(.x >= 0 & .x <= 1, na.rm = TRUE))))
  out %>% dplyr::select(id, all_of(G4_COLS))
}
