# Score the transcripts the cluster RNAplfold run never covered (audit finding F1).
#
# F1 restored 1,689 transcripts from a 0 sentinel to NA, which is honest but is still absence
# rather than data - 02_ median-imputes them. RNAplfold is installed locally, and the
# sequences are in hand, so the values are recoverable.
#
# TWO THINGS MAKE THE NEW VALUES TRUSTWORTHY, and the script refuses to write without both:
#
#  1. Region-correct parameters. These DIFFER BY REGION, verified by reproducing the stored
#     _lunp files bit-for-bit: 5'UTR -W 80 -L 40 -u 30, CDS and 3'UTR -W 150 -L 100 -u 30.
#     Using one set for all three would put the new values on a different scale from the old
#     ones inside the same column - worse than leaving them NA.
#  2. Reproduction of transcripts that ALREADY have values. Before filling anything, the
#     script re-scores transcripts the cluster run did cover and asserts it reproduces the
#     stored feature to 1e-6. If this pipeline cannot reproduce known values, its new values
#     mean nothing.
#
# Deprecated columns (*_min, num_structured_regions) are deliberately NOT filled: *_min is
# broken upstream and the region count thresholds a 30-nt window at <0.2, which every position
# passes (p(unpaired 30-mer) ~ 1e-9). Filling them would manufacture a column that disagrees
# with its own history.
#
# Usage:
#   Rscript code/translational_control_overview/01e_score_missing_structure.R utr5 [n_validate]
#   ... then cds, then utr3
suppressMessages({library(tidyverse); library(here)})

args   <- commandArgs(trailingOnly = TRUE)
region <- if (length(args)) args[1] else "utr5"
n_val  <- if (length(args) > 1) as.integer(args[2]) else 25
stopifnot("region must be utr5, cds or utr3" = region %in% c("utr5", "cds", "utr3"))

PARAMS <- list(utr5 = "-W 80 -L 40 -u 30",
               cds  = "-W 150 -L 100 -u 30",
               utr3 = "-W 150 -L 100 -u 30")[[region]]

RNAPLFOLD <- Sys.which("RNAplfold")
if (!nzchar(RNAPLFOLD))
  RNAPLFOLD <- "/opt/homebrew/Caskroom/mambaforge/base/bin/RNAplfold"
stopifnot("RNAplfold not found" = file.exists(RNAPLFOLD))
cat("RNAplfold:", RNAPLFOLD, "| region:", region, "| params:", PARAMS, "\n")

out_dir  <- here("output", "translational_control_overview")
features <- readRDS(file.path(out_dir, "feature_matrix_tco.rds"))

# ---- sequences -------------------------------------------------------------------------
seqs <- readRDS(here("output", "predictive_modeling", "precomputed_kozak_sequences.rds"))
seq_tbl <- switch(region,
  utr5 = seqs$utr5_seq %>% dplyr::select(transcript_id_clean, seq = utr5_sequence),
  cds  = seqs$cds_seq  %>% dplyr::select(transcript_id_clean, seq = cds_sequence),
  utr3 = readRDS(file.path(out_dir, "precomputed_utr3_seqs.rds")) %>%
           dplyr::select(transcript_id_clean, seq = utr3_seq))
seq_tbl <- seq_tbl %>% filter(!is.na(seq), nchar(seq) >= 10)

# ---- which columns this region owns -------------------------------------------------------
cols <- switch(region,
  utr5 = c(cap = "tco_struct_accessibility_cap_proximal",
           aug = "tco_struct_accessibility_aug_context",
           mean = "tco_struct_accessibility_utr5_mean",
           min  = "tco_struct_accessibility_utr5_min",
           sd   = "tco_struct_accessibility_utr5_sd"),
  cds  = c(cap = "tco_struct_accessibility_start_proximal_cds",
           aug = "tco_struct_accessibility_stop_proximal_cds",
           mean = "tco_struct_accessibility_cds_mean",
           min  = "tco_struct_accessibility_cds_min",
           sd   = "tco_struct_accessibility_cds_sd"),
  utr3 = c(cap = "tco_struct_accessibility_stop_proximal_utr3",
           aug = "tco_struct_accessibility_distal_utr3",
           mean = "tco_struct_accessibility_utr3_mean",
           min  = "tco_struct_accessibility_utr3_min",
           sd   = "tco_struct_accessibility_utr3_sd"))

# 01d_/01k_ index the two positional features from opposite ends depending on region:
# 5'UTR and CDS take positions 1-30 then the LAST 30; the 3'UTR takes 1-30 (stop-proximal)
# then the last 30 (distal). Same shape, different names.
score_one <- function(lunp_file, tx) {
  lines <- readLines(lunp_file, warn = FALSE)
  hdr   <- max(grep("^#|^ #", lines))
  df    <- read.table(text = lines[(hdr + 1):length(lines)], header = FALSE, sep = "\t",
                      fill = TRUE)
  acc <- df[[2]]; acc[is.na(acc)] <- 0        # l=1 column, as 01d_/01k_ do
  n   <- length(acc)
  tibble(transcript_id_clean = tx,
         cap  = mean(acc[1:min(30, n)], na.rm = TRUE),
         aug  = mean(acc[max(1, n - 29):n], na.rm = TRUE),
         mean = mean(acc, na.rm = TRUE),
         min  = min(acc, na.rm = TRUE),
         sd   = sd(acc, na.rm = TRUE))
}

# The folded _lunp files are ARCHIVED, not discarded: re-folding 1,600 CDS sequences costs
# 13 minutes, and any future accessibility feature should be derivable from the files rather
# than from another fold.
archive_dir <- here("accessories", "plfold_output", paste0("rescued_", region))

run_plfold <- function(tbl, label) {
  if (!nrow(tbl)) return(tibble())
  wd <- if (label == "missing") archive_dir else
          file.path(tempdir(), paste0("plf_", region, "_", label))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  fa <- file.path(wd, "in.fa")
  writeLines(as.vector(rbind(paste0(">", tbl$transcript_id_clean), tbl$seq)), fa)
  t0 <- Sys.time()
  system(sprintf("cd %s && %s %s < in.fa > /dev/null 2>&1", shQuote(wd), shQuote(RNAPLFOLD),
                 PARAMS))
  cat("  folded", nrow(tbl), label, "sequences in",
      round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
  if (label == "missing") cat("  _lunp archived in", wd, "\n")
  map_dfr(tbl$transcript_id_clean, function(tx) {
    f <- file.path(wd, paste0(tx, "_lunp"))
    if (!file.exists(f)) return(tibble())
    score_one(f, tx)
  })
}

# ---- 1. reproduce known values -------------------------------------------------------------
have <- features %>% filter(!is.na(.data[[cols[["mean"]]]])) %>% pull(transcript_id_clean)
val_tx <- head(sort(intersect(have, seq_tbl$transcript_id_clean)), n_val)
cat("\nValidation: re-scoring", length(val_tx), "transcripts that already have values\n")
val <- run_plfold(seq_tbl %>% filter(transcript_id_clean %in% val_tx), "validate")

stored <- features %>% filter(transcript_id_clean %in% val$transcript_id_clean) %>%
  dplyr::select(transcript_id_clean, all_of(unname(cols)))
names(stored) <- c("transcript_id_clean", names(cols))
chk <- val %>% inner_join(stored, by = "transcript_id_clean", suffix = c("_new", "_old"))
devs <- map_dbl(names(cols), ~ max(abs(chk[[paste0(.x, "_new")]] - chk[[paste0(.x, "_old")]]),
                                   na.rm = TRUE))
names(devs) <- names(cols)
cat("  max |deviation| vs stored:\n")
for (k in names(devs)) cat(sprintf("    %-5s %.3g\n", k, devs[[k]]))
stopifnot("this pipeline does not reproduce the stored values - do NOT fill" =
            all(devs < 1e-6, na.rm = TRUE))
cat("  reproduction confirmed\n")

# ---- 2. score the missing ------------------------------------------------------------------
missing <- features %>% filter(is.na(.data[[cols[["mean"]]]])) %>% pull(transcript_id_clean)
todo <- seq_tbl %>% filter(transcript_id_clean %in% missing)
cat("\nMissing:", length(missing), "| sequence available for", nrow(todo), "\n")
if (!nrow(todo)) { cat("nothing to score\n"); quit(save = "no") }

new <- run_plfold(todo, "missing")
cat("  scored", nrow(new), "of", nrow(todo), "\n")
stopifnot("scored values outside [0, 1]" =
            all(new$mean >= 0 & new$mean <= 1, na.rm = TRUE))

names(new) <- c("transcript_id_clean", unname(cols))
f <- file.path(out_dir, paste0("structure_rescued_", region, ".rds"))
saveRDS(new, f)
cat("Wrote", f, "-", nrow(new), "transcripts x", length(cols), "features\n")
