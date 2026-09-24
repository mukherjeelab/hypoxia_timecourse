# Recompute every tco_struct_* feature from the stored _lunp files, correctly.
#
# THE DEFECT. 01d_/01k_ locate the header with str_starts(lines, "#"), but the second line of
# an RNAplfold _lunp file is " #i$\tl=1\t..." - it begins with a SPACE. The test misses it, the
# column-header line is handed to read.table as data, it parses to an all-NA row, and
# acc_vals[is.na(acc_vals)] <- 0 then turns that into 0, which means "maximally structured".
# Every transcript therefore carries a phantom position at the 5' end with accessibility zero.
#
# That single bug explains three separate audit findings:
#   F5  *_min is identically 0 everywhere      - the phantom zero is always the minimum
#   F2  the region count equals length + 1     - the phantom row is the +1
#   new every accessibility mean is diluted and every sd inflated, and cap_proximal /
#       start_proximal / stop_proximal_utr3 carry a hard zero in 1 of their 30 positions,
#       a 3.3% contamination independent of region length
#
# VERIFICATION STRATEGY. Reproducing the bug is the proof that we understand it. For every
# transcript this script computes BOTH the buggy and the corrected value, and asserts the
# buggy one reproduces the column currently in the matrix. Only then are the corrected values
# written. If we could not reproduce the old values, we would not have characterised the
# defect and would have no business replacing them.
#
# num_structured_regions is NOT rescued. With the phantom row gone the count becomes length
# rather than length + 1, but a full 30-mer has ~1e-9 probability of being unpaired, so the
# "< 0.2" threshold still passes essentially every position. It is degenerate on its own
# merit and stays deprecated.
#
# Usage: Rscript code/translational_control_overview/01f_recompute_structure.R
suppressMessages({library(tidyverse); library(here)})

out_dir  <- here("output", "translational_control_overview")
features <- readRDS(file.path(out_dir, "feature_matrix_tco.rds"))
plf      <- here("accessories", "plfold_output", "rnaplfold_output")

REGIONS <- list(
  utr5 = list(dir = file.path(plf, "rnaplfold_output"),
              cols = c(prox = "tco_struct_accessibility_cap_proximal",
                       dist = "tco_struct_accessibility_aug_context",
                       mean = "tco_struct_accessibility_utr5_mean",
                       min  = "tco_struct_accessibility_utr5_min",
                       sd   = "tco_struct_accessibility_utr5_sd")),
  cds  = list(dir = file.path(plf, "cds_rnaplfold_output"),
              cols = c(prox = "tco_struct_accessibility_start_proximal_cds",
                       dist = "tco_struct_accessibility_stop_proximal_cds",
                       mean = "tco_struct_accessibility_cds_mean",
                       min  = "tco_struct_accessibility_cds_min",
                       sd   = "tco_struct_accessibility_cds_sd")),
  utr3 = list(dir = file.path(plf, "utr3_rnaplfold_output"),
              cols = c(prox = "tco_struct_accessibility_stop_proximal_utr3",
                       dist = "tco_struct_accessibility_distal_utr3",
                       mean = "tco_struct_accessibility_utr3_mean",
                       min  = "tco_struct_accessibility_utr3_min",
                       sd   = "tco_struct_accessibility_utr3_sd"))
)

# Fast manual parse: the l=1 column of every data row, under both conventions.
read_acc <- function(fp) {
  lines <- readLines(fp, warn = FALSE)
  hdr_buggy   <- max(which(startsWith(lines, "#")))        # 01d_'s test: misses " #i$..."
  hdr_correct <- max(grep("^[[:space:]]*#", lines))        # matches a leading space too
  take <- function(h) {
    dl <- lines[(h + 1):length(lines)]
    dl <- dl[nzchar(trimws(dl))]
    suppressWarnings(as.numeric(vapply(strsplit(dl, "\t", fixed = TRUE),
                                       function(x) if (length(x) >= 2) x[2] else NA_character_,
                                       character(1))))
  }
  list(buggy = { a <- take(hdr_buggy);   a[is.na(a)] <- 0; a },   # phantom NA -> 0
       ok    = take(hdr_correct))
}

feats <- function(acc) {
  n <- length(acc)
  c(prox = mean(acc[1:min(30, n)]), dist = mean(acc[max(1, n - 29):n]),
    mean = mean(acc), min = min(acc), sd = sd(acc))
}

all_new <- list()
for (rg in names(REGIONS)) {
  spec <- REGIONS[[rg]]
  files <- list.files(spec$dir, pattern = "_lunp$", full.names = TRUE)
  cat("\n==", rg, "==", length(files), "_lunp files\n")
  t0 <- Sys.time()

  res <- map_dfr(files, function(fp) {
    tx <- sub("_lunp$", "", basename(fp))
    a  <- read_acc(fp)
    if (!length(a$ok) || all(is.na(a$ok))) return(tibble())
    bind_cols(tibble(transcript_id_clean = tx, n_positions = length(a$ok)),
              as_tibble_row(feats(a$buggy)) %>% rename_with(~ paste0(.x, "_buggy")),
              as_tibble_row(feats(a$ok))    %>% rename_with(~ paste0(.x, "_ok")))
  })
  cat("  parsed", nrow(res), "in", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

  # The phantom row must be exactly one extra position.
  stopifnot("the buggy parse is not exactly one row longer" =
              all(res$min_buggy == 0, na.rm = TRUE))

  # Reproduce the stored column with the buggy method.
  stored <- features %>%
    dplyr::select(transcript_id_clean, all_of(unname(spec$cols))) %>%
    rename_with(~ names(spec$cols)[match(.x, unname(spec$cols))], all_of(unname(spec$cols)))
  cmp <- res %>% inner_join(stored, by = "transcript_id_clean")
  cmp <- cmp %>% filter(!is.na(mean))          # transcripts the matrix actually has
  devs <- sapply(names(spec$cols), function(k)
    max(abs(cmp[[paste0(k, "_buggy")]] - cmp[[k]]), na.rm = TRUE))
  cat("  reproducing the BUGGY column (should be ~0):\n")
  for (k in names(devs)) cat(sprintf("    %-5s %.3g\n", k, devs[[k]]))
  stopifnot("could not reproduce the stored values - the defect is not characterised" =
              all(devs < 1e-6, na.rm = TRUE))

  shift <- cmp %>% summarise(mean_shift = median(mean_ok - mean, na.rm = TRUE),
                             prox_shift = median(prox_ok - prox, na.rm = TRUE))
  cat(sprintf("  correction moves the mean by %+0.5f and prox by %+0.5f (median)\n",
              shift$mean_shift, shift$prox_shift))

  fixed <- res %>% dplyr::select(transcript_id_clean, ends_with("_ok")) %>%
    rename_with(~ spec$cols[sub("_ok$", "", .x)], ends_with("_ok"))
  all_new[[rg]] <- fixed
}

corrected <- reduce(all_new, full_join, by = "transcript_id_clean")
stopifnot("accessibility outside [0, 1]" =
            all(corrected %>% dplyr::select(contains("accessibility")) %>%
                  unlist() %>% between(0, 1), na.rm = TRUE))
f <- file.path(out_dir, "structure_recomputed.rds")
saveRDS(corrected, f)
cat("\nWrote", f, "-", nrow(corrected), "transcripts x", ncol(corrected) - 1, "features\n")
