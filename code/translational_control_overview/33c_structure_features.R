# 33c_structure_features.R - the 15 RNAplfold accessibility features, computed from sequence.
#
# Each region is folded IN ISOLATION, with the parameters that 02e_validate_rnaplfold_params.R
# recovered from the stored output and verified byte-identical on ViennaRNA 2.7.2:
#     5'UTR      -W 80  -L 40  -u 30
#     CDS, 3'UTR -W 150 -L 100 -u 30
# They are NOT the same for all three regions. A single global setting reproduces the CDS and
# 3'UTR and silently mis-generates every 5'UTR, so the per-region table is the definition.
#
# PARSING. An RNAplfold _lunp file's second line is " #i$\tl=1\t..." - it begins with a SPACE.
# The frozen pipeline tested for a header with str_starts(lines, "#"), missed that line, parsed
# it as data, got an all-NA row and replaced the NA with 0, which reads as "maximally
# structured". That phantom zero is the F1/F5 defect. Here the header test allows leading
# whitespace, and a row that still parses to NA is an error rather than a zero.
#
# Validated by 33f_, which requires these functions to reproduce the stored matrix.
suppressPackageStartupMessages({ library(tidyverse); library(here); library(Biostrings) })

PLFOLD_FLAGS <- c(utr5 = "-W 80 -L 40 -u 30",
                  cds  = "-W 150 -L 100 -u 30",
                  utr3 = "-W 150 -L 100 -u 30")
# region -> the five column names, in the order feats() returns them.
STRUCT_COLS <- list(
  utr5 = c(prox = "tco_struct_accessibility_cap_proximal",
           dist = "tco_struct_accessibility_aug_context",
           mean = "tco_struct_accessibility_utr5_mean",
           min  = "tco_struct_accessibility_utr5_min",
           sd   = "tco_struct_accessibility_utr5_sd"),
  cds  = c(prox = "tco_struct_accessibility_start_proximal_cds",
           dist = "tco_struct_accessibility_stop_proximal_cds",
           mean = "tco_struct_accessibility_cds_mean",
           min  = "tco_struct_accessibility_cds_min",
           sd   = "tco_struct_accessibility_cds_sd"),
  utr3 = c(prox = "tco_struct_accessibility_stop_proximal_utr3",
           dist = "tco_struct_accessibility_distal_utr3",
           mean = "tco_struct_accessibility_utr3_mean",
           min  = "tco_struct_accessibility_utr3_min",
           sd   = "tco_struct_accessibility_utr3_sd"))

#' Read the l=1 (single-nucleotide unpaired probability) column of an RNAplfold _lunp file.
read_lunp_l1 <- function(fp) {
  lines <- readLines(fp, warn = FALSE)
  hdr <- max(grep("^[[:space:]]*#", lines))      # leading space allowed - see the note above
  dl  <- lines[(hdr + 1):length(lines)]
  dl  <- dl[nzchar(trimws(dl))]
  v <- suppressWarnings(as.numeric(vapply(strsplit(dl, "\t", fixed = TRUE),
                                          function(x) if (length(x) >= 2) x[2] else NA_character_,
                                          character(1))))
  stopifnot("an accessibility value failed to parse - do not silently zero-fill it" = !anyNA(v),
            "accessibility outside [0,1]" = all(v >= 0 & v <= 1))
  v
}

#' prox = mean of the first 30 positions, dist = mean of the last 30, plus mean/min/sd.
lunp_feats <- function(acc) {
  n <- length(acc)
  c(prox = mean(acc[1:min(30, n)]), dist = mean(acc[max(1, n - 29):n]),
    mean = mean(acc), min = min(acc), sd = sd(acc))
}

#' Fold one region for many sequences and return the five summaries per sequence.
#' RNAplfold writes <name>_lunp into the working directory, so it is run inside a temp dir.
fold_region <- function(ids, seqs, region, keep_dir = NULL) {
  stopifnot("unknown region" = region %in% names(PLFOLD_FLAGS))
  d <- keep_dir %||% file.path(tempdir(), paste0("plf_", region, "_", as.integer(Sys.time())))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  fa <- file.path(d, "in.fa")
  writeXStringSet(DNAStringSet(setNames(toupper(seqs), ids)), fa)
  old <- setwd(d); on.exit(setwd(old), add = TRUE)
  st <- system(paste("RNAplfold", PLFOLD_FLAGS[[region]], "<", shQuote(fa)),
               ignore.stdout = TRUE, ignore.stderr = TRUE)
  stopifnot("RNAplfold failed - is ViennaRNA on the PATH?" = st == 0)
  out <- map_dfr(ids, function(i) {
    fp <- file.path(d, paste0(i, "_lunp"))
    stopifnot("RNAplfold wrote no _lunp for this sequence" = file.exists(fp))
    as_tibble_row(lunp_feats(read_lunp_l1(fp))) %>% mutate(id = i, .before = 1)
  })
  setNames(out, c("id", unname(STRUCT_COLS[[region]][c("prox","dist","mean","min","sd")])))
}

#' Compute all 15 structure features.
#' @param d tibble with columns id, utr5, cds, utr3
structure_features <- function(d) {
  stopifnot(all(c("id", "utr5", "cds", "utr3") %in% names(d)))
  res <- purrr::reduce(
    purrr::imap(list(utr5 = d$utr5, cds = d$cds, utr3 = d$utr3),
                function(sq, rg) fold_region(d$id, sq, rg)),
    dplyr::full_join, by = "id")
  cols <- unname(unlist(STRUCT_COLS))
  stopifnot("a structure feature is missing" = all(cols %in% names(res)),
            "accessibility outside [0,1]" =
              all(map_lgl(res[cols], ~ all(.x >= 0 & .x <= 1, na.rm = TRUE))),
            # The phantom-zero signature: a min of exactly 0 everywhere was the original defect.
            "every min is exactly 0 - the phantom-zero defect has reappeared" =
              !all(res[[STRUCT_COLS$cds[["min"]]]] == 0))
  res %>% dplyr::select(id, all_of(cols))
}
