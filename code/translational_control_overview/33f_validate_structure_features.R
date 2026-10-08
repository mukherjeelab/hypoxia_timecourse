# 33f_validate_structure_features.R - does 33c_ reproduce the matrix's structure features?
#
# Two independent things have to hold, and they fail in different ways, so both are checked:
#   A. REGENERATION - running RNAplfold here, with the per-region flags, must reproduce the
#      stored _lunp file for the same transcript. This tests the parameters and the local
#      ViennaRNA build. A version or flag difference shows up here.
#   B. AGGREGATION  - parsing and summarising must reproduce the stored matrix columns. This
#      tests the header handling and the prox/dist/mean/min/sd definitions.
#
# Check A is run against the ORIGINAL stored output; check B against feature_matrix_tco.rds,
# which holds the corrected (phantom-zero repaired) values.
#   Rscript code/translational_control_overview/33f_validate_structure_features.R [n]
suppressPackageStartupMessages({ library(tidyverse); library(here); library(Biostrings) })
source(here("code", "translational_control_overview", "33c_structure_features.R"))

args <- commandArgs(trailingOnly = TRUE)
N    <- if (length(args)) as.integer(args[1]) else 40L
TOL  <- 1e-6
OUT_DIR <- here("output", "translational_control_overview")
PLF <- here("accessories", "plfold_output", "rnaplfold_output")
STORED_DIR <- c(utr5 = "rnaplfold_output", cds = "cds_rnaplfold_output", utr3 = "utr3_rnaplfold_output")

fa <- function(f) { s <- readDNAStringSet(here("g4mer", "fasta", f)); tibble(id = names(s), seq = as.character(s)) }
seqs <- purrr::reduce(list(fa("utr5_sequences.fa") %>% dplyr::rename(utr5 = seq),
                           fa("cds_sequences.fa")  %>% dplyr::rename(cds  = seq),
                           fa("utr3_sequences.fa") %>% dplyr::rename(utr3 = seq)),
                      dplyr::inner_join, by = "id")
mat  <- readRDS(file.path(OUT_DIR, "feature_matrix_tco.rds"))
cols <- unname(unlist(STRUCT_COLS))

set.seed(9)
samp <- seqs %>% filter(id %in% mat$transcript_id_clean) %>% slice_sample(n = min(N, nrow(.)))
cat("validating structure on", nrow(samp), "transcripts\n\n")

# ---- A. regeneration: our RNAplfold output vs the stored _lunp ----------------------
cat("A. regenerating stored _lunp files\n")
regen <- map_dfr(names(STORED_DIR), function(rg) {
  sq <- samp[[rg]]; ids <- samp$id
  d <- file.path(tempdir(), paste0("regen_", rg)); unlink(d, recursive = TRUE)
  invisible(fold_region(ids, sq, rg, keep_dir = d))
  map_dfr(ids, function(i) {
    stored <- file.path(PLF, STORED_DIR[[rg]], paste0(i, "_lunp"))
    ours   <- file.path(d, paste0(i, "_lunp"))
    if (!file.exists(stored)) return(tibble(region = rg, id = i, status = "no stored file", dev = NA_real_))
    a <- read_lunp_l1(ours); b <- read_lunp_l1(stored)
    if (length(a) != length(b))
      return(tibble(region = rg, id = i, status = sprintf("length %d vs %d", length(a), length(b)), dev = NA_real_))
    tibble(region = rg, id = i, status = "ok", dev = max(abs(a - b)))
  })
})
print(regen %>% group_by(region) %>%
        summarise(n = n(), compared = sum(status == "ok"), missing = sum(status == "no stored file"),
                  max_dev = max(dev, na.rm = TRUE), .groups = "drop") %>% as.data.frame(), digits = 4)
regen_ok <- regen %>% filter(status == "ok")
stopifnot("no stored _lunp files were available to compare against" = nrow(regen_ok) > 0)
if (any(regen_ok$dev > TOL) || any(regen$status != "ok" & regen$status != "no stored file")) {
  cat("\nFAILED regeneration:\n"); print(as.data.frame(regen %>% filter(status != "ok" | dev > TOL) %>% head(10)))
  quit(status = 1)
}
cat(sprintf("   PASS: %d files regenerate to within %.0e\n\n", nrow(regen_ok), TOL))

# ---- B. aggregation: our summaries vs the stored matrix -----------------------------
cat("B. reproducing the matrix columns\n")
got <- structure_features(samp %>% dplyr::select(id, utr5, cds, utr3))
ref <- mat %>% filter(transcript_id_clean %in% got$id) %>%
  dplyr::select(id = transcript_id_clean, all_of(cols))
cmp <- got %>% inner_join(ref, by = "id", suffix = c("_new", "_ref"))

res <- map_dfr(cols, function(f) {
  a <- cmp[[paste0(f, "_new")]]; b <- cmp[[paste0(f, "_ref")]]
  both <- !is.na(a) & !is.na(b)
  tibble(feature = f, n = sum(both), na_mismatch = sum(xor(is.na(a), is.na(b))),
         max_abs_dev = if (any(both)) max(abs(a[both] - b[both])) else NA_real_)
}) %>% mutate(ok = !is.na(max_abs_dev) & max_abs_dev < TOL & na_mismatch == 0) %>% arrange(ok, desc(max_abs_dev))
print(as.data.frame(res), digits = 4)

if (any(!res$ok)) {
  cat("\nFAILED aggregation\n")
  f <- res$feature[!res$ok][1]
  a <- cmp[[paste0(f, "_new")]]; b <- cmp[[paste0(f, "_ref")]]
  i <- order(-abs(a - b))[1:3]
  print(tibble(id = cmp$id[i], new = a[i], ref = b[i], dev = a[i] - b[i]))
  quit(status = 1)
}
cat(sprintf("\nPASS: all %d structure features reproduce the stored matrix (tol %.0e)\n", nrow(res), TOL))
