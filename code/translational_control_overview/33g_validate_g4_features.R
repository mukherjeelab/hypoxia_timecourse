# 33g_validate_g4_features.R - does 33d_ reproduce the matrix's G4 features?
#
# Split into two checks, because they fail for different reasons:
#   A. AGGREGATION/RESIDUAL - take the STORED G4mer summaries and run only the residual step.
#      This tests the saved coefficients, the log2(length + 1) term and the GC column used,
#      with no model inference involved. It is exact, and covers every transcript.
#   B. INFERENCE - re-score a handful of real transcripts with G4mer here and require the
#      summary to match the stored one. This tests the window, the stride and the model.
#      B is slow, so it runs on few transcripts; A is the one that covers the arithmetic.
#   Rscript code/translational_control_overview/33g_validate_g4_features.R [n_infer]
suppressPackageStartupMessages({ library(tidyverse); library(here); library(Biostrings) })
source(here("code", "translational_control_overview", "33d_g4_features.R"))

args <- commandArgs(trailingOnly = TRUE)
N_INFER <- if (length(args)) as.integer(args[1]) else 6L
TOL <- 1e-6
OUT_DIR <- here("output", "translational_control_overview")
mat <- readRDS(file.path(OUT_DIR, "feature_matrix_tco.rds"))
coefs <- readRDS(file.path(OUT_DIR, "g4_resid_coefs_mdamb231.rds"))

# ---- A. residual arithmetic, on the stored summaries --------------------------------
cat("A. residual step against the stored G4mer summaries\n")
res <- map_dfr(names(G4_STRIDE), function(rg) {
  s <- read_tsv(here("output", "g4mer", paste0("g4mer_summary_", rg, ".tsv")), show_col_types = FALSE)
  j <- mat %>% dplyr::select(transcript_id_clean, gc = all_of(paste0(rg, "_gc")),
                             ref_resid = all_of(paste0("g4mer_max_resid_", rg)),
                             ref_mean  = all_of(paste0("g4mer_mean_", rg)),
                             ref_frac  = all_of(paste0("g4mer_frac_above_", rg))) %>%
    inner_join(s, by = "transcript_id_clean")
  j$new_resid <- g4_residual(j$g4mer_max, j$region_length, j$gc, rg, coefs)
  tibble(region = rg, n = sum(!is.na(j$new_resid) & !is.na(j$ref_resid)),
         resid_dev = max(abs(j$new_resid - j$ref_resid), na.rm = TRUE),
         mean_dev  = max(abs(j$g4mer_mean - j$ref_mean), na.rm = TRUE),
         frac_dev  = max(abs(j$g4mer_frac_above - j$ref_frac), na.rm = TRUE))
})
print(as.data.frame(res), digits = 4)
if (any(c(res$resid_dev, res$mean_dev, res$frac_dev) > TOL)) { cat("\nFAILED residual check\n"); quit(status = 1) }
cat(sprintf("   PASS: residual, mean and frac_above reproduce the matrix (tol %.0e)\n\n", TOL))

# ---- B. model inference, on a few transcripts ---------------------------------------
cat("B. re-running G4mer inference on", N_INFER, "transcripts per region\n")
fa <- function(f) { s <- readDNAStringSet(here("g4mer", "fasta", f)); tibble(id = names(s), seq = as.character(s)) }
set.seed(9)
infer <- map_dfr(names(G4_STRIDE), function(rg) {
  src <- fa(paste0(rg, "_sequences.fa")) %>% filter(id %in% mat$transcript_id_clean) %>%
    slice_sample(n = N_INFER)
  got <- run_g4mer(src$id, src$seq, rg)
  ref <- read_tsv(here("output", "g4mer", paste0("g4mer_summary_", rg, ".tsv")), show_col_types = FALSE)
  j <- got %>% inner_join(ref, by = "transcript_id_clean", suffix = c("_new", "_ref"))
  tibble(region = rg, n = nrow(j),
         windows_match = all(j$n_windows_new == j$n_windows_ref),
         max_dev  = max(abs(j$g4mer_max_new  - j$g4mer_max_ref)),
         mean_dev = max(abs(j$g4mer_mean_new - j$g4mer_mean_ref)))
})
print(as.data.frame(infer), digits = 4)
if (!all(infer$windows_match)) { cat("\nFAILED: window count differs - the stride is wrong\n"); quit(status = 1) }
if (any(c(infer$max_dev, infer$mean_dev) > 1e-4)) { cat("\nFAILED: inference does not reproduce stored scores\n"); quit(status = 1) }
cat("\nPASS: G4 features reproduce the stored matrix and the stored inference.\n")
