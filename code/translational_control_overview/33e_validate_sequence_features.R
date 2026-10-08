# 33e_validate_sequence_features.R - does 33b_ reproduce the feature matrix from sequence alone?
#
# 33b_ re-implements the model's sequence-derived features so they can be computed for a plasmid,
# which has no transcript annotation. A re-implementation nobody has checked is a second opinion,
# not a measurement. This runs 33b_ on real transcripts - whose region sequences are the same
# FASTAs the pipeline itself scored - and requires the output to match the stored matrix.
#
# Any column that does not match is a defect in 33b_, and the plasmid numbers must not be read
# until it does. The script reports per-feature deviation and exits non-zero on failure.
#   Rscript code/translational_control_overview/33e_validate_sequence_features.R [n_transcripts]
suppressPackageStartupMessages({ library(tidyverse); library(here); library(Biostrings) })
source(here("code", "translational_control_overview", "33b_sequence_features.R"))

args <- commandArgs(trailingOnly = TRUE)
N    <- if (length(args)) as.integer(args[1]) else 400L
TOL  <- 1e-8

fa <- function(f) { s <- readDNAStringSet(here("g4mer", "fasta", f)); tibble(id = names(s), seq = as.character(s)) }
seqs <- purrr::reduce(list(fa("utr5_sequences.fa") %>% dplyr::rename(utr5 = seq),
                           fa("cds_sequences.fa")  %>% dplyr::rename(cds  = seq),
                           fa("utr3_sequences.fa") %>% dplyr::rename(utr3 = seq)),
                      dplyr::inner_join, by = "id")
cat("transcripts with all three regions:", nrow(seqs), "\n")

mat <- readRDS(file.path(OUT_DIR, "feature_matrix_tco.rds"))
feat_cols <- readRDS(file.path(OUT_DIR,
  "rf_model_data_eif33d_promotes_hypoxia_1hr_lfc0.5_tco_intrinsic_corrected_both_clipaggregate.rds"))$feature_cols
check_cols <- intersect(feat_cols, names(seq_features(seqs[1:2, ])))
cat("sequence-only features to check:", length(check_cols), "of", length(feat_cols), "modelled\n")

# A fixed random sample, so the check is the same set every run.
set.seed(9)
samp <- seqs %>% filter(id %in% mat$transcript_id_clean) %>% slice_sample(n = min(N, nrow(.)))
# The pipeline computes codon features only on in-frame CDS; out-of-frame transcripts are a
# property of the annotation, not of this code, so they are excluded rather than asserted on.
samp <- samp %>% filter(nchar(cds) %% 3 == 0)
cat("validating on", nrow(samp), "transcripts\n\n")

got <- seq_features(samp %>% dplyr::select(id, utr5, cds, utr3))
ref <- mat %>% filter(transcript_id_clean %in% got$id) %>%
  dplyr::select(id = transcript_id_clean, all_of(check_cols))
cmp <- got %>% dplyr::select(id, all_of(check_cols)) %>%
  inner_join(ref, by = "id", suffix = c("_new", "_ref"))

res <- map_dfr(check_cols, function(f) {
  a <- cmp[[paste0(f, "_new")]]; b <- cmp[[paste0(f, "_ref")]]
  both <- !is.na(a) & !is.na(b)
  tibble(feature = f, n_compared = sum(both),
         na_mismatch = sum(xor(is.na(a), is.na(b))),
         max_abs_dev = if (any(both)) max(abs(a[both] - b[both])) else NA_real_,
         spearman = if (sum(both) > 2 && sd(a[both]) > 0 && sd(b[both]) > 0)
                      cor(a[both], b[both], method = "spearman") else NA_real_)
}) %>% mutate(ok = !is.na(max_abs_dev) & max_abs_dev < TOL & na_mismatch == 0) %>%
  arrange(ok, desc(max_abs_dev))

print(as.data.frame(res), digits = 4)
bad <- res %>% filter(!ok)
cat(sprintf("\n%d of %d features reproduce the stored matrix exactly (tol %.0e)\n",
            sum(res$ok), nrow(res), TOL))
if (nrow(bad)) {
  cat("\nFAILED:\n"); print(as.data.frame(bad), digits = 4)
  for (f in head(bad$feature, 4)) {
    a <- cmp[[paste0(f, "_new")]]; b <- cmp[[paste0(f, "_ref")]]
    i <- which(!is.na(a) & !is.na(b) & abs(a - b) > TOL)
    if (length(i)) { cat("\n", f, "- worst examples:\n", sep = "")
      print(head(tibble(id = cmp$id[i], new = a[i], ref = b[i], dev = a[i] - b[i])[order(-abs(a[i] - b[i])), ], 3)) }
  }
  quit(status = 1)
}
cat("PASS: the from-sequence path reproduces every sequence-only feature.\n")
