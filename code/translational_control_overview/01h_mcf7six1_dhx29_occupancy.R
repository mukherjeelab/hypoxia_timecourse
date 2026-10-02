# 01h_mcf7six1_dhx29_occupancy.R
# Recompute the Hia 2026 DHX29 occupancy score on the MCF7-SIX1 dominant transcripts.
#
# WHY THIS IS A RECOMPUTE AND NOT A JOIN. The score is not a measurement per transcript; it is
# a codon-weighted average of the Table S8 codon efficiency values over the CDS:
#
#     occupancy = sum(CES[codon] * count) / n_sense_codons
#
# so it is computable for any transcript with a CDS. The saved
# output/predictive_modeling/dhx29_riboseq_occupancy.rds covers only 49.5% of MCF7-SIX1's
# dominant transcripts, but that is not a data limit: notebook 23 restricts scoring to
# precomputed_most_abundant_tx.rds, which is the MDA-MB-231 transcript set. Joining it onto
# MCF7-SIX1 would leave the feature half missing and leaning on train-median imputation for no
# reason, so it is recomputed on MCF7-SIX1's own transcripts instead.
#
# REPRODUCTION FIRST. Before computing anything new, the same function is run over the 231
# transcript set and required to reproduce the stored values exactly. Without that, a
# recomputed column could differ from the 231 column for reasons that have nothing to do with
# the transcript set, and the two cell lines would not be comparable.
#
# Usage: Rscript code/translational_control_overview/01h_mcf7six1_dhx29_occupancy.R [cell]
#   cell = mcf7six1 (default) | hela | dhx29. hela / dhx29 score the transcript set in
#   precomputed_teleman_tx.rds / precomputed_dhx29_tx.rds and write
#   dhx29_riboseq_occupancy_<cell>.rds. The reproduction check against the 231 file runs first
#   in every case.

suppressPackageStartupMessages({
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings); library(dplyr); library(readr); library(here)
})

# ---- 1. Table S8, loaded and replicated exactly as notebook 23 does it -------------------
s8_path <- here("accessories", "csc_data", "hia_dhx29",
                "science.adw0288_tables_s1_to_s11", "science.adw0288_table_s8.csv")
stopifnot("Table S8 not found" = file.exists(s8_path))
stop_codons <- c("TAA", "TAG", "TGA")

s8 <- read_csv(s8_path, show_col_types = FALSE) %>%
  mutate(codon = as.character(codon)) %>%
  filter(!codon %in% stop_codons)

# The same replication assertions notebook 23 applies at load. Nothing is transcribed: the
# AU3 bias is recomputed from the file and checked, so a wrong sheet or a stale download fires
# here rather than producing a plausible-looking column.
s8_au3 <- substr(s8$codon, 3, 3) %in% c("A", "T")
stopifnot(
  "S8 must contain all 61 sense codons" = nrow(s8) == 61,
  "CES is a fold change and must be positive" = all(s8$Fold_Change > 0),
  "CES range implausible - check the S8 file/column" =
    between(max(s8$Fold_Change), 1.1, 3.0) && between(min(s8$Fold_Change), 0.3, 0.95),
  "Table S8 replication failed: top CES codons are not AU3-biased" =
    mean(s8_au3[order(-s8$Fold_Change)][1:10]) >= 0.8,
  "Table S8 replication failed: bottom CES codons are not GC3-biased" =
    mean(s8_au3[order(s8$Fold_Change)][1:10]) <= 0.4
)
ces_vec <- setNames(s8$Fold_Change, s8$codon)
cat(sprintf("Table S8: %d sense codons, CES %.3f-%.3f\n",
            nrow(s8), min(s8$Fold_Change), max(s8$Fold_Change)))

# ---- 2. The scoring function, character-for-character from notebook 23 -------------------
compute_occupancy <- function(seq_char) {
  n_codons <- nchar(seq_char) / 3
  if (n_codons < 10) return(NA_real_)
  starts <- seq(1, nchar(seq_char) - 2, by = 3)
  codons <- substring(seq_char, starts, starts + 2)
  codons <- codons[!codons %in% stop_codons]
  if (length(codons) == 0) return(NA_real_)
  sum(ces_vec[codons], na.rm = TRUE) / length(codons)
}

txdb <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
cds_all <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38,
                                 cdsBy(txdb, by = "tx", use.names = TRUE))
names(cds_all) <- sub("\\..*", "", names(cds_all))

score_set <- function(keep_tx, label) {
  s <- cds_all[names(cds_all) %in% keep_tx]
  s <- s[width(s) %% 3 == 0]                       # notebook 23's frame filter
  v <- vapply(as.character(s), compute_occupancy, numeric(1))
  names(v) <- names(s)
  cat(sprintf("%-10s requested %6d | with in-frame CDS %6d | scored %6d\n",
              label, length(keep_tx), length(s), sum(!is.na(v))))
  tibble(transcript_id_clean = names(v), dhx29_occupancy_score = unname(v)) %>%
    filter(!is.na(dhx29_occupancy_score))
}

# ---- 3. Reproduction gate against the stored 231 values -----------------------------------
stored <- readRDS(here("output", "predictive_modeling", "dhx29_riboseq_occupancy.rds"))
tx231  <- as.character(unlist(readRDS(here("output", "predictive_modeling",
                          "precomputed_most_abundant_tx.rds"))$transcript_id_clean))
repro <- score_set(tx231, "231")
chk <- stored %>%
  dplyr::select(transcript_id_clean, stored = dhx29_occupancy_score) %>%
  inner_join(repro, by = "transcript_id_clean")
dev <- max(abs(chk$stored - chk$dhx29_occupancy_score))
cat(sprintf("\nreproduction of stored 231 scores: n=%d, max|dev| %.3e\n", nrow(chk), dev))
stopifnot(
  "the recompute does not reproduce the stored 231 occupancy scores - do not trust it on MCF7-SIX1" =
    dev < 1e-12,
  "the reproduction check covered implausibly few transcripts" = nrow(chk) > 7000
)

# ---- 4. the target transcript set ------------------------------------------------------------
CELL <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "mcf7six1" }
stopifnot("cell must be mcf7six1, hela, dhx29 or mane" = CELL %in% c("mcf7six1", "hela", "dhx29", "mane", "mdamb231"))
m6 <- if (CELL == "mcf7six1") readRDS(here("output", "predictive_modeling",
                                            "feature_matrix_mcf7six1_external_stability.rds")) else
  readRDS(here("output", "predictive_modeling",
               c(hela = "precomputed_teleman_tx.rds", dhx29 = "precomputed_dhx29_tx.rds",
                 mane = "precomputed_mane_tx.rds",
           mdamb231 = "precomputed_mdamb231_tx.rds")[[CELL]]))
stopifnot("the MCF7-SIX1 matrix is not one transcript per gene" =
            !any(duplicated(m6$gene_id_clean)))
mcf <- score_set(m6$transcript_id_clean, CELL)

out <- mcf %>%
  left_join(m6 %>% transmute(transcript_id_clean,
                             ensembl_gene = gene_id_clean, symbol), by = "transcript_id_clean")
stopifnot("a scored MCF7-SIX1 transcript lost its gene mapping" = !any(is.na(out$ensembl_gene)),
          "scores must be positive - CES is a fold change" = all(out$dhx29_occupancy_score > 0))

# Same schema as the 231 file, minus `group`, which notebook 23 adds for its own plots.
saveRDS(out, here("output", "predictive_modeling",
                  paste0("dhx29_riboseq_occupancy_", CELL, ".rds")))

cov_before <- mean(m6$transcript_id_clean %in% sub("[.].*", "", stored$transcript_id_clean))
cov_after  <- mean(m6$transcript_id_clean %in% out$transcript_id_clean)
cat(sprintf("\n%s coverage of the occupancy feature:\n  joining the 231 table: %.1f%%\n  recomputed here      : %.1f%%\n", CELL,
            100 * cov_before, 100 * cov_after))
cat(sprintf("for reference, the 231 matrix carries it at 83.8%%\n"))
cat(sprintf("\nscore distribution - 231 stored: median %.4f | recomputed: median %.4f\n",
            median(stored$dhx29_occupancy_score), median(out$dhx29_occupancy_score)))
cat("wrote", paste0("dhx29_riboseq_occupancy_", CELL, ".rds"), "\n")
