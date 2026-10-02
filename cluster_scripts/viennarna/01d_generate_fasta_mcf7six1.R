# 01d_generate_fasta_mcf7six1.R
# Region FASTAs for the MCF7-SIX1 dominant-transcript set, for RNAplfold/RNAfold.
#
# Adapted from 01c_generate_fasta_teleman.R so the MCF7-SIX1 structure-feature universe is
# built exactly the way the MDA-MB-231 and Teleman ones were. Only the transcript set differs.
#
# The MCF7-SIX1 dominant transcript per gene was already chosen in
# code/predictive_modeling/43_feature_extraction_mcf7six1_codon.Rmd, from the MCF7-SIX1 siCTRL
# normoxia Salmon samples (NM2023_0104/0114/0124), with longest-CDS as the tie-break. That
# choice is read back out of its saved matrix rather than re-derived, so this cannot disagree
# with the matrix the codon and stability features were computed on.
#
# ONLY THE MISSING TRANSCRIPTS ARE WRITTEN. An RNAplfold _lunp file is a property of the
# sequence, not of the cell line that selected it, and regeneration was verified byte-identical
# against the stored MDA-MB-231 files, so a transcript dominant in both cell lines is reused
# rather than recomputed. That takes the job from 11,565 transcripts to about 5,900.
#
# Usage: Rscript cluster_scripts/viennarna/01d_generate_fasta_mcf7six1.R [cell]
#   cell = mcf7six1 (default) | hela. HeLa reads precomputed_teleman_tx.rds and never writes it.

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(dplyr)
  library(here)
})

CELL <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "mcf7six1" }
stopifnot("cell must be mcf7six1, hela, dhx29 or mane" = CELL %in% c("mcf7six1", "hela", "dhx29", "mane", "mdamb231"))
txdb_path <- here("accessories", "human", "txdb.gencode49.sqlite")
m6_path   <- here("output", "predictive_modeling",
                  "feature_matrix_mcf7six1_external_stability.rds")
plf       <- here("accessories", "plfold_output", "rnaplfold_output")
out_dir   <- here("accessories", "plfold_output", paste0(CELL, "_fasta"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot("TxDb not found - run 00_build_txdb_v49.R first" = file.exists(txdb_path))

txdb   <- loadDb(txdb_path)
genome <- BSgenome.Hsapiens.UCSC.hg38

if (CELL %in% c("hela", "dhx29", "mane", "mdamb231")) {
  # hela: Roiuk-dominant HeLa transcript. dhx29: the representative transcript per gene from
  # code/dhx29_reanalysis/03 (most footprints, then MANE, CDS length, ID), HEK293T.
  # mane: MANE Select for the Weber DAP5 genes (10w_weber_dap5_outcome_prep.Rmd).
  rds <- c(hela = "precomputed_teleman_tx.rds", dhx29 = "precomputed_dhx29_tx.rds",
           mane = "precomputed_mane_tx.rds",
           mdamb231 = "precomputed_mdamb231_tx.rds")[[CELL]]
  tx_all <- unique(readRDS(here("output", "predictive_modeling", rds))$transcript_id_clean)
  cat(CELL, "transcripts:", length(tx_all), "\n")
} else {
stopifnot("MCF7-SIX1 feature matrix not found" = file.exists(m6_path))
m6 <- readRDS(m6_path)
stopifnot("the MCF7-SIX1 matrix is not one transcript per gene" =
            !any(duplicated(m6$gene_id_clean)),
          "the MCF7-SIX1 matrix has no transcript_id_clean" =
            "transcript_id_clean" %in% names(m6))
tx_all <- unique(m6$transcript_id_clean)
cat("MCF7-SIX1 dominant transcripts:", length(tx_all), "\n")

# Save the mapping in the same shape as precomputed_most_abundant_tx.rds, so any later step can
# use it the way the 231 and Teleman steps use theirs.
tx_map <- m6 %>% dplyr::select(transcript_id_clean, gene_id_clean, gene_id, symbol)
saveRDS(tx_map, here("output", "predictive_modeling", "precomputed_mcf7six1_tx.rds"))
cat("wrote precomputed_mcf7six1_tx.rds\n")
}

region_dir <- c(utr5 = "rnaplfold_output", cds = "cds_rnaplfold_output",
                utr3 = "utr3_rnaplfold_output")

write_region <- function(grl, label) {
  seqs <- extractTranscriptSeqs(genome, grl)
  names(seqs) <- sub("\\..*", "", names(seqs))
  seqs <- seqs[names(seqs) %in% tx_all]
  seqs <- seqs[width(seqs) > 0]

  have <- sub("_lunp$", "", list.files(file.path(plf, region_dir[[label]]), pattern = "_lunp$"))
  need <- setdiff(names(seqs), have)
  reuse <- length(seqs) - length(need)
  seqs_need <- seqs[names(seqs) %in% need]

  out <- file.path(out_dir, paste0(CELL, "_", label, "_sequences.fa"))
  writeXStringSet(seqs_need, out)
  cat(sprintf("%-5s extractable %6d | already have _lunp %6d | to compute %6d | %10.0f nt -> %s\n",
              label, length(seqs), reuse, length(seqs_need),
              sum(as.numeric(width(seqs_need))), basename(out)))
  invisible(seqs_need)
}

cds  <- write_region(cdsBy(txdb, by = "tx", use.names = TRUE),      "cds")
utr5 <- write_region(fiveUTRsByTranscript(txdb, use.names = TRUE),  "utr5")
utr3 <- write_region(threeUTRsByTranscript(txdb, use.names = TRUE), "utr3")

cat("\nNext: bash cluster_scripts/viennarna/02d_run_rnaplfold_mcf7six1.sh\n")
