# 01c_generate_fasta_teleman.R
# Generate CDS, 5'UTR and 3'UTR FASTA for the Teleman (HeLa) dominant transcript set,
# for ViennaRNA structure prediction.
#
# Input:  output/predictive_modeling/precomputed_teleman_tx.rds  (written by notebook 78)
# Output: <work_dir>/teleman_{cds,utr5,utr3}_sequences.fa
#
# Needs BSgenome.Hsapiens.UCSC.hg38, which is often absent from cluster envs. It is
# usually easiest to run this LOCALLY and copy the resulting FASTA to the cluster;
# RNAplfold/RNAfold need only ViennaRNA, no R.
#
#   # locally, from the repo root (paths below line up with the repo layout):
#   VIENNA_WORK_DIR="$PWD" VIENNA_OUT_DIR="$PWD/output/viennarna" \
#     Rscript cluster_scripts/viennarna/01c_generate_fasta_teleman.R
#
#   # or on the cluster, if BSgenome is available there:
#   Rscript 01c_generate_fasta_teleman.R

library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

# Defaults target the cluster; override with env vars to run locally.
work_dir <- Sys.getenv("VIENNA_WORK_DIR", "/beevol/home/matlinka/timecourse/viennarna")
out_dir  <- Sys.getenv("VIENNA_OUT_DIR",  work_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("reading inputs from:", work_dir, "\n")
cat("writing FASTA to   :", out_dir, "\n")

txdb_path       <- file.path(work_dir, "accessories", "human", "txdb.gencode49.sqlite")
tx_mapping_path <- file.path(work_dir, "output", "predictive_modeling",
                             "precomputed_teleman_tx.rds")

txdb   <- loadDb(txdb_path)
genome <- BSgenome.Hsapiens.UCSC.hg38

tx_map <- readRDS(tx_mapping_path)
cat("Teleman transcript set:", nrow(tx_map), "transcripts\n")
stopifnot("mapping must carry transcript_id_clean" =
            "transcript_id_clean" %in% names(tx_map))

# ---------------------------------------------------------------------------
# Helper: extract a region, restrict to the transcript set, write FASTA
# ---------------------------------------------------------------------------
write_region <- function(grl, label) {
  seqs <- extractTranscriptSeqs(genome, grl)
  names(seqs) <- sub("\\..*", "", names(seqs))
  keep <- names(seqs) %in% tx_map$transcript_id_clean
  seqs <- seqs[keep]
  seqs <- seqs[width(seqs) > 0]
  # ViennaRNA reads RNA; T -> U is not required by RNAplfold but keeps the FASTA honest
  out <- file.path(out_dir, paste0("teleman_", label, "_sequences.fa"))
  writeXStringSet(seqs, out)
  cat(sprintf("%-5s %6d sequences  median width %5.0f  -> %s\n",
              label, length(seqs), median(width(seqs)), basename(out)))
  invisible(seqs)
}

cds_seqs  <- write_region(cdsBy(txdb, by = "tx", use.names = TRUE),        "cds")
utr5_seqs <- write_region(fiveUTRsByTranscript(txdb, use.names = TRUE),    "utr5")
utr3_seqs <- write_region(threeUTRsByTranscript(txdb, use.names = TRUE),   "utr3")

cat("\nCoverage of the requested transcript set:\n")
cat("  CDS  :", round(100 * length(cds_seqs)  / nrow(tx_map), 1), "%\n")
cat("  5'UTR:", round(100 * length(utr5_seqs) / nrow(tx_map), 1), "%\n")
cat("  3'UTR:", round(100 * length(utr3_seqs) / nrow(tx_map), 1), "%\n")
cat("\nDone. Next: 02c_run_rnaplfold_teleman.sh and 03c_run_rnafold_teleman.sh\n")
