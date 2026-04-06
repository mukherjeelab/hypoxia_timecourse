# 01b_generate_fasta_cds_utr3.R
# Generate FASTA files of CDS and 3'UTR sequences for ViennaRNA structure prediction
# Run on cluster: bsub < 01b_generate_fasta_cds_utr3.bsub

library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

project_dir <- "/beevol/home/matlinka/timecourse"
work_dir    <- file.path(project_dir, "viennarna")

txdb_path   <- file.path(work_dir, "accessories", "human", "txdb.gencode49.sqlite")
tx_mapping_path <- file.path(work_dir, "output", "predictive_modeling", "precomputed_most_abundant_tx.rds")

# Load TxDb and genome
txdb   <- loadDb(txdb_path)
genome <- BSgenome.Hsapiens.UCSC.hg38

# Load precomputed most abundant transcript mapping (same set as 5'UTR run)
precomputed_tx <- readRDS(tx_mapping_path)
cat("Loaded most abundant transcript mapping:", nrow(precomputed_tx), "transcripts\n")

# ---------------------------------------------------------------------------
# CDS sequences
# ---------------------------------------------------------------------------
cds_all      <- cdsBy(txdb, by = "tx", use.names = TRUE)
cds_seqs_all <- extractTranscriptSeqs(genome, cds_all)
cat("Total CDS sequences in TxDb:", length(cds_seqs_all), "\n")

cds_names_clean <- sub("\\..*", "", names(cds_seqs_all))
keep_idx        <- cds_names_clean %in% precomputed_tx$transcript_id_clean
cds_seqs        <- cds_seqs_all[keep_idx]
names(cds_seqs) <- cds_names_clean[keep_idx]
cat("CDS sequences matching transcript set:", length(cds_seqs), "\n")

# Remove sequences shorter than 10 nt
cds_seqs <- cds_seqs[width(cds_seqs) >= 10]
cat("After length filter (>= 10 nt):", length(cds_seqs), "\n")

cat("\nCDS length summary:\n")
print(summary(as.integer(width(cds_seqs))))

# Write CDS FASTA — all sequences (for RNAplfold)
cds_fa_path <- file.path(work_dir, "cds_sequences.fa")
fasta_lines <- character(length(cds_seqs) * 2)
for (i in seq_along(cds_seqs)) {
  fasta_lines[(i - 1) * 2 + 1] <- paste0(">", names(cds_seqs)[i])
  fasta_lines[(i - 1) * 2 + 2] <- as.character(cds_seqs[[i]])
}
writeLines(fasta_lines, cds_fa_path)
cat("Wrote CDS FASTA:", cds_fa_path, "\n")

# Write CDS FASTA filtered to <= 10,000 nt (for RNAfold — excludes extreme outliers)
cds_seqs_rnafold <- cds_seqs[width(cds_seqs) <= 10000]
n_excluded <- length(cds_seqs) - length(cds_seqs_rnafold)
cat("CDS sequences > 10,000 nt excluded from RNAfold FASTA:", n_excluded, "\n")

cds_rnafold_fa_path <- file.path(work_dir, "cds_sequences_rnafold.fa")
fasta_lines_rf <- character(length(cds_seqs_rnafold) * 2)
for (i in seq_along(cds_seqs_rnafold)) {
  fasta_lines_rf[(i - 1) * 2 + 1] <- paste0(">", names(cds_seqs_rnafold)[i])
  fasta_lines_rf[(i - 1) * 2 + 2] <- as.character(cds_seqs_rnafold[[i]])
}
writeLines(fasta_lines_rf, cds_rnafold_fa_path)
cat("Wrote CDS RNAfold FASTA:", cds_rnafold_fa_path,
    "(", length(cds_seqs_rnafold), "sequences)\n")

# Save metadata
saveRDS(
  data.frame(
    transcript_id_clean = names(cds_seqs),
    cds_length          = as.integer(width(cds_seqs))
  ),
  file.path(work_dir, "cds_seq_metadata.rds")
)

# ---------------------------------------------------------------------------
# 3'UTR sequences
# ---------------------------------------------------------------------------
utr3_all      <- threeUTRsByTranscript(txdb, use.names = TRUE)
utr3_seqs_all <- extractTranscriptSeqs(genome, utr3_all)
cat("\nTotal 3'UTR sequences in TxDb:", length(utr3_seqs_all), "\n")

utr3_names_clean <- sub("\\..*", "", names(utr3_seqs_all))
keep_idx_utr3    <- utr3_names_clean %in% precomputed_tx$transcript_id_clean
utr3_seqs        <- utr3_seqs_all[keep_idx_utr3]
names(utr3_seqs) <- utr3_names_clean[keep_idx_utr3]
cat("3'UTR sequences matching transcript set:", length(utr3_seqs), "\n")

n_no_utr3 <- sum(!precomputed_tx$transcript_id_clean %in% names(utr3_seqs))
cat("Transcripts without annotated 3'UTR:", n_no_utr3,
    sprintf("(%.1f%%)\n", 100 * n_no_utr3 / nrow(precomputed_tx)))

# Remove sequences shorter than 10 nt
utr3_seqs <- utr3_seqs[width(utr3_seqs) >= 10]
cat("After length filter (>= 10 nt):", length(utr3_seqs), "\n")

cat("\n3'UTR length summary:\n")
print(summary(as.integer(width(utr3_seqs))))

# Write 3'UTR FASTA
utr3_fa_path <- file.path(work_dir, "utr3_sequences.fa")
fasta_lines_utr3 <- character(length(utr3_seqs) * 2)
for (i in seq_along(utr3_seqs)) {
  fasta_lines_utr3[(i - 1) * 2 + 1] <- paste0(">", names(utr3_seqs)[i])
  fasta_lines_utr3[(i - 1) * 2 + 2] <- as.character(utr3_seqs[[i]])
}
writeLines(fasta_lines_utr3, utr3_fa_path)
cat("Wrote 3'UTR FASTA:", utr3_fa_path, "\n")

# Save metadata
saveRDS(
  data.frame(
    transcript_id_clean = names(utr3_seqs),
    utr3_length         = as.integer(width(utr3_seqs))
  ),
  file.path(work_dir, "utr3_seq_metadata.rds")
)

cat("\nDone.\n")
