# 01_generate_fasta.R
# Generate FASTA files of 5'UTR, CDS, and 3'UTR sequences for G4mer rG4 scoring.
#
# Adapted from cluster_scripts/viennarna/01b_generate_fasta_cds_utr3.R so that the
# G4mer transcript universe is identical to the ViennaRNA structure-feature universe:
# both subset to output/predictive_modeling/precomputed_most_abundant_tx.rds.
#
# Local:   Rscript g4mer/01_generate_fasta.R
# Cluster: Rscript 01_generate_fasta.R --workdir /beevol/home/matlinka/timecourse/g4mer

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
})

# ---------------------------------------------------------------------------
# Paths — default to repo-relative via here(), override with --workdir on cluster
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
workdir_arg <- if ("--workdir" %in% args) args[which(args == "--workdir") + 1] else NA

if (is.na(workdir_arg)) {
  suppressPackageStartupMessages(library(here))
  txdb_path       <- here("accessories", "human", "txdb.gencode49.sqlite")
  tx_mapping_path <- here("output", "predictive_modeling", "precomputed_most_abundant_tx.rds")
  out_dir         <- here("g4mer", "fasta")
} else {
  txdb_path       <- file.path(workdir_arg, "accessories", "human", "txdb.gencode49.sqlite")
  tx_mapping_path <- file.path(workdir_arg, "output", "predictive_modeling",
                               "precomputed_most_abundant_tx.rds")
  out_dir         <- file.path(workdir_arg, "fasta")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot("TxDb not found — run 00_build_txdb_v49.R first" = file.exists(txdb_path))
stopifnot("precomputed_most_abundant_tx.rds not found" = file.exists(tx_mapping_path))

txdb   <- loadDb(txdb_path)
genome <- BSgenome.Hsapiens.UCSC.hg38

precomputed_tx <- readRDS(tx_mapping_path)
cat("Loaded most abundant transcript mapping:", nrow(precomputed_tx), "transcripts\n")

# Minimum length: G4mer uses 6-mer tokens, so anything under ~10 nt is uninformative.
MIN_LEN <- 10L

# ---------------------------------------------------------------------------
# Helper: extract region, subset to transcript universe, write FASTA + metadata
# ---------------------------------------------------------------------------
gc_percent <- function(dna) {
  100 * as.numeric(letterFrequency(dna, letters = "GC", as.prob = TRUE))
}

# G-richness of the mRNA-sense strand. rG4s form from G-runs, so this is a more
# targeted covariate than GC alone when interpreting G4mer scores downstream.
g_percent <- function(dna) {
  100 * as.numeric(letterFrequency(dna, letters = "G", as.prob = TRUE))
}

write_region <- function(region_grl, region_name) {
  # extractTranscriptSeqs is strand-corrected and exon-spliced — gives mRNA-sense
  # sequence directly. Do not substitute bedtools getfasta on genomic intervals.
  seqs_all <- extractTranscriptSeqs(genome, region_grl)
  cat("\n[", region_name, "] total sequences in TxDb:", length(seqs_all), "\n")

  names_clean <- sub("\\..*", "", names(seqs_all))
  keep_idx    <- names_clean %in% precomputed_tx$transcript_id_clean
  seqs        <- seqs_all[keep_idx]
  names(seqs) <- names_clean[keep_idx]
  cat("[", region_name, "] matching transcript set:", length(seqs), "\n")

  n_missing <- sum(!precomputed_tx$transcript_id_clean %in% names(seqs))
  cat("[", region_name, "] transcripts with no annotated region:", n_missing,
      sprintf("(%.1f%%)\n", 100 * n_missing / nrow(precomputed_tx)))

  seqs <- seqs[width(seqs) >= MIN_LEN]
  cat("[", region_name, "] after length filter (>=", MIN_LEN, "nt):", length(seqs), "\n")

  cat("[", region_name, "] length summary:\n")
  print(summary(as.integer(width(seqs))))

  fa_path <- file.path(out_dir, paste0(region_name, "_sequences.fa"))
  writeXStringSet(seqs, fa_path)
  cat("[", region_name, "] wrote FASTA:", fa_path, "\n")

  # Metadata carries the covariates the association analysis needs (length, GC,
  # G-content) so 01u_ never has to re-extract sequences.
  meta <- data.frame(
    transcript_id_clean = names(seqs),
    region              = region_name,
    region_length       = as.integer(width(seqs)),
    gc_pct              = gc_percent(seqs),
    g_pct               = g_percent(seqs),
    stringsAsFactors    = FALSE
  )
  meta <- merge(
    meta,
    precomputed_tx[, c("transcript_id_clean", "gene_id_clean", "symbol")],
    by = "transcript_id_clean", all.x = TRUE
  )

  # Validation: GC and G percentages are bounded by definition
  stopifnot("GC% out of range" = all(meta$gc_pct >= 0 & meta$gc_pct <= 100))
  stopifnot("G% out of range"  = all(meta$g_pct  >= 0 & meta$g_pct  <= 100))
  stopifnot("region_length below filter" = all(meta$region_length >= MIN_LEN))

  meta_path <- file.path(out_dir, paste0(region_name, "_seq_metadata.tsv"))
  write.table(meta, meta_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[", region_name, "] wrote metadata:", meta_path, "\n")

  invisible(meta)
}

# ---------------------------------------------------------------------------
# Extract all three regions
# ---------------------------------------------------------------------------
meta_utr5 <- write_region(fiveUTRsByTranscript(txdb, use.names = TRUE), "utr5")
meta_cds  <- write_region(cdsBy(txdb, by = "tx", use.names = TRUE),     "cds")
meta_utr3 <- write_region(threeUTRsByTranscript(txdb, use.names = TRUE), "utr3")

# ---------------------------------------------------------------------------
# Biological spot-check: ACTB CDS GC% should be 50-70% (matches 01_ convention)
# ---------------------------------------------------------------------------
actb <- meta_cds[which(meta_cds$symbol == "ACTB"), ]
if (nrow(actb) == 1) {
  cat("\nSpot-check ACTB CDS GC%:", round(actb$gc_pct, 1), "\n")
  stopifnot("ACTB CDS GC% outside expected 50-70% range" =
              actb$gc_pct >= 50 && actb$gc_pct <= 70)
} else {
  warning("ACTB not found in CDS metadata — spot-check skipped")
}

cat("\nDone. FASTA + metadata written to:", out_dir, "\n")
