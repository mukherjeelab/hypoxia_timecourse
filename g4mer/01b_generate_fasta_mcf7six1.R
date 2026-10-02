# 01b_generate_fasta_mcf7six1.R
# Region FASTAs for the MCF7-SIX1 dominant transcripts that do NOT already have a G4mer score.
#
# Only the missing transcripts are scored. A G4mer score is a property of the sequence, and
# regenerating the 5'UTR summary reproduced the stored values BIT-IDENTICALLY, so a transcript
# dominant in both cell lines needs no second inference. That takes the job from ~11,200
# transcripts per region to ~5,600.
#
# Usage: Rscript g4mer/01b_generate_fasta_mcf7six1.R [cell]
#   cell = mcf7six1 (default) | hela. HeLa uses precomputed_teleman_tx.rds (the Roiuk-dominant
#   transcript per gene). A transcript already scored in ANY existing summary is reused.

suppressPackageStartupMessages({
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings); library(dplyr); library(readr); library(here)
})
MIN_LEN <- 10L                       # same filter as 01_generate_fasta.R: 6-mer tokenisation
out_dir <- here("g4mer", "fasta"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
g4_dir  <- here("output", "g4mer")

CELL   <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "mcf7six1" }
TX_RDS <- c(mcf7six1 = "precomputed_mcf7six1_tx.rds", hela = "precomputed_teleman_tx.rds",
            dhx29 = "precomputed_dhx29_tx.rds", mane = "precomputed_mane_tx.rds",
           mdamb231 = "precomputed_mdamb231_tx.rds")
stopifnot("unknown cell set" = CELL %in% names(TX_RDS))
tx_map <- readRDS(here("output", "predictive_modeling", TX_RDS[[CELL]]))
want   <- unique(tx_map$transcript_id_clean)
cat(CELL, "transcripts:", length(want), "\n")

txdb   <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
genome <- BSgenome.Hsapiens.UCSC.hg38
REG <- list(utr5 = function() fiveUTRsByTranscript(txdb, use.names = TRUE),
            cds  = function() cdsBy(txdb, by = "tx", use.names = TRUE),
            utr3 = function() threeUTRsByTranscript(txdb, use.names = TRUE))

for (rg in names(REG)) {
  # Subset the ranges BEFORE extracting: extracting all 233,995 transcripts exhausts memory on
  # this machine and is most of the work thrown away.
  grl <- REG[[rg]]()
  names(grl) <- sub("\\..*", "", names(grl))
  grl <- grl[names(grl) %in% want]
  seqs <- extractTranscriptSeqs(genome, grl)
  names(seqs) <- sub("\\..*", "", names(seqs))
  seqs <- seqs[width(seqs) >= MIN_LEN]
  rm(grl); invisible(gc(FALSE))

  # Every existing summary for this region: the 231 base plus any other cell line's merged one.
  fs <- c(file.path(g4_dir, paste0("g4mer_summary_", rg, ".tsv")),
          list.files(g4_dir, pattern = paste0("^g4mer_summary_", rg, "_.*\\.tsv$"), full.names = TRUE))
  have <- unique(unlist(lapply(fs[file.exists(fs)],
                               function(f) read_tsv(f, show_col_types = FALSE)$transcript_id_clean)))
  need <- setdiff(names(seqs), have)
  out  <- seqs[names(seqs) %in% need]
  p <- file.path(out_dir, paste0(CELL, "_", rg, "_sequences.fa"))
  writeXStringSet(out, p)
  cat(sprintf("%-5s extractable %6d | already scored %6d | to score %6d | %9.0f nt -> %s\n",
              rg, length(seqs), length(seqs) - length(need), length(out),
              sum(as.numeric(width(out))), basename(p)))
  rm(seqs, out); invisible(gc(FALSE))
}
cat("\nNext: g4mer/02_run_g4mer.py per region (stride 10 for utr5, 20 for cds/utr3)\n")
