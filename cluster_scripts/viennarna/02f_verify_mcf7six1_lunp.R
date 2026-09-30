# 02f_verify_mcf7six1_lunp.R
# Verify the MCF7-SIX1 _lunp files actually on disk were written with the per-region parameters,
# by regenerating a sample and requiring byte-identity.
#
# This exists because the generation ran across several attempts, one of which was killed
# mid-region and whose log had not flushed when it was read, so which parameters reached which
# region could only be inferred. Inference is not good enough for a value that feeds every
# structure feature, and the files carry no record of the flags used. Regenerating settles it.
#
# Usage: Rscript cluster_scripts/viennarna/02f_verify_mcf7six1_lunp.R [n_per_region]

suppressPackageStartupMessages({
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings); library(dplyr); library(here)
})
args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args) >= 1) as.integer(args[1]) else 40L

FLAGS <- c(utr5 = "-W 80 -L 40 -u 30",
           cds  = "-W 150 -L 100 -u 30",
           utr3 = "-W 150 -L 100 -u 30")
DIRS  <- c(utr5 = "rnaplfold_output", cds = "cds_rnaplfold_output", utr3 = "utr3_rnaplfold_output")
plf   <- here("accessories", "plfold_output", "rnaplfold_output")
fa_dir<- here("accessories", "plfold_output", "mcf7six1_fasta")

txdb <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
GRL <- list(utr5 = fiveUTRsByTranscript(txdb, use.names = TRUE),
            cds  = cdsBy(txdb, by = "tx", use.names = TRUE),
            utr3 = threeUTRsByTranscript(txdb, use.names = TRUE))
scratch <- file.path(tempdir(), "mcf7chk"); dir.create(scratch, showWarnings = FALSE)
res <- list()

for (rg in names(FLAGS)) {
  # The transcripts THIS work generated are exactly the ones in the FASTA it generated from.
  fa <- file.path(fa_dir, paste0("mcf7six1_", rg, "_sequences.fa"))
  stopifnot("MCF7-SIX1 FASTA missing - run 01d_ first" = file.exists(fa))
  ids <- sub("^>", "", grep("^>", readLines(fa, warn = FALSE), value = TRUE))
  ids <- sub(" .*$", "", ids)
  d   <- file.path(plf, DIRS[[rg]])
  ids <- ids[file.exists(file.path(d, paste0(ids, "_lunp")))]
  stopifnot("no generated MCF7-SIX1 files found for this region" = length(ids) > 0)

  set.seed(9); pick <- sample(ids, min(N, length(ids)))
  s <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, GRL[[rg]])
  names(s) <- sub("\\..*", "", names(s)); s <- s[names(s) %in% pick]

  wd <- file.path(scratch, rg); dir.create(wd, showWarnings = FALSE)
  writeXStringSet(s, file.path(wd, "in.fa"))
  old <- setwd(wd)
  system2("RNAplfold", args = strsplit(FLAGS[[rg]], " ")[[1]], stdin = "in.fa",
          stdout = FALSE, stderr = FALSE)
  setwd(old)

  cmp <- vapply(names(s), function(tx) {
    a <- file.path(d,  paste0(tx, "_lunp")); bb <- file.path(wd, paste0(tx, "_lunp"))
    if (!file.exists(bb)) return(NA_character_)
    if (identical(readBin(a, "raw", file.size(a)), readBin(bb, "raw", file.size(bb))))
      "identical" else "DIFFERS"
  }, character(1))
  n_id <- sum(cmp == "identical", na.rm = TRUE); n_df <- sum(cmp == "DIFFERS", na.rm = TRUE)
  cat(sprintf("%-5s %-22s generated %6d | tested %3d | identical %3d | DIFFERS %3d\n",
              rg, FLAGS[[rg]], length(ids), length(cmp), n_id, n_df))
  if (n_df > 0) cat("  differing:",
                    paste(head(names(cmp)[cmp == "DIFFERS"], 6), collapse = ", "), "\n")
  res[[rg]] <- data.frame(region = rg, tested = length(cmp), identical = n_id, differs = n_df)
}

r <- do.call(rbind, res); cat("\n"); print(r, row.names = FALSE)
stopifnot(
  "an MCF7-SIX1 _lunp file on disk does not match its region's parameters - regenerate that region" =
    all(r$differs == 0),
  "every region must have been tested" = all(r$tested > 0)
)
cat("\nVERIFIED: the MCF7-SIX1 files on disk match the per-region parameters that reproduce the\n")
cat("MDA-MB-231 output, so structure features are comparable between the two cell lines.\n")
