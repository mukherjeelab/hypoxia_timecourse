# 02e_validate_rnaplfold_params.R
# Prove the RNAplfold parameters used for the MCF7-SIX1 transcripts are the ones the MDA-MB-231
# output was produced with, by regenerating stored 231 files and requiring a BYTE-IDENTICAL match.
#
# Why this is needed rather than reading the scripts: the committed 231 and Teleman runners call
# "RNAplfold -W 150 -L 100" with no -u, which on ViennaRNA 2.7.2 writes no _lunp file at all, so
# the scripts as committed cannot be what produced the stored data. The parameters therefore have
# to be recovered from the output and then verified, not taken from the source.
#
# Only files written BEFORE today are sampled, so this tests original 231 output rather than
# files this work has just added to the same directories.
#
# Usage: Rscript cluster_scripts/viennarna/02e_validate_rnaplfold_params.R [n_per_region]

suppressPackageStartupMessages({
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings); library(dplyr); library(here)
})

args <- commandArgs(trailingOnly = TRUE)
N    <- if (length(args) >= 1) as.integer(args[1]) else 60L

# The parameters under test, PER REGION. They are not the same for all three, which is the whole
# reason this script exists: the 5'UTR output was generated in an earlier run with a smaller
# window, and running -W 150 -L 100 over it reproduced only 12 of 60 files - the 12 being short
# 5'UTRs where the window exceeds the sequence and so has no effect. On the rest, values differed
# by up to 0.94 on a 0-1 scale. A single global setting here would have passed the CDS and 3'UTR
# checks and silently mis-generated every 5'UTR.
FLAGS <- c(utr5 = "-W 80 -L 40 -u 30",
           cds  = "-W 150 -L 100 -u 30",
           utr3 = "-W 150 -L 100 -u 30")

txdb <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
plf  <- here("accessories", "plfold_output", "rnaplfold_output")
tx231 <- readRDS(here("output", "predictive_modeling",
                      "precomputed_most_abundant_tx.rds"))$transcript_id_clean
scratch <- file.path(tempdir(), "plfval"); dir.create(scratch, showWarnings = FALSE)

regions <- list(
  utr5 = list(dir = "rnaplfold_output",      grl = fiveUTRsByTranscript(txdb, use.names = TRUE)),
  cds  = list(dir = "cds_rnaplfold_output",  grl = cdsBy(txdb, by = "tx", use.names = TRUE)),
  utr3 = list(dir = "utr3_rnaplfold_output", grl = threeUTRsByTranscript(txdb, use.names = TRUE))
)

cutoff <- as.POSIXct(format(Sys.Date()))    # midnight today
summary_rows <- list()

for (rg in names(regions)) {
  spec <- regions[[rg]]
  d    <- file.path(plf, spec$dir)
  # list.files, not Sys.glob: these directories exceed ARG_MAX for a shell glob.
  f    <- list.files(d, pattern = "_lunp$", full.names = TRUE)
  info <- file.info(f)
  orig <- f[info$mtime < cutoff]                       # written before today = original 231 run
  ids  <- sub("_lunp$", "", basename(orig))
  cand <- ids[ids %in% tx231]
  stopifnot("no pre-existing 231 files to validate against" = length(cand) > 0)

  set.seed(9)
  pick <- sample(cand, min(N, length(cand)))

  seqs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, spec$grl)
  names(seqs) <- sub("\\..*", "", names(seqs))
  seqs <- seqs[names(seqs) %in% pick]
  stopifnot("sampled transcripts are not all extractable" = length(seqs) == length(pick))

  wd <- file.path(scratch, rg); dir.create(wd, showWarnings = FALSE)
  fa <- file.path(wd, "in.fa"); writeXStringSet(seqs, fa)
  old <- setwd(wd)
  system2("RNAplfold", args = strsplit(FLAGS[[rg]], " ")[[1]], stdin = "in.fa",
          stdout = FALSE, stderr = FALSE)
  setwd(old)

  cmp <- vapply(names(seqs), function(tx) {
    a <- file.path(d, paste0(tx, "_lunp")); b <- file.path(wd, paste0(tx, "_lunp"))
    if (!file.exists(b)) return(NA_character_)
    ra <- readBin(a, "raw", file.size(a)); rb <- readBin(b, "raw", file.size(b))
    if (identical(ra, rb)) "identical" else "DIFFERS"
  }, character(1))

  n_id <- sum(cmp == "identical", na.rm = TRUE); n_df <- sum(cmp == "DIFFERS", na.rm = TRUE)
  n_na <- sum(is.na(cmp))
  cat(sprintf("%-5s tested %3d | identical %3d | DIFFERS %3d | not produced %3d | pool of originals %5d\n",
              rg, length(cmp), n_id, n_df, n_na, length(cand)))
  if (n_df > 0) {
    bad <- names(cmp)[which(cmp == "DIFFERS")]
    cat("  differing:", paste(head(bad, 8), collapse = ", "), "\n")
    for (tx in head(bad, 3)) {
      rd <- function(p) { L <- readLines(p); m <- do.call(rbind, lapply(strsplit(L[-(1:2)], "\t"),
              function(z) suppressWarnings(as.numeric(z)))); m }
      A <- rd(file.path(d, paste0(tx, "_lunp"))); B <- rd(file.path(wd, paste0(tx, "_lunp")))
      cat(sprintf("    %s dims %s vs %s | max|dev| %.3e\n", tx,
                  paste(dim(A), collapse = "x"), paste(dim(B), collapse = "x"),
                  if (identical(dim(A), dim(B))) max(abs(A - B), na.rm = TRUE) else NA_real_))
    }
  }
  summary_rows[[rg]] <- data.frame(region = rg, tested = length(cmp), identical = n_id,
                                   differs = n_df, missing = n_na)
}

res <- do.call(rbind, summary_rows)
cat("\n"); print(res, row.names = FALSE)
stopifnot(
  "a regenerated file is not byte-identical to the stored MDA-MB-231 output - the parameters differ" =
    all(res$differs == 0),
  "RNAplfold did not produce a file for some sampled transcript" = all(res$missing == 0),
  "every region must actually have been tested" = all(res$tested > 0)
)
cat("\nVALIDATED byte-for-byte against stored MDA-MB-231 output, per region:\n")
for (rg in names(FLAGS)) cat(sprintf("  %-5s %s\n", rg, FLAGS[[rg]]))
cat(sprintf("over %d transcripts total. The MCF7-SIX1 files are generated with these same\n",
            sum(res$tested)))
cat("per-region parameters, so accessibility is comparable across the two cell lines within a\n")
cat("region. It is NOT comparable ACROSS regions: the 5'UTR window is smaller.\n")
