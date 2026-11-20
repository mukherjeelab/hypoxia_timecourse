# Length and GC% per region for a set of transcripts (GENCODE v49, hg38).
# Usage: Rscript code/finding_positive_controls/transcript_length_gc.R ENST... ENST...
# Writes/updates output/finding_positive_controls/transcript_length_gc.csv
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38); library(Biostrings)
})
ids <- sub("\\..*", "", commandArgs(trailingOnly = TRUE))
stopifnot(length(ids) > 0)
txdb <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
g <- BSgenome.Hsapiens.UCSC.hg38

region_seqs <- function(grl) {
  grl <- grl[sub("\\..*", "", names(grl)) %in% ids]
  s <- extractTranscriptSeqs(g, grl); names(s) <- sub("\\..*", "", names(s)); s
}
seqs <- list(
  transcript = region_seqs(exonsBy(txdb, "tx", use.names = TRUE)),
  utr5 = region_seqs(fiveUTRsByTranscript(txdb, use.names = TRUE)),
  cds  = region_seqs(cdsBy(txdb, "tx", use.names = TRUE)),
  utr3 = region_seqs(threeUTRsByTranscript(txdb, use.names = TRUE))
)
stopifnot("transcript not in GENCODE v49" = all(ids %in% names(seqs$transcript)))

gc_pct <- function(s) 100 * rowSums(letterFrequency(s, c("G", "C"))) / width(s)
res <- imap_dfr(seqs, function(s, region) tibble(
  transcript_id_clean = names(s), region = region,
  length_nt = width(s), gc_pct = round(gc_pct(s), 1))) %>%
  complete(transcript_id_clean = ids, region = names(seqs))

# checks: regions tile the transcript; GC in range; CDS in frame and starts ATG
w <- res %>% dplyr::select(-gc_pct) %>% pivot_wider(names_from = region, values_from = length_nt)
stopifnot(all(w$transcript == rowSums(w[c("utr5", "cds", "utr3")], na.rm = TRUE)),
          all(res$gc_pct >= 0 & res$gc_pct <= 100, na.rm = TRUE),
          all(width(seqs$cds) %% 3 == 0),
          all(as.character(subseq(seqs$cds, 1, 3)) == "ATG"))

out <- here("output", "finding_positive_controls", "transcript_length_gc.csv")
if (file.exists(out)) res <- read_csv(out, show_col_types = FALSE) %>%
  filter(!transcript_id_clean %in% ids) %>% bind_rows(res)
write_csv(res, out)
print(as.data.frame(res %>% filter(transcript_id_clean %in% ids) %>% arrange(transcript_id_clean)))
