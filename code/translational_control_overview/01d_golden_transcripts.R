# Golden-transcript snapshot: 20 transcripts, every modelled feature, committed to the repo.
#
# Range checks pass on a scrambled join; family-level statistics pass on a value that has
# drifted. A frozen per-transcript snapshot is the check that catches "this number is not what
# it was last week" regardless of where the change came from.
#
# The file is REGENERATED ON PURPOSE, never automatically. When a feature is legitimately
# corrected - as seven were by the audit - re-run this script and the diff in the commit shows
# exactly which transcripts and columns moved. A silent change fails 01_'s assertion instead.
#
# Transcripts are chosen deterministically and span the awkward cases, not just easy ones:
# the named biological spot-check genes, transcripts with NA structure (the F1 group),
# truncated CDS, and non-AUG starts.
#
# Usage: Rscript code/translational_control_overview/01d_golden_transcripts.R
suppressMessages({library(tidyverse); library(here)})

out_dir <- here("output", "translational_control_overview")
features <- readRDS(file.path(out_dir, "feature_matrix_tco.rds"))
ann      <- read_csv(file.path(out_dir, "feature_annotation.csv"), show_col_types = FALSE)

modelled <- ann %>% filter(block %in% c("intrinsic", "external"), !deprecated) %>% pull(feature)
modelled <- modelled[map_lgl(modelled, ~ is.numeric(features[[.x]]))]

pick_first <- function(idx, n) head(sort(features$transcript_id_clean[idx]), n)

chosen <- unique(c(
  # named genes the biological spot-checks depend on
  features$transcript_id_clean[features$symbol %in%
    c("ACTB", "COL1A1", "RPL13A", "RPS20", "EIF3D", "EIF3E", "DHX29")],
  # the F1 group: structure restored to NA
  pick_first(is.na(features$tco_struct_accessibility_cds_mean), 4),
  # truncated CDS (cds_end_NF) and non-AUG starts
  pick_first(!features$cds_has_stop, 3),
  pick_first(!features$cds_has_start & features$cds_has_stop, 3),
  # plain transcripts, for a baseline
  pick_first(features$cds_complete & !is.na(features$tco_struct_accessibility_cds_mean), 3)
))
chosen <- head(chosen, 20)

golden <- features %>%
  filter(transcript_id_clean %in% chosen) %>%
  dplyr::select(transcript_id_clean, symbol, all_of(modelled)) %>%
  arrange(transcript_id_clean) %>%
  pivot_longer(all_of(modelled), names_to = "feature", values_to = "value")


# ---------------------------------------------------------------------------------------
# Independent verification of the snapshot
# ---------------------------------------------------------------------------------------
# A snapshot on its own is a REGRESSION detector: it freezes whatever the build produced.
# Had this file been made before the audit it would have frozen the 1,689 zero-filled
# structure values and asserted them forever. So each value carries how far it has actually
# been checked:
#
#   recomputed - independently derived here from GENCODE v49 + hg38 + the published reference
#                tables, by code that does not read the builder's intermediates
#   traced     - matched by key to the external file it came from; provenance verified, the
#                underlying measurement is not
#   unverified - neither is possible from here (computed off-repo, or needs a tool we would
#                simply be re-running rather than independently reproducing)
#
# Caveat worth keeping in view: a recomputation written by the same author shares that
# author's assumptions. It proves internal consistency, not truth. Independent re-derivation
# by someone else is what catches a shared misreading.
suppressMessages({library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)
                  library(Biostrings)})

txdb   <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
genome <- BSgenome.Hsapiens.UCSC.hg38
ref_dir <- here("accessories", "collaborator_reference")

grab <- function(fn) {
  g <- fn(txdb, use.names = TRUE)
  keep <- names(g)[sub("\\..*", "", names(g)) %in% chosen]
  sq <- extractTranscriptSeqs(genome, g[keep])
  setNames(as.character(sq), sub("\\..*", "", names(sq)))
}
cds_s  <- grab(function(db, ...) cdsBy(db, by = "tx", ...))
utr5_s <- grab(fiveUTRsByTranscript)
utr3_s <- grab(threeUTRsByTranscript)

gc_pct <- function(x) ifelse(is.na(x) | nchar(x) == 0, NA_real_,
                             100 * (str_count(x, "[GCgc]") / nchar(x)))
stops <- c("TAA", "TAG", "TGA")
next_stop <- function(u) {
  if (is.na(u) || nchar(u) < 3) return(NA_integer_)
  n <- nchar(u) %/% 3L
  cod <- substring(u, seq(1, n * 3L, 3L), seq(3, n * 3L, 3L))
  w <- which(cod %in% stops); if (!length(w)) NA_integer_ else w[1]
}
charge_max <- function(aa, w = 30L) {
  ch <- strsplit(aa, "", fixed = TRUE)[[1]]; if (!length(ch)) return(NA_real_)
  v <- integer(length(ch)); v[ch %in% c("K","R")] <- 1L; v[ch %in% c("D","E")] <- -1L
  if (length(v) <= w) return(as.numeric(sum(v)))
  cs <- cumsum(c(0L, v)); as.numeric(max(cs[(w+1L):length(cs)] - cs[seq_len(length(cs)-w)]))
}
codon_run <- function(s, codon) {
  n <- nchar(s) %/% 3L; if (n < 1L) return(0L)
  r <- rle(substring(s, seq(1, n*3L, 3L), seq(3, n*3L, 3L)) == codon)
  if (!any(r$values)) 0L else max(r$lengths[r$values])
}

cridge <- read_tsv(file.path(ref_dir, "cridge2018_termination_contexts.tsv"), show_col_types = FALSE)
noderer <- read_tsv(file.path(ref_dir, "noderer_2014_tis_efficiency.tsv"), skip = 1,
                    show_col_types = FALSE) %>% mutate(context = chartr("U", "T", sequence))

# A named character vector errors on a missing name rather than returning NULL, and 134
# transcripts have no annotated 3'UTR at all - which is exactly why the snapshot includes
# awkward transcripts.
getseq <- function(v, tx) if (tx %in% names(v)) str_to_upper(unname(v[[tx]])) else NA_character_

recomp <- map_dfr(sort(unique(golden$transcript_id_clean)), function(tx) {
  cd <- getseq(cds_s, tx); u5 <- getseq(utr5_s, tx); u3 <- getseq(utr3_s, tx)
  if (is.na(cd)) return(tibble())
  n <- nchar(cd); ncod <- n %/% 3L
  stop_cod <- substr(cd, (ncod-1)*3L+1L, ncod*3L)
  has_stop <- stop_cod %in% stops
  s_int <- substr(cd, if (substr(cd,1,3) == "ATG") 4L else 1L, if (has_stop) (ncod-1)*3L else ncod*3L)
  aa <- if (nchar(s_int) >= 3 && nchar(s_int) %% 3 == 0)
    as.character(Biostrings::translate(DNAStringSet(s_int), if.fuzzy.codon = "solve")) else NA_character_
  p4 <- if (!is.na(u3) && nchar(u3) >= 1) substr(u3, 1, 1) else NA_character_
  ctx6 <- if (has_stop && !is.na(u3) && nchar(u3) >= 3) paste0(stop_cod, substr(u3,1,3)) else NA_character_
  ctx9 <- if (has_stop && !is.na(u3) && nchar(u3) >= 6) paste0(stop_cod, substr(u3,1,6)) else NA_character_
  nod_ctx <- if (!is.na(u5) && nchar(u5) >= 6 && substr(cd,1,3) == "ATG")
    paste0(substr(u5, nchar(u5)-5, nchar(u5)), substr(cd,1,5)) else NA_character_
  k <- if (is.na(ctx6)) NA_real_ else cridge$readthrough_measure[match(ctx6,
         cridge$termination_context_plus1_to_plus6)]
  gc3v <- { kk <- seq.int(3L, ncod*3L, 3L); 100*sum(substr(rep.int(cd,ncod), kk, kk) %in% c("G","C"))/ncod }
  ni <- nchar(s_int) %/% 3L
  gc3i <- if (ni > 0) { kk <- seq.int(3L, ni*3L, 3L)
    100*sum(substr(rep.int(s_int,ni), kk, kk) %in% c("G","C"))/ni } else NA_real_
  tibble(transcript_id_clean = tx, feature = c(
    "utr5_length","cds_length","utr3_length",
    "utr5_gc","cds_gc","utr3_gc","gc3","gc3_internal",
    "term_stop_TAA","term_stop_TAG","term_stop_TGA",
    "term_plus4_A","term_plus4_C","term_plus4_G","term_plus4_T",
    "term_window_gc","term_interval_codons","term_readthrough_cridge",
    "proline_fraction","ppp_motif_density","max_net_charge_30aa",
    "max_consecutive_aaa","max_consecutive_aag","noderer_tis_efficiency"),
    recomputed = c(
      if (is.na(u5)) 0 else nchar(u5), n, if (is.na(u3)) 0 else nchar(u3),
      gc_pct(u5), gc_pct(cd), gc_pct(u3), gc3v, if (has_stop) gc3i else NA_real_,
      as.integer(stop_cod=="TAA"), as.integer(stop_cod=="TAG"), as.integer(stop_cod=="TGA"),
      as.integer(p4=="A"), as.integer(p4=="C"), as.integer(p4=="G"), as.integer(p4=="T"),
      if (is.na(ctx9)) NA_real_ else gc_pct(ctx9),
      if (is.na(u3)) NA_integer_ else next_stop(u3), if (length(k)) k else NA_real_,
      if (is.na(aa)) NA_real_ else as.numeric(letterFrequency(AAStringSet(aa),"P",as.prob=TRUE)),
      if (is.na(aa)) NA_real_ else 100*vcountPattern("PPP", AAStringSet(aa))/nchar(aa),
      if (is.na(aa)) NA_real_ else charge_max(aa),
      codon_run(s_int,"AAA"), codon_run(s_int,"AAG"),
      if (is.na(nod_ctx)) NA_real_ else noderer$efficiency[match(nod_ctx, noderer$context)]))
})

# Mask the codon-derived recomputations on transcripts the builder masks (no terminal stop).
trunc_tx <- features$transcript_id_clean[!features$cds_has_stop]
recomp <- recomp %>%
  mutate(recomputed = if_else(transcript_id_clean %in% trunc_tx &
                                feature %in% c("gc3_internal"), NA_real_, recomputed))

cmp <- golden %>% inner_join(recomp, by = c("transcript_id_clean", "feature")) %>%
  mutate(dev = abs(value - recomputed))

# Value disagreement, where both sides have a number.
bad <- cmp %>% filter(!is.na(dev), dev > 1e-6)
cat("\nIndependent recomputation:", nrow(cmp), "values across",
    n_distinct(cmp$feature), "features\n")
cat("  numeric disagreements >1e-6:", nrow(bad), "\n")
if (nrow(bad)) print(as.data.frame(bad %>% arrange(desc(dev)) %>% head(10)))

# The two NA directions mean opposite things, so they are counted separately.
#   builder has a value, recomputation does not -> the builder may be scoring something that
#       cannot be computed (a truncated CDS, a missing 3'UTR). This is the dangerous direction.
#   builder is NA, recomputation has a value    -> the builder is masking conservatively.
#       Expected: codon and peptide features are deliberately masked wherever cds_has_stop is
#       FALSE, and utr3-derived columns are NA where GENCODE annotates no 3'UTR.
fabricated <- cmp %>% filter(!is.na(value), is.na(recomputed))
masked     <- cmp %>% filter(is.na(value), !is.na(recomputed))
cat("  builder has a value where recomputation cannot:", nrow(fabricated), "\n")
if (nrow(fabricated)) print(as.data.frame(head(fabricated, 8)))
cat("  builder masks where recomputation can (expected: truncated CDS / no 3'UTR):",
    nrow(masked), "across", n_distinct(masked$transcript_id_clean), "transcripts\n")

no_stop <- features$transcript_id_clean[!features$cds_has_stop]
unexplained <- masked %>% filter(!transcript_id_clean %in% no_stop,
                                 !str_detect(feature, "utr3|term_"))
cat("  of those, NOT explained by a truncated CDS or a 3'UTR feature:", nrow(unexplained), "\n")
if (nrow(unexplained)) print(as.data.frame(head(unexplained, 8)))

stopifnot(
  "recomputation disagrees with the stored feature" = nrow(bad) == 0,
  "the builder scores a feature the primary data cannot support" = nrow(fabricated) == 0,
  "a masked value is not explained by CDS truncation or a missing 3'UTR" = nrow(unexplained) == 0
)

# Traced: provenance confirmed against the source file, value not independently derived.
traced <- c(grep("^g4mer_", unique(golden$feature), value = TRUE),
            "hia2026_dhx29_occupancy_tx")
golden <- golden %>%
  mutate(status = case_when(feature %in% recomp$feature ~ "recomputed",
                            feature %in% traced          ~ "traced",
                            TRUE                         ~ "unverified"))
cat("\nVerification status of the snapshot:\n")
print(golden %>% distinct(feature, status) %>% count(status) %>% as.data.frame())

f <- here("code", "translational_control_overview", "golden_transcripts.csv")
write_csv(golden, f)
cat("Wrote", f, "-", n_distinct(golden$transcript_id_clean), "transcripts x",
    n_distinct(golden$feature), "features =", nrow(golden), "values\n")
cat("  recomputed:", sum(golden$status == "recomputed"),
    " traced:", sum(golden$status == "traced"),
    " unverified:", sum(golden$status == "unverified"), "(values)\n")
cat("Transcripts:", paste(sort(unique(golden$transcript_id_clean)), collapse = ", "), "\n")
