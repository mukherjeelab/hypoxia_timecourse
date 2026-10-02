# 01i_build_mcf7six1_baseline.R
# SUPERSEDED 2026-10-01 for MCF7-SIX1 by 01l_build_celltx_baseline.R mcf7six1 (via 01m_). This
# builder assumed 31 columns in every _lunp file (dropping 650 short 5'UTRs) and carried CDS GC /
# Kozak from 43_, which blanks them for CDS lengths not divisible by 3. Kept for the record.
# Build the MCF7-SIX1 counterpart of feature_matrix_dhx29_kd_feature.rds, so that
# 01_build_feature_matrix.Rmd can produce the MCF7-SIX1 matrix from the SAME code that produces
# the MDA-MB-231 one. The alternative - a second implementation of the feature definitions - is
# how two cell lines end up measuring different quantities under identical column names.
#
# WHAT IS COMPUTED HERE vs CARRIED OVER
#
# Computed on MCF7-SIX1's own dominant transcripts, because these are properties of the
# transcript and the MCF7-SIX1 choice differs from 231's for 34.8% of genes:
#   * struct_* (18)              from the _lunp files, per-region RNAplfold parameters
#   * csc_first75, csc_q1..q4    iCodon positional codon stability
#
# Carried over from the 231 baseline by GENE or SYMBOL, because that is the key the source data
# has - CLIP joins on gene_id_clean, and the DAP5/eIF4E/Teleman columns join on symbol (see
# 01b_feature_extraction_targeted.Rmd). Joining by gene is not an approximation here; it is the
# resolution at which the measurement exists:
#   * clip_* (23), dap5_*/eif4e_*/te_4e_resistance_lfc (6), hia2026_dhx29_kd_rna_lfc
#
# Carried over as LABELS, never as features:
#   * te_lfc, te_padj, te_lfc_bin, indiv_te_* (6), delta_te_* (2)
# These are MDA-MB-231 quantities. They are block = "label" in feature_annotation.csv and are
# asserted out of every model's feature set. They are kept only because the annotation-parity
# gate requires both cell lines to carry the same columns, and because 10_'s provenance check
# uses te_lfc and indiv_te_hypoxia_1hr as the two things the outcome must NOT be. Their presence
# in an MCF7-SIX1 matrix is deliberate and inert.
#
# BOTH COMPUTED BLOCKS CARRY A REPRODUCTION GATE: the same code is first run over the 231
# transcript set and required to reproduce the stored 231 values. Without that, a difference
# between the two cell lines could be this script rather than the biology.
#
# RUN IN STAGES, ONE PROCESS EACH:
#   Rscript 01i_build_mcf7six1_baseline.R struct    # _lunp -> cache   (no sequences held)
#   Rscript 01i_build_mcf7six1_baseline.R csc       # iCodon -> cache  (CDS only)
#   Rscript 01i_build_mcf7six1_baseline.R kozak     # sequence cache   (CDS + 5'UTR)
#   Rscript 01i_build_mcf7six1_baseline.R assemble  # join + parity gate + write
#   Rscript 01i_build_mcf7six1_baseline.R all       # all four in this process (needs headroom)
#
# Staging is not tidiness. Run as one process this was killed twice with no R error and no
# output - the signature of an out-of-memory kill, whose buffered stdout dies with it - because
# it held every CDS and 5'UTR sequence for its whole life while also running iCodon. Each stage
# now holds only what it needs and checkpoints its result, so a kill costs one stage.

suppressPackageStartupMessages({
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38); library(Biostrings)
  library(dplyr); library(tidyr); library(purrr); library(readr); library(here); library(iCodon)
})

STAGES <- c("struct", "csc", "kozak", "assemble")
.a <- commandArgs(trailingOnly = TRUE)
STAGE <- if (length(.a) >= 1) .a[1] else "all"
stopifnot("stage must be one of struct, csc, kozak, assemble, all" =
            STAGE %in% c(STAGES, "all"))
do_stage <- function(x) STAGE == "all" || STAGE == x
cat("stage:", STAGE, "\n"); flush.console()

PD  <- here("output", "predictive_modeling")
plf <- here("accessories", "plfold_output", "rnaplfold_output")
b   <- readRDS(file.path(PD, "feature_matrix_dhx29_kd_feature.rds"))
m6  <- readRDS(file.path(PD, "feature_matrix_mcf7six1_external_stability.rds"))
stopifnot("MCF7-SIX1 matrix is not one transcript per gene" = !any(duplicated(m6$gene_id_clean)))

txdb   <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
genome <- BSgenome.Hsapiens.UCSC.hg38
grab <- function(grl) { s <- extractTranscriptSeqs(genome, grl)
                        names(s) <- sub("\\..*", "", names(s)); s }
# Names only, no sequence extraction: the universe filter needs transcript IDs, and pulling
# every CDS and 5'UTR sequence just to take names() is what made this run out of memory.
ids_of <- function(grl) unique(sub("\\..*", "", names(grl)))
cds_ids  <- ids_of(cdsBy(txdb, by = "tx", use.names = TRUE))
utr5_ids <- ids_of(fiveUTRsByTranscript(txdb, use.names = TRUE))

# 01_ asserts a 5'UTR and a CDS for EVERY transcript in the matrix, so the universe is
# restricted here rather than failing there.
tx <- m6$transcript_id_clean
tx <- tx[tx %in% cds_ids & tx %in% utr5_ids]
cat(sprintf("MCF7-SIX1 transcripts: %d -> %d with both a 5'UTR and a CDS\n",
            nrow(m6), length(tx)))
m6 <- m6 %>% filter(transcript_id_clean %in% tx)
rm(cds_ids, utr5_ids); invisible(gc(verbose = FALSE))
# Extracted lazily, only by the stages that need them, and ONLY for the transcripts in play.
#
# extractTranscriptSeqs() over the whole TxDb builds sequences for all 233,995 transcripts, which
# is what was killing these stages: the process died with no R error and its buffered output lost,
# at a different point each run. Subsetting the GRangesList first cuts the work about tenfold.
#
# The universe is the MCF7-SIX1 transcripts PLUS the 231 baseline's, because the reproduction
# gates score 231 transcripts and would otherwise silently find nothing to compare.
SEQ_TX <- union(tx, b$transcript_id_clean)
cds_all <- NULL; utr5_all <- NULL
grab_subset <- function(grl) {
  names(grl) <- sub("\\..*", "", names(grl))
  grl <- grl[names(grl) %in% SEQ_TX]
  s <- extractTranscriptSeqs(genome, grl)
  names(s) <- sub("\\..*", "", names(s))
  s
}
need_cds  <- function() { if (is.null(cds_all)) {
  cds_all  <<- grab_subset(cdsBy(txdb, by = "tx", use.names = TRUE))
  cat("  extracted", length(cds_all), "CDS sequences\n"); flush.console() }; invisible(NULL) }
need_utr5 <- function() { if (is.null(utr5_all)) {
  utr5_all <<- grab_subset(fiveUTRsByTranscript(txdb, use.names = TRUE))
  cat("  extracted", length(utr5_all), "5'UTR sequences\n"); flush.console() }; invisible(NULL) }

# ---------------------------------------------------------------------------------------------
# 1. struct_* from the _lunp files
# ---------------------------------------------------------------------------------------------
# read_acc/feats are the correct-parse half of 01f_recompute_structure.R, repeated rather than
# sourced because that file is a script with no functions to import. The reproduction gate below
# is what keeps the copy honest: it must reproduce 01f_'s stored 231 output exactly.
# NCOL is 31 for these files: a position index plus the unpaired probabilities for l = 1..30.
# It is read off the file rather than assumed, and asserted, because -u is what sets it.
read_acc <- function(fp) {
  hdr <- readLines(fp, n = 2L, warn = FALSE)
  h   <- sum(grepl("^[[:space:]]*#", hdr))       # matches " #i$..." - it begins with a space
  # scan(), not readLines()+strsplit(). strsplit interns every field in R's global string cache,
  # which is never freed: over ~34,000 files that is millions of permanent CHARSXPs, and the
  # process died partway through the loop at a different index on each attempt. scan() reads
  # doubles straight out of the file and allocates nothing that outlives the call.
  v <- scan(fp, what = double(), skip = h, na.strings = "NA", quiet = TRUE)
  if (!length(v)) return(numeric(0))
  if (length(v) %% 31L != 0L) return(numeric(0))  # malformed; treated as no data
  # Column 2 of each row is l = 1, which is what 01f_ takes as x[2].
  v[seq(2L, length(v), by = 31L)]
}
feats <- function(acc) { n <- length(acc)
  c(prox = mean(acc[1:min(30, n)]), dist = mean(acc[max(1, n - 29):n]),
    mean = mean(acc), min = min(acc), sd = stats::sd(acc)) }

STRUCT <- list(
  utr5 = list(dir = "rnaplfold_output", cols = c(
    prox = "struct_accessibility_cap_proximal", dist = "struct_accessibility_aug_context",
    mean = "struct_accessibility_utr5_mean",    min  = "struct_accessibility_utr5_min",
    sd   = "struct_accessibility_utr5_sd"),     nsr = "struct_num_structured_regions"),
  cds  = list(dir = "cds_rnaplfold_output", cols = c(
    prox = "struct_accessibility_start_proximal_cds", dist = "struct_accessibility_stop_proximal_cds",
    mean = "struct_accessibility_cds_mean",     min  = "struct_accessibility_cds_min",
    sd   = "struct_accessibility_cds_sd"),      nsr = "struct_num_structured_regions_cds"),
  utr3 = list(dir = "utr3_rnaplfold_output", cols = c(
    prox = "struct_accessibility_stop_proximal_utr3", dist = "struct_accessibility_distal_utr3",
    mean = "struct_accessibility_utr3_mean",    min  = "struct_accessibility_utr3_min",
    sd   = "struct_accessibility_utr3_sd"),     nsr = "struct_num_structured_regions_utr3"))

# Cached PER REGION. Reading ~34,000 files takes tens of minutes and the process has been killed
# partway through more than once; without per-region caching each kill threw away every region
# that had already finished.
struct_for <- function(want_tx, cache_prefix = NULL) {
  out <- tibble(transcript_id_clean = want_tx)
  for (rg in names(STRUCT)) {
    sp <- STRUCT[[rg]]
    rgc <- if (!is.null(cache_prefix)) paste0(cache_prefix, "_", rg, ".rds") else NULL
    if (!is.null(rgc) && file.exists(rgc)) {
      add <- readRDS(rgc)
      stopifnot("a cached region does not match the requested transcripts" =
                  identical(nrow(add), length(want_tx)))
      out <- bind_cols(out, add); cat(sprintf("    %s loaded from cache\n", rg)); next
    }
    fp <- file.path(plf, sp$dir, paste0(want_tx, "_lunp"))
    ok <- file.exists(fp)
    mat <- matrix(NA_real_, nrow = length(want_tx), ncol = 6,
                  dimnames = list(NULL, c(names(sp$cols), "nsr")))
    wi <- which(ok)
    for (j in seq_along(wi)) {
      i <- wi[j]
      if (j %% 2000 == 0) { cat(sprintf("\r    %s %d/%d", rg, j, length(wi))); flush.console() }
      a <- read_acc(fp[i])
      a <- a[!is.na(a)]
      if (!length(a)) next
      mat[i, 1:5] <- feats(a)[names(sp$cols)]
      # num_structured_regions, as 01f_ defines it after the phantom row is gone.
      mat[i, 6] <- length(a)
    }
    cat(sprintf("\r    %s %d/%d done\n", rg, length(wi), length(wi)))
    add <- as_tibble(mat)
    names(add) <- c(unname(sp$cols), sp$nsr)
    if (!is.null(rgc)) saveRDS(add, rgc)
    out <- bind_cols(out, add)
    invisible(gc(verbose = FALSE))
  }
  out
}

if (do_stage("struct")) {
# --- reproduction gate: 231 -------------------------------------------------------------------
rec_f <- here("output", "translational_control_overview", "structure_recomputed.rds")
stopifnot("01f_'s structure_recomputed.rds not found - run 01f_recompute_structure.R first" =
            file.exists(rec_f))
rec <- readRDS(rec_f)
set.seed(9)
probe <- sample(rec$transcript_id_clean, min(300, nrow(rec)))
mine  <- struct_for(probe)
cmp_cols <- intersect(sub("^tco_", "", colnames(rec)), colnames(mine))
cmp_cols <- setdiff(cmp_cols, "transcript_id_clean")
ref <- rec %>% filter(transcript_id_clean %in% probe) %>%
  rename_with(~ sub("^tco_", "", .x)) %>% arrange(transcript_id_clean)
got <- mine %>% arrange(transcript_id_clean)
devs <- map_dbl(cmp_cols, function(c_) {
  x <- ref[[c_]]; y <- got[[c_]]; k <- is.finite(x) & is.finite(y)
  if (!any(k)) return(NA_real_); max(abs(x[k] - y[k])) })
names(devs) <- cmp_cols
cat("\nstructure reproduction gate (vs 01f_ stored 231 values), max|dev| per column:\n")
print(round(devs, 12))
stopifnot("the structure computation here does not reproduce 01f_'s 231 values" =
            all(devs < 1e-10, na.rm = TRUE),
          "the structure gate compared too few columns" = length(cmp_cols) >= 12)
cat("structure reproduction gate PASSED on", length(probe), "transcripts\n")

# Checkpointed. This step reads ~34,000 _lunp files and two earlier attempts were killed here
# with no R error at all, which is what an out-of-memory kill looks like through a pipe. Caching
# makes the script resumable and, more usefully, isolates which stage is responsible.
struct_cache <- file.path(PD, "mcf7six1_struct_from_lunp.rds")
if (file.exists(struct_cache)) {
  struct6 <- readRDS(struct_cache)
  cat("loaded cached MCF7-SIX1 struct columns\n")
} else {
  cat("computing MCF7-SIX1 struct columns from _lunp files...\n"); flush.console()
  struct6 <- struct_for(m6$transcript_id_clean,
                        cache_prefix = file.path(PD, "mcf7six1_struct_region"))
  saveRDS(struct6, struct_cache)
  cat("cached MCF7-SIX1 struct columns\n")
}
invisible(gc(verbose = FALSE))
cat(sprintf("MCF7-SIX1 struct coverage: utr5 %.1f%% | cds %.1f%% | utr3 %.1f%%\n",
  100*mean(!is.na(struct6$struct_accessibility_utr5_mean)),
  100*mean(!is.na(struct6$struct_accessibility_cds_mean)),
  100*mean(!is.na(struct6$struct_accessibility_utr3_mean))))
}  # end struct stage

# ---------------------------------------------------------------------------------------------
# 2. positional CSC (iCodon), as 01j_ computes it
# ---------------------------------------------------------------------------------------------
if (do_stage("csc")) {
need_cds()
predict_human <- iCodon::predict_stability("human")
# Chunked. An unchunked call over ~11,200 sequences x 5 CDS segments died silently on an 8 GB
# machine - no error, just a dead process - so the batch size is a memory ceiling, not a
# preference. Progress is printed so a stall is distinguishable from a crash next time.
CSC_CHUNK <- 1000L
csc_of <- function(seqs_chr, label = "") {
  n <- length(seqs_chr)
  if (!n) return(numeric(0))
  # iCodon's preprocess_secuences() ERRORS on duplicate input ("Input sequences should be
  # unique, repeats found"), and duplicates are common here: many transcripts share an identical
  # CDS quarter. Scoring the unique set and mapping back is both required and cheaper. The
  # 200-transcript reproduction gate happened to contain no duplicates, which is why this only
  # surfaced on the full set - a reminder that a gate on a sample does not exercise every path.
  u <- unique(seqs_chr)
  out_u <- numeric(length(u))
  idx <- split(seq_along(u), ceiling(seq_along(u) / CSC_CHUNK))
  for (i in seq_along(idx)) {
    out_u[idx[[i]]] <- suppressWarnings(as.numeric(predict_human(u[idx[[i]]])))
    cat(sprintf("\r    %s %d/%d unique", label, min(i * CSC_CHUNK, length(u)), length(u)))
    flush.console(); gc(verbose = FALSE)
  }
  cat(sprintf("\r    %s %d unique of %d done\n", label, length(u), n))
  out_u[match(seqs_chr, u)]
}
positional_csc <- function(want_tx) {
  s <- cds_all[names(cds_all) %in% want_tx]
  ch <- as.character(s); nm <- names(s)
  w  <- nchar(ch); ncod <- w %/% 3
  sub_of <- function(from_cod, to_cod) {
    a <- pmax(1, from_cod) * 3 - 2; b <- pmin(ncod, to_cod) * 3
    ifelse(b >= a, substring(ch, a, b), NA_character_)
  }
  q <- function(i) { lo <- floor((i - 1) * ncod / 4) + 1; hi <- floor(i * ncod / 4); sub_of(lo, hi) }
  parts <- list(csc_first75 = ifelse(w >= 225, substring(ch, 1, 225), NA_character_),
                csc_q1 = q(1), csc_q2 = q(2), csc_q3 = q(3), csc_q4 = q(4))
  res <- tibble(transcript_id_clean = nm)
  for (nmx in names(parts)) {
    v <- rep(NA_real_, length(nm)); k <- !is.na(parts[[nmx]]) & nchar(parts[[nmx]]) >= 30
    if (any(k)) v[k] <- csc_of(parts[[nmx]][k], nmx)
    res[[nmx]] <- v
  }
  res
}

# --- reproduction gate: 231 -------------------------------------------------------------------
# Only defined by the struct stage; suppress the no-such-object warnings in other stages.
suppressWarnings(rm(list = intersect(c("rec", "mine", "ref", "got"), ls())))
invisible(gc(verbose = FALSE))
set.seed(9)
probe2 <- sample(b$transcript_id_clean[!is.na(b$csc_q1)], min(200, sum(!is.na(b$csc_q1))))
pc <- positional_csc(probe2)
ref2 <- b %>% filter(transcript_id_clean %in% probe2) %>%
  dplyr::select(transcript_id_clean, csc_first75, csc_q1, csc_q2, csc_q3, csc_q4) %>%
  arrange(transcript_id_clean)
got2 <- pc %>% arrange(transcript_id_clean)
devs2 <- map_dbl(c("csc_first75","csc_q1","csc_q2","csc_q3","csc_q4"), function(c_) {
  x <- ref2[[c_]]; y <- got2[[c_]]; k <- is.finite(x) & is.finite(y)
  if (!any(k)) return(NA_real_); max(abs(x[k] - y[k])) })
names(devs2) <- c("csc_first75","csc_q1","csc_q2","csc_q3","csc_q4")
cat("\npositional CSC reproduction gate (vs the 231 baseline), max|dev| per column:\n")
print(signif(devs2, 6))
stopifnot("the positional CSC computation here does not reproduce the 231 baseline values" =
            all(devs2 < 1e-6, na.rm = TRUE))
cat("positional CSC reproduction gate PASSED on", length(probe2), "transcripts\n")

csc_cache <- file.path(PD, "mcf7six1_positional_csc.rds")
if (file.exists(csc_cache)) {
  csc6 <- readRDS(csc_cache); cat("loaded cached MCF7-SIX1 positional CSC\n")
} else {
  cat("computing MCF7-SIX1 positional CSC...\n"); flush.console()
  csc6 <- positional_csc(m6$transcript_id_clean)
  saveRDS(csc6, csc_cache); cat("cached MCF7-SIX1 positional CSC\n")
}
invisible(gc(verbose = FALSE))
}  # end csc stage

# ---------------------------------------------------------------------------------------------
# 3. gene- and symbol-level carry-over from the 231 baseline
# ---------------------------------------------------------------------------------------------
if (do_stage("assemble")) {
# Read the checkpoints rather than recomputing, so assembly is cheap and repeatable.
struct_cache <- file.path(PD, "mcf7six1_struct_from_lunp.rds")
csc_cache    <- file.path(PD, "mcf7six1_positional_csc.rds")
stopifnot("run the struct stage first" = file.exists(struct_cache),
          "run the csc stage first"    = file.exists(csc_cache))
struct6 <- readRDS(struct_cache); csc6 <- readRDS(csc_cache)

missing_cols <- setdiff(colnames(b), colnames(m6))
computed     <- c(colnames(struct6), colnames(csc6), "transcript_id")
carry        <- setdiff(missing_cols, computed)

# FORCED CARRY-OVER. These already exist in the MCF7-SIX1 matrix, but were computed there by
# 43_feature_extraction_mcf7six1_codon.Rmd under a different convention from the 231 pipeline's,
# so keeping them would put two different quantities under one column name.
#
# cnot3_weighted_codon_score is the clear case: notebook 22 defines it per GENE and asserts it
# "must be non-negative", and the MCF7-SIX1 value for DHX29 is -0.0295 - it fails the 231
# definition's own invariant. Measured across the 5,809 transcripts dominant in both cell lines,
# the two versions differ by up to 0.208 on EVERY one. Notebook 22 keys it on ensembl_gene, so a
# gene-level carry-over is the definition rather than an approximation.
#
# The cost is coverage: genes absent from the 231 baseline get NA and are imputed like any other
# missing feature. That is the right trade - an identical-but-sparser column is comparable, a
# fully-populated but differently-defined one is not.
FORCE_CARRY <- c("cnot3_weighted_codon_score")
force_ok <- intersect(FORCE_CARRY, colnames(b))
carry <- union(carry, force_ok)
if (length(force_ok)) {
  cat("forced carry-over from the 231 baseline (definition differs in the MCF7 matrix):",
      paste(force_ok, collapse = ", "), "\n")
  m6 <- m6 %>% dplyr::select(-any_of(force_ok))
}

by_gene <- b %>% dplyr::select(gene_id_clean, all_of(carry)) %>%
  group_by(gene_id_clean) %>% slice(1) %>% ungroup()
stopifnot("the 231 baseline is not one row per gene - a gene-level carry-over would be ambiguous" =
            nrow(by_gene) == n_distinct(b$gene_id_clean))

out <- m6 %>%
  dplyr::select(any_of(colnames(b))) %>%
  left_join(struct6, by = "transcript_id_clean") %>%
  left_join(csc6,    by = "transcript_id_clean") %>%
  left_join(by_gene, by = "gene_id_clean") %>%
  mutate(transcript_id = transcript_id_clean)

cat(sprintf("\ncarried over by gene: %d columns | coverage of MCF7-SIX1 genes in the 231 baseline: %.1f%%\n",
            length(carry), 100 * mean(m6$gene_id_clean %in% by_gene$gene_id_clean)))

# ---------------------------------------------------------------------------------------------
# 4. the parity assertion this whole script exists to satisfy
# ---------------------------------------------------------------------------------------------
only_here <- setdiff(colnames(out), colnames(b))
only_231  <- setdiff(colnames(b), colnames(out))
if (length(only_here)) cat("only in MCF7-SIX1:", paste(only_here, collapse = ", "), "\n")
if (length(only_231))  cat("only in 231:",      paste(only_231,  collapse = ", "), "\n")
out <- out %>% dplyr::select(all_of(colnames(b)))       # same columns, same ORDER
stopifnot(
  "the MCF7-SIX1 baseline does not carry exactly the 231 baseline's columns" =
    identical(colnames(out), colnames(b)),
  "transcript_id_clean is not unique" = !any(duplicated(out$transcript_id_clean)),
  "the gene identifier trio is incomplete" =
    all(c("gene_id", "gene_id_clean", "symbol") %in% colnames(out)),
  "accessibility outside [0,1]" =
    all(out$struct_accessibility_cds_mean >= 0 & out$struct_accessibility_cds_mean <= 1,
        na.rm = TRUE),
  "GC content outside [0,100]" = all(out$cds_gc >= 0 & out$cds_gc <= 100, na.rm = TRUE)
)
saveRDS(out, file.path(PD, "feature_matrix_mcf7six1_dhx29_kd_feature.rds"))
cat(sprintf("\nwrote feature_matrix_mcf7six1_dhx29_kd_feature.rds - %d x %d (231 baseline: %d x %d)\n",
            nrow(out), ncol(out), nrow(b), ncol(b)))
}  # end assemble stage

# ---------------------------------------------------------------------------------------------
# 5. Kozak sequence cache, in the shape 01_ expects
# ---------------------------------------------------------------------------------------------
if (do_stage("kozak")) {
need_cds(); need_utr5()
ref_seqs <- readRDS(file.path(PD, "precomputed_kozak_sequences.rds"))
stopifnot("the 231 Kozak cache is not the expected two-table list" =
            all(c("utr5_seq", "cds_seq") %in% names(ref_seqs)))
u5col <- setdiff(colnames(ref_seqs$utr5_seq), "transcript_id_clean")
cdcol <- setdiff(colnames(ref_seqs$cds_seq),  "transcript_id_clean")
# Keyed on m6, not on `out`: the two carry identical rows by construction, and depending on
# `out` would tie this stage to the assemble stage for no reason.
want <- m6$transcript_id_clean
mk <- function(sset, col) {
  s <- sset[names(sset) %in% want]
  tibble(transcript_id_clean = names(s), !!col := as.character(s))
}
seqs6 <- list(utr5_seq = mk(utr5_all, u5col), cds_seq = mk(cds_all, cdcol))
stopifnot("Kozak cache does not cover every baseline transcript" =
            nrow(seqs6$utr5_seq) == length(want) && nrow(seqs6$cds_seq) == length(want))
bl <- file.path(PD, "feature_matrix_mcf7six1_dhx29_kd_feature.rds")
if (file.exists(bl))
  stopifnot("the Kozak cache does not match the written baseline's transcripts" =
              setequal(want, readRDS(bl)$transcript_id_clean))
saveRDS(seqs6, file.path(PD, "precomputed_mcf7six1_kozak_sequences.rds"))
cat(sprintf("wrote precomputed_mcf7six1_kozak_sequences.rds - utr5 %d, cds %d rows\n",
            nrow(seqs6$utr5_seq), nrow(seqs6$cds_seq)))
}  # end kozak stage
cat("\nNext: G4mer summaries for these transcripts, then 01_build_feature_matrix.Rmd with\n")
cat("  cell_line='MCF7-SIX1', baseline_rds='feature_matrix_mcf7six1_dhx29_kd_feature.rds',\n")
cat("  kozak_seqs_rds='precomputed_mcf7six1_kozak_sequences.rds',\n")
cat("  dhx29_occ_rds='dhx29_riboseq_occupancy_mcf7six1.rds', g4mer_tag='_mcf7six1',\n")
cat("  cache_tag='_mcf7six1', out_rds='feature_matrix_tco_mcf7six1.rds'\n")
