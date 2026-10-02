# 01l_build_celltx_baseline.R
# Build the counterpart of feature_matrix_dhx29_kd_feature.rds for a dataset's OWN transcript set
# (one transcript per gene), so 01_build_feature_matrix.Rmd can produce that dataset's matrix from
# the same code that produces the MDA-MB-231 one.
#
#   Rscript 01l_build_celltx_baseline.R <cell> <stage>
#     cell  = hela  (precomputed_teleman_tx.rds: Roiuk-dominant HeLa transcript)
#           | dhx29 (precomputed_dhx29_tx.rds: HEK293T representative, dhx29_reanalysis 03)
#     stage = struct | csc | seq | assemble | kozak   (one process each; see 01i_ on why)
#
# WHY A NEW SCRIPT RATHER THAN 01i_. 01i_ starts from an existing MCF7-SIX1 matrix built by
# notebook 43 with every sequence column already present. HeLa and HEK293T have no such matrix
# (the Teleman matrix disagrees with the 231 baseline on kozak_optimal for half of the shared
# transcripts), so every TRANSCRIPT-PROPERTY column is computed here from sequence.
#
# EVERY COMPUTED COLUMN IS GATED. Before the target set is computed, the same code is run on a
# sample of the 231 baseline's own transcripts and must reproduce the stored values (structure:
# against 01f_'s recomputation; everything else: against the baseline). A second implementation
# of a feature definition is how two cell lines end up measuring different quantities under one
# column name; the gate is what makes this one the same definition.
#
# Definitions reproduced (source in brackets):
#   lengths, log2 lengths, regional GC, transcript_gc  [predictive_modeling/01_feature_extraction]
#   gc3 (whole CDS, stop included)                     [01_feature_extraction calculate_gc3]
#   tai                                                [01g_ tAI chunk]
#   csc (whole CDS, length %% 3 == 0)                  [01g_ CSC chunk]
#   csc_first75, csc_q1..q4                            [01j_, via 01i_ positional_csc]
#   kozak_score, kozak_optimal, kozak_pwm_score        [01h_]
#   struct_* (18)                                      [01f_recompute_structure / 01i_]
# cai / fop are left NA: 01_ REPLACES them for every non-231 build with the shared-weight values.
# Everything else (231 TE labels, gene-level external data) is carried over by GENE from the 231
# baseline, exactly as 01i_ does for MCF7-SIX1; none of it enters an intrinsic model.

suppressPackageStartupMessages({
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38); library(Biostrings)
  library(dplyr); library(tidyr); library(purrr); library(readr); library(stringr); library(here)
})
.a <- commandArgs(trailingOnly = TRUE)
CELL  <- if (length(.a) >= 1) .a[1] else stop("usage: 01l_build_celltx_baseline.R <hela|dhx29|mane|mcf7six1> <stage>")
STAGE <- if (length(.a) >= 2) .a[2] else stop("give a stage: struct | csc | seq | assemble | kozak")
# mane = MANE Select for the Weber DAP5 genes (10w_weber_dap5_outcome_prep.Rmd), the one set not
# chosen by abundance in its own data.
TX_RDS <- c(hela = "precomputed_teleman_tx.rds", dhx29 = "precomputed_dhx29_tx.rds",
            mane = "precomputed_mane_tx.rds",
            # MCF7-SIX1 moved onto this builder 2026-10-01: its older one (01i_) shared the
            # short-5'UTR reading bug (650 transcripts lost) and took CDS GC/Kozak from 43_,
            # which blanks them when the CDS length is not divisible by 3 (105 transcripts).
            mcf7six1 = "precomputed_mcf7six1_tx.rds",
            # mdamb231: struct stage only, read by 01_ (late_struct_rds) to fill MDA-MB-231
            # transcripts whose _lunp files were written after 01f_ ran.
            mdamb231 = "precomputed_mdamb231_tx.rds")
stopifnot("unknown cell set" = CELL %in% names(TX_RDS),
          "unknown stage" = STAGE %in% c("struct", "csc", "seq", "assemble", "kozak"))
cat("cell:", CELL, "| stage:", STAGE, "\n"); flush.console()

PD  <- here("output", "predictive_modeling")
plf <- here("accessories", "plfold_output", "rnaplfold_output")
b   <- readRDS(file.path(PD, "feature_matrix_dhx29_kd_feature.rds"))
txm <- readRDS(file.path(PD, TX_RDS[[CELL]])) %>% dplyr::select(transcript_id_clean, gene_id_clean, symbol)
stopifnot("transcript map is not one transcript per gene" = !any(duplicated(txm$gene_id_clean)),
          "transcript map has duplicated transcripts" = !any(duplicated(txm$transcript_id_clean)))
cache <- function(x) file.path(PD, paste0(CELL, "_", x, ".rds"))

txdb   <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
genome <- BSgenome.Hsapiens.UCSC.hg38
ids_of <- function(grl) unique(sub("\\..*", "", names(grl)))
# 01_ asserts a 5'UTR and a CDS for every transcript in a matrix, so the universe is restricted
# here, as in 01i_.
tx <- txm$transcript_id_clean
tx <- tx[tx %in% ids_of(cdsBy(txdb, by = "tx", use.names = TRUE)) &
         tx %in% ids_of(fiveUTRsByTranscript(txdb, use.names = TRUE))]
cat(sprintf("%s transcripts: %d -> %d with both a 5'UTR and a CDS\n", CELL, nrow(txm), length(tx)))
txm <- txm %>% filter(transcript_id_clean %in% tx)

# The 231 probe used by every gate: drawn once, so each stage tests the same transcripts.
set.seed(9)
PROBE <- sample(b$transcript_id_clean, 300)
SEQ_TX <- union(tx, PROBE)
grab <- function(grl) {
  names(grl) <- sub("\\..*", "", names(grl))
  grl <- grl[names(grl) %in% SEQ_TX]
  s <- extractTranscriptSeqs(genome, grl); names(s) <- sub("\\..*", "", names(s)); s
}
gate <- function(got, cols, tol, label) {
  ref <- b %>% filter(transcript_id_clean %in% got$transcript_id_clean) %>%
    dplyr::select(transcript_id_clean, all_of(cols))
  j <- inner_join(got, ref, by = "transcript_id_clean", suffix = c(".new", ".ref"))
  devs <- map_dbl(cols, function(cc) {
    x <- as.numeric(unlist(j[[paste0(cc, ".new")]])); y <- as.numeric(unlist(j[[paste0(cc, ".ref")]]))
    na_mismatch <- sum(is.na(x) != is.na(y)); k <- is.finite(x) & is.finite(y)
    if (na_mismatch > 0) return(Inf)
    if (!any(k)) return(NA_real_); max(abs(x[k] - y[k])) })
  names(devs) <- cols
  cat(sprintf("\n%s reproduction gate vs the 231 baseline (%d transcripts), max|dev|:\n", label, nrow(j)))
  print(signif(devs, 4))
  stopifnot("a gate compared too few transcripts" = nrow(j) >= 200,
            "a computed column does not reproduce the 231 baseline (Inf = NA pattern differs)" =
              all(devs <= tol, na.rm = TRUE))
  cat(label, "gate PASSED\n")
}

# =============================================================================================
# struct: 18 struct_* columns from the _lunp files (01i_ code, unchanged)
# =============================================================================================
if (STAGE == "struct") {
read_acc <- function(fp) {
  hdr <- readLines(fp, n = 2L, warn = FALSE)
  h   <- sum(grepl("^[[:space:]]*#", hdr))
  # Width is read from the file, not assumed: RNAplfold caps -u at the region length, so a
  # 5'UTR shorter than 30 nt has fewer than 31 columns. Assuming 31 silently dropped them (590
  # HEK293T 5'UTRs). 01f_ takes column 2 of each line whatever the width; so does this.
  first <- readLines(fp, n = h + 1L, warn = FALSE)[h + 1L]
  if (is.na(first) || !nzchar(trimws(first))) return(numeric(0))
  w <- length(strsplit(first, "\t", fixed = FALSE)[[1]])
  v <- scan(fp, what = double(), skip = h, na.strings = "NA", quiet = TRUE)
  if (!length(v) || length(v) %% w != 0L) return(numeric(0))
  v[seq(2L, length(v), by = w)]
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
struct_for <- function(want_tx, cache_prefix = NULL) {
  out <- tibble(transcript_id_clean = want_tx)
  for (rg in names(STRUCT)) {
    sp <- STRUCT[[rg]]
    rgc <- if (!is.null(cache_prefix)) paste0(cache_prefix, "_", rg, ".rds") else NULL
    if (!is.null(rgc) && file.exists(rgc)) {
      add <- readRDS(rgc); stopifnot("cached region does not match" = nrow(add) == length(want_tx))
      out <- bind_cols(out, add); cat(sprintf("    %s loaded from cache\n", rg)); next
    }
    fp <- file.path(plf, sp$dir, paste0(want_tx, "_lunp")); ok <- file.exists(fp)
    mat <- matrix(NA_real_, nrow = length(want_tx), ncol = 6, dimnames = list(NULL, c(names(sp$cols), "nsr")))
    wi <- which(ok)
    for (j in seq_along(wi)) {
      i <- wi[j]
      if (j %% 2000 == 0) { cat(sprintf("\r    %s %d/%d", rg, j, length(wi))); flush.console() }
      a <- read_acc(fp[i]); a <- a[!is.na(a)]; if (!length(a)) next
      mat[i, 1:5] <- feats(a)[names(sp$cols)]; mat[i, 6] <- length(a)
    }
    cat(sprintf("\r    %s %d/%d done\n", rg, length(wi), length(wi)))
    add <- as_tibble(mat); names(add) <- c(unname(sp$cols), sp$nsr)
    if (!is.null(rgc)) saveRDS(add, rgc)
    out <- bind_cols(out, add); invisible(gc(verbose = FALSE))
  }
  out
}
rec <- readRDS(here("output", "translational_control_overview", "structure_recomputed.rds"))
probe <- intersect(PROBE, rec$transcript_id_clean)
mine <- struct_for(probe)
cmp <- setdiff(intersect(sub("^tco_", "", colnames(rec)), colnames(mine)), "transcript_id_clean")
ref <- rec %>% filter(transcript_id_clean %in% probe) %>% rename_with(~ sub("^tco_", "", .x)) %>%
  arrange(transcript_id_clean)
got <- mine %>% arrange(transcript_id_clean)
devs <- map_dbl(cmp, function(c_) { x <- ref[[c_]]; y <- got[[c_]]; k <- is.finite(x) & is.finite(y)
  if (!any(k)) NA_real_ else max(abs(x[k] - y[k])) }); names(devs) <- cmp
cat("structure gate (vs 01f_), max|dev|:\n"); print(round(devs, 12))
# One-directional on purpose. NA here but a value in 01f_ = this reader dropped a file (the
# short-5'UTR bug). A value here but NA in 01f_ = a _lunp file written after 01f_ ran (the
# MCF7-SIX1/HeLa/HEK293T precomputes add files to the same directories), i.e. extra coverage.
dropped <- map_int(cmp, function(c_) sum(is.na(got[[c_]]) & !is.na(ref[[c_]]))); names(dropped) <- cmp
added   <- map_int(cmp, function(c_) sum(!is.na(got[[c_]]) & is.na(ref[[c_]])))
cat(sprintf("NA pattern vs 01f_: %d values dropped by this reader | up to %d transcripts with files newer than 01f_\n",
            sum(dropped), max(added)))
stopifnot("structure does not reproduce 01f_" = all(devs < 1e-10, na.rm = TRUE), length(cmp) >= 12,
          "this reader left NA where 01f_ has a value (a file was read as empty)" = all(dropped == 0))
cat("structure gate PASSED on", length(probe), "transcripts\n")
st <- struct_for(tx, cache_prefix = file.path(PD, paste0(CELL, "_struct_region")))
saveRDS(st, cache("struct_from_lunp"))
cat(sprintf("%s struct coverage: utr5 %.1f%% | cds %.1f%% | utr3 %.1f%%\n", CELL,
  100*mean(!is.na(st$struct_accessibility_utr5_mean)), 100*mean(!is.na(st$struct_accessibility_cds_mean)),
  100*mean(!is.na(st$struct_accessibility_utr3_mean))))
}

# =============================================================================================
# csc: whole-CDS csc (01g_) and positional csc (01j_ via 01i_)
# =============================================================================================
if (STAGE == "csc") {
suppressPackageStartupMessages(library(iCodon))
cds_all <- grab(cdsBy(txdb, by = "tx", use.names = TRUE))
predict_human <- iCodon::predict_stability("human")
csc_of <- function(seqs_chr, label = "") {
  if (!length(seqs_chr)) return(numeric(0))
  u <- unique(seqs_chr); out_u <- numeric(length(u))
  idx <- split(seq_along(u), ceiling(seq_along(u) / 1000L))
  for (i in seq_along(idx)) {
    out_u[idx[[i]]] <- suppressWarnings(as.numeric(predict_human(u[idx[[i]]])))
    cat(sprintf("\r    %s %d/%d unique", label, min(i * 1000L, length(u)), length(u))); flush.console(); gc(verbose = FALSE)
  }
  cat("\n"); out_u[match(seqs_chr, u)]
}
# Exactly 01g_ (whole CDS) and 01j_ (positional), including their NA rules: only CDS whose length
# is a multiple of 3 are scored at all; csc_first75 needs >= 225 nt; quarter boundaries are 01j_'s
# codon-aligned ones and a segment needs >= 1 codon. (01i_'s positional_csc used a >= 30 nt rule
# and no frame filter; its gate compared only where both values existed, so the NA pattern was
# never tested. The positional columns are not modelled, so no model was affected.)
csc_set <- function(want_tx) {
  s <- cds_all[names(cds_all) %in% want_tx]
  res <- tibble(transcript_id_clean = names(s), csc = NA_real_, csc_first75 = NA_real_,
                csc_q1 = NA_real_, csc_q2 = NA_real_, csc_q3 = NA_real_, csc_q4 = NA_real_)
  s3 <- s[width(s) %% 3 == 0]; ch <- as.character(s3); nm <- names(s3); w <- nchar(ch)
  put <- function(col, vals, ids) res[[col]][match(ids, res$transcript_id_clean)] <<- vals
  put("csc", csc_of(ch, "csc"), nm)
  long <- w >= 225
  if (any(long)) put("csc_first75", csc_of(substring(ch[long], 1, 225), "csc_first75"), nm[long])
  n_cod <- as.integer(w / 3)
  b1 <- pmax(floor(n_cod * 0.25) * 3L, 3L); b2 <- pmax(floor(n_cod * 0.50) * 3L, b1 + 3L)
  b3 <- pmax(floor(n_cod * 0.75) * 3L, b2 + 3L); b4 <- as.integer(w)
  segs <- list(csc_q1 = list(rep(1L, length(ch)), b1), csc_q2 = list(b1 + 1L, b2),
               csc_q3 = list(b2 + 1L, b3), csc_q4 = list(b3 + 1L, b4))
  for (q in names(segs)) {
    st <- segs[[q]][[1]]; en <- segs[[q]][[2]]; ok <- (en - st + 1L) >= 3L
    if (any(ok)) put(q, csc_of(substring(ch[ok], st[ok], en[ok]), q), nm[ok])
  }
  res
}
gate(csc_set(PROBE), c("csc", "csc_first75", "csc_q1", "csc_q2", "csc_q3", "csc_q4"), 1e-6, "CSC")
saveRDS(csc_set(tx), cache("csc"))
cat("cached", CELL, "csc\n")
}

# =============================================================================================
# seq: lengths, GC, transcript_gc, gc3, tai, Kozak
# =============================================================================================
if (STAGE == "seq") {
suppressPackageStartupMessages(library(coRdon))
grl5 <- fiveUTRsByTranscript(txdb, use.names = TRUE); grlc <- cdsBy(txdb, by = "tx", use.names = TRUE)
grl3 <- threeUTRsByTranscript(txdb, use.names = TRUE)
utr5 <- grab(grl5); cds <- grab(grlc); utr3 <- grab(grl3)
len_of <- function(grl) { names(grl) <- sub("\\..*", "", names(grl)); grl <- grl[names(grl) %in% SEQ_TX]
  tibble(transcript_id_clean = names(grl), len = as.integer(sum(width(grl)))) }
gc_content <- function(seq) { seq <- str_to_upper(seq); g <- str_count(seq, "[GC]"); n <- str_length(seq)
  r <- (g / n) * 100; r[n == 0] <- NA; r }
calculate_gc3 <- function(seqs) { seqs <- str_to_upper(seqs)
  vapply(seqs, function(seq) { if (is.na(seq) || nchar(seq) < 3L) return(NA_real_)
    n_codons <- nchar(seq) %/% 3L; starts <- seq.int(3L, n_codons * 3L, 3L)
    thirds <- substr(rep.int(seq, n_codons), starts, starts); sum(thirds %in% c("G", "C")) * 100 / n_codons
  }, numeric(1L), USE.NAMES = FALSE) }
sq <- function(s) tibble(transcript_id_clean = names(s), seq = as.character(s))

# tAI exactly as 01g_
gen_code <- Biostrings::getGeneticCode("SGC0")
sense_codons <- names(gen_code[gen_code != "*"])
trna_gcn <- read_csv(here("accessories", "human", "human_trna_gcn.csv"), comment = "#", show_col_types = FALSE)
decode_anticodon <- function(ac) {
  x1 <- substr(ac, 1, 1); x2 <- substr(ac, 2, 2); x3 <- substr(ac, 3, 3)
  wc <- c(A = "T", T = "A", G = "C", C = "G"); y1 <- unname(wc[x3]); y2 <- unname(wc[x2])
  if (x1 == "A") { y3 <- c("T", "C", "A"); svs <- c(0, 0, 0) } else if (x1 == "G") { y3 <- c("C", "T"); svs <- c(0, 0.41) } else if (x1 == "C") { y3 <- "G"; svs <- 0 } else { y3 <- "A"; svs <- 0 }
  data.frame(codon = paste0(y1, y2, y3), s = svs, stringsAsFactors = FALSE) }
aa_of_codon <- setNames(as.character(gen_code), names(gen_code))
ws_tai <- setNames(rep(0, length(gen_code)), names(gen_code))
for (i in seq_len(nrow(trna_gcn))) { d <- decode_anticodon(trna_gcn$anticodon[i])
  for (j in seq_len(nrow(d))) { cd <- d$codon[j]
    if (!is.na(aa_of_codon[cd]) && aa_of_codon[cd] == trna_gcn$amino_acid[i]) ws_tai[cd] <- ws_tai[cd] + (1 - d$s[j]) * trna_gcn$gcn[i] } }
w_tai <- ws_tai / max(ws_tai[sense_codons])
tai_of <- function(cds_set) {
  s <- cds_set[width(cds_set) %% 3 == 0]
  cm <- codonCounts(codonTable(s))
  use <- sense_codons[!is.na(w_tai[sense_codons]) & w_tai[sense_codons] > 0]
  v <- apply(cm, 1, function(gc_) { c_ <- gc_[use]; n <- sum(c_); if (n == 0) NA_real_ else exp(sum(c_ * log(w_tai[use])) / n) })
  tibble(transcript_id_clean = names(s), tai = unname(v))
}

# Kozak exactly as 01h_
match_kozak <- function(seq) {
  if (nchar(seq) < 10) return("unknown_short")
  if (grepl("^(GCCACC|GCCGCC)ATGG$", seq)) return("optimal")
  if (grepl("^...[AG]..ATGG$", seq))        return("strong")
  if (grepl("^...[AG]..ATG[ACT]$", seq))    return("moderate")
  if (grepl("^...[CT]..ATGG$", seq))        return("moderate")
  return("unknown") }
score_map <- c(unknown = 0L, moderate = 1L, strong = 2L, optimal = 3L, unknown_short = NA_integer_)
kozak_pwm <- matrix(c(0.23, 0.22, 0.27, 0.28, 0.25, 0.28, 0.23, 0.24, 0.21, 0.34, 0.23, 0.22,
  0.71, 0.09, 0.13, 0.07, 0.27, 0.28, 0.19, 0.26, 0.25, 0.39, 0.15, 0.21, 0.17, 0.12, 0.54, 0.17),
  nrow = 7, ncol = 4, byrow = TRUE, dimnames = list(c("m6","m5","m4","m3","m2","m1","p4"), c("A","C","G","T")))
# This is 01h_'s matrix as written there, INCLUDING its known errors at -3 and -6 (CLAUDE.md):
# kozak_pwm_score is the superseded 'current' variant and must reproduce the 231 column, which
# the gate below checks. The corrected kozak_pwm_score_v2 is computed by 01_ from the file.
log_odds <- log2((kozak_pwm + 0.001) / 0.25)
pos_idx <- c(1, 2, 3, 4, 5, 6, 10); pos_nm <- c("m6","m5","m4","m3","m2","m1","p4")
score_pwm <- function(seq) { if (nchar(seq) < 10) return(NA_real_); nts <- strsplit(seq, "")[[1]]; s <- 0
  for (i in seq_along(pos_idx)) { nt <- nts[pos_idx[i]]; if (!nt %in% c("A","C","G","T")) return(NA_real_); s <- s + log_odds[pos_nm[i], nt] }; s }

seq_set <- function(want_tx) {
  L5 <- len_of(grl5) %>% dplyr::rename(utr5_length = len); LC <- len_of(grlc) %>% dplyr::rename(cds_length = len)
  L3 <- len_of(grl3) %>% dplyr::rename(utr3_length = len)
  G5 <- sq(utr5) %>% transmute(transcript_id_clean, utr5_gc = gc_content(seq), u5 = seq)
  GC <- sq(cds)  %>% transmute(transcript_id_clean, cds_gc = gc_content(seq), gc3 = calculate_gc3(seq), cd = seq)
  G3 <- sq(utr3) %>% transmute(transcript_id_clean, utr3_gc = gc_content(seq))
  out <- tibble(transcript_id_clean = want_tx) %>%
    left_join(L5, by = "transcript_id_clean") %>% left_join(L3, by = "transcript_id_clean") %>%
    left_join(LC, by = "transcript_id_clean") %>%
    left_join(G5, by = "transcript_id_clean") %>% left_join(G3, by = "transcript_id_clean") %>%
    left_join(GC, by = "transcript_id_clean") %>%
    mutate(log2_utr5_length = log2(utr5_length + 1), log2_utr3_length = log2(utr3_length + 1),
           log2_cds_length = log2(cds_length + 1),
           total_length = utr5_length + cds_length + utr3_length, log2_total_length = log2(total_length + 1),
           transcript_gc = (utr5_gc * utr5_length + cds_gc * cds_length + utr3_gc * utr3_length) /
                           (utr5_length + cds_length + utr3_length))
  kz <- out %>% filter(!is.na(u5), !is.na(cd)) %>%
    mutate(u5len = nchar(u5), win = paste0(str_sub(u5, start = pmax(1L, u5len - 5L)), str_sub(cd, 1, 4)),
           cat_ = vapply(win, match_kozak, character(1)),
           kozak_score = unname(score_map[cat_]), kozak_optimal = as.integer(kozak_score == 3L),
           kozak_pwm_score = vapply(win, score_pwm, numeric(1))) %>%
    dplyr::select(transcript_id_clean, kozak_score, kozak_optimal, kozak_pwm_score)
  out %>% dplyr::select(-u5, -cd) %>% left_join(kz, by = "transcript_id_clean") %>%
    left_join(tai_of(cds[names(cds) %in% want_tx]), by = "transcript_id_clean")
}
SEQ_COLS <- c("utr5_length", "utr3_length", "cds_length", "log2_utr5_length", "log2_utr3_length",
              "log2_cds_length", "total_length", "log2_total_length", "utr5_gc", "utr3_gc", "cds_gc",
              "transcript_gc", "gc3", "tai", "kozak_score", "kozak_optimal", "kozak_pwm_score")
gate(seq_set(PROBE), SEQ_COLS, 1e-9, "sequence features")
saveRDS(seq_set(tx), cache("seqfeat"))
cat("cached", CELL, "sequence features\n")
}

# =============================================================================================
# assemble: computed transcript columns + gene-level carry-over, in the 231 column order
# =============================================================================================
if (STAGE == "assemble") {
for (f in c("struct_from_lunp", "csc", "seqfeat"))
  stopifnot(setNames(list(file.exists(cache(f))), paste("run the stage that writes", cache(f)))[[1]])
st <- readRDS(cache("struct_from_lunp")); cs <- readRDS(cache("csc")); sf <- readRDS(cache("seqfeat"))
# A cache built for an earlier, smaller transcript set is skipped by 01m_ as "cached" and would
# leave the new transcripts silently NA after the left joins below (543 MANE transcripts, 2026-10-01).
for (cc in list(st = st, cs = cs, sf = sf))
  stopifnot("a stage cache does not cover this cell's transcript set - delete it and rerun that stage" =
              all(txm$transcript_id_clean %in% cc$transcript_id_clean))
computed <- c(setdiff(colnames(st), "transcript_id_clean"), setdiff(colnames(cs), "transcript_id_clean"),
              setdiff(colnames(sf), "transcript_id_clean"), "cai", "fop")
idc <- c("transcript_id", "transcript_id_clean", "gene_id", "gene_id_clean", "symbol")
carry <- setdiff(colnames(b), c(computed, idc))
by_gene <- b %>% dplyr::select(gene_id_clean, all_of(carry)) %>% group_by(gene_id_clean) %>% slice(1) %>% ungroup()
stopifnot("231 baseline is not one row per gene" = nrow(by_gene) == n_distinct(b$gene_id_clean))
anno <- read_tsv(here("accessories", "human", "gene_anno_hs_dm_v49_r111.tsv"),
                 col_names = c("gene_id", "a_symbol", "biotype", "species", "release"), show_col_types = FALSE) %>%
  transmute(gene_id, gene_id_clean = sub("\\..*", "", gene_id)) %>% distinct(gene_id_clean, .keep_all = TRUE)
# Genes whose ID is absent from the GENCODE v49 annotation are dropped (Kate, 2026-10-01: 13 HeLa
# genes from the Teleman transcript map, e.g. SOD2, SCO2, AKAP2, retired or re-IDed since the
# annotation that data was quantified on). They are listed so the loss is visible.
# Also, and more generally: the map's gene ID must be the gene that owns this transcript in v49.
# The HeLa map comes from a Salmon run on an older annotation, and 22 of its transcripts carry a
# stale gene ID - the 12 above (gene ID since re-assigned, e.g. SOD2 -> ENSG00000291237) plus 10
# whose old gene still exists but no longer owns the transcript (e.g. CYHR1). Outcomes join on
# gene, so such a row could pair one gene's measurement with another gene's transcript. Dropped.
.k <- AnnotationDbi::select(txdb, keys = keys(txdb, "TXNAME"), keytype = "TXNAME", columns = "GENEID")
v49_gene <- setNames(sub("[.].*", "", .k$GENEID), sub("[.].*", "", .k$TXNAME))
gone <- txm %>% filter(!gene_id_clean %in% anno$gene_id_clean |
                         is.na(v49_gene[transcript_id_clean]) | v49_gene[transcript_id_clean] != gene_id_clean)
cat(sprintf("dropped %d of %d genes absent from v49 or not owning their transcript there%s
", nrow(gone), nrow(txm),
            if (nrow(gone)) paste0(": ", paste(gone$symbol, collapse = ", ")) else ""))
stopifnot("more than 1% of genes are absent from v49 - this is not a few retired IDs" =
            nrow(gone) <= 0.01 * nrow(txm))
txm <- txm %>% filter(!transcript_id_clean %in% gone$transcript_id_clean)
out <- txm %>%
  left_join(anno, by = "gene_id_clean") %>%
  left_join(sf, by = "transcript_id_clean") %>% left_join(cs, by = "transcript_id_clean") %>%
  left_join(st, by = "transcript_id_clean") %>%
  mutate(cai = NA_real_, fop = NA_real_, transcript_id = transcript_id_clean) %>%
  left_join(by_gene, by = "gene_id_clean")
cat(sprintf("carried by gene from the 231 baseline: %d columns | %s genes present in it: %.1f%%\n",
            length(carry), CELL, 100 * mean(txm$gene_id_clean %in% b$gene_id_clean)))
miss <- setdiff(colnames(b), colnames(out)); if (length(miss)) cat("MISSING:", paste(miss, collapse = ", "), "\n")
out <- out %>% dplyr::select(all_of(colnames(b)))
stopifnot("does not carry exactly the 231 baseline's columns" = identical(colnames(out), colnames(b)),
          "transcript_id_clean not unique" = !any(duplicated(out$transcript_id_clean)),
          "gene_id missing for some genes" = !any(is.na(out$gene_id)),
          "accessibility outside [0,1]" = all(out$struct_accessibility_cds_mean >= 0 & out$struct_accessibility_cds_mean <= 1, na.rm = TRUE),
          "GC outside [0,100]" = all(out$cds_gc >= 0 & out$cds_gc <= 100, na.rm = TRUE),
          "tAI outside [0,1]" = all(out$tai >= 0 & out$tai <= 1, na.rm = TRUE),
          "Kozak tier outside 0-3" = all(out$kozak_score %in% c(0:3, NA)))
fn <- paste0("feature_matrix_", CELL, "_dhx29_kd_feature.rds")
saveRDS(out, file.path(PD, fn))
cat(sprintf("wrote %s - %d x %d (231 baseline: %d x %d)\n", fn, nrow(out), ncol(out), nrow(b), ncol(b)))
}

# =============================================================================================
# kozak: the sequence cache 01_ expects (as 01i_)
# =============================================================================================
if (STAGE == "kozak") {
utr5 <- grab(fiveUTRsByTranscript(txdb, use.names = TRUE)); cds <- grab(cdsBy(txdb, by = "tx", use.names = TRUE))
ref_seqs <- readRDS(file.path(PD, "precomputed_kozak_sequences.rds"))
# Follow the assembled baseline: assemble drops genes absent from v49 (12 for HeLa), so the
# cache must cover exactly its transcripts, which must all be in this cell's universe.
bl <- file.path(PD, paste0("feature_matrix_", CELL, "_dhx29_kd_feature.rds"))
if (file.exists(bl)) {
  bl_tx <- readRDS(bl)$transcript_id_clean
  stopifnot("baseline has transcripts outside this cell's universe" = all(bl_tx %in% tx))
  tx <- bl_tx
}
u5col <- setdiff(colnames(ref_seqs$utr5_seq), "transcript_id_clean")
cdcol <- setdiff(colnames(ref_seqs$cds_seq),  "transcript_id_clean")
mk <- function(sset, col) { s <- sset[names(sset) %in% tx]; tibble(transcript_id_clean = names(s), !!col := as.character(s)) }
seqs <- list(utr5_seq = mk(utr5, u5col), cds_seq = mk(cds, cdcol))
stopifnot("Kozak cache does not cover every transcript" = nrow(seqs$utr5_seq) == length(tx) && nrow(seqs$cds_seq) == length(tx))
bl <- file.path(PD, paste0("feature_matrix_", CELL, "_dhx29_kd_feature.rds"))
if (file.exists(bl)) stopifnot("Kozak cache and baseline disagree" = setequal(tx, readRDS(bl)$transcript_id_clean))
saveRDS(seqs, file.path(PD, paste0("precomputed_", CELL, "_kozak_sequences.rds")))
cat(sprintf("wrote precomputed_%s_kozak_sequences.rds\n", CELL))
}
