# 33b_sequence_features.R - compute the model's intrinsic features from SEQUENCE ALONE.
#
# The feature matrix is built by joining columns that upstream notebooks computed from a
# transcript annotation. A plasmid has no annotation, so its features have to be computed from
# the sequence. This file re-implements the 33 sequence-only features of the 57-feature intrinsic
# model, following the definitions in 01_build_feature_matrix.Rmd exactly. Structure (15) needs
# RNAplfold and G4 (9) needs the G4mer model; they live in 33c_ and 33d_.
#
# EVERY definition here is reproduced, not reinvented, and 33e_ asserts that: it runs these
# functions on real transcripts and requires the output to match the stored matrix. A
# re-implementation that is not checked that way is just a second opinion.
#
# Reference tables are imported, never transcribed (accessories/collaborator_reference/), and
# the codon weights are LOADED from codon_weights_mdamb231.rds rather than refitted - refitting
# on 4 plasmids would make the features mean something different from the models' training data.
#
# Usage:  source("code/translational_control_overview/33b_sequence_features.R")
#         f <- seq_features(tibble(id = ..., utr5 = ..., cds = ..., utr3 = ...))
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(Biostrings); library(coRdon); library(iCodon)
})

REF_DIR <- here("accessories", "collaborator_reference")
OUT_DIR <- here("output", "translational_control_overview")

GEN_CODE     <- Biostrings::getGeneticCode("SGC0")
STOPS        <- names(GEN_CODE[GEN_CODE == "*"])
SENSE_CODONS <- names(GEN_CODE[GEN_CODE != "*"])
AA_TO_CODONS <- split(SENSE_CODONS, GEN_CODE[SENSE_CODONS])

codons_of <- function(s) {
  n <- nchar(s) %/% 3L
  if (n < 1L) return(character(0))
  substring(s, seq.int(1L, n * 3L, 3L), seq.int(3L, n * 3L, 3L))
}
pct_gc <- function(s) if (!nchar(s)) NA_real_ else
  100 * mean(strsplit(s, "", fixed = TRUE)[[1]] %in% c("G", "C"))

# ---- reference tables, loaded once --------------------------------------------------
.kozak_logodds <- local({
  k <- read_tsv(file.path(REF_DIR, "kozak_1987_table1_frequencies.tsv"), show_col_types = FALSE)
  win <- c("-6" = 1, "-5" = 2, "-4" = 3, "-3" = 4, "-2" = 5, "-1" = 6, "+4" = 10)
  stopifnot("unexpected Kozak context_index set" = setequal(k$context_index, c(1:6, 10)))
  m <- k %>% arrange(match(context_index, win)) %>%
    dplyr::select(A = A_frequency, C = C_frequency, G = G_frequency, T = T_frequency) %>% as.matrix()
  rownames(m) <- names(win)
  stopifnot("published Kozak rows do not sum to 1" = all(abs(rowSums(m) - 1) < 1e-8))
  list(lo = log2((m + 0.001) / 0.25), win = win)
})

.noderer <- local({
  n <- read_tsv(file.path(REF_DIR, "noderer_2014_tis_efficiency.tsv"), skip = 1, show_col_types = FALSE) %>%
    mutate(context = str_replace_all(sequence, "U", "T"))
  stopifnot("Noderer table is not 65,536 contexts"  = nrow(n) == 65536L,
            "Noderer contexts are not AUG-anchored" = all(str_sub(n$context, 7, 9) == "ATG"),
            "Noderer efficiency outside 12-150"     = all(between(n$efficiency, 12, 150)))
  setNames(n$efficiency, n$context)
})

.cridge <- local({
  cr <- read_tsv(file.path(REF_DIR, "cridge2018_termination_contexts.tsv"), show_col_types = FALSE)
  stopifnot("Cridge table is not 192 contexts" = nrow(cr) == 192L,
            "Cridge median readthrough off"    = between(median(cr$readthrough_measure), 0.08, 0.095))
  setNames(cr$readthrough_measure, cr$termination_context_plus1_to_plus6)
})

.tai_w <- local({
  gcn <- read_tsv(file.path(REF_DIR, "human_trna_gene_copy_hg38.tsv"), comment = "#", show_col_types = FALSE)
  wob <- read_tsv(file.path(REF_DIR, "trna_wobble_scheme_v1.tsv"), show_col_types = FALSE) %>%
    transmute(w = chartr("U", "T", anticodon_wobble_base_rna),
              n3 = chartr("U", "T", codon_third_base_rna), s = penalty_s)
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  W <- map_dbl(SENSE_CODONS, function(cod) {
    r <- wob %>% filter(n3 == substr(cod, 3, 3))
    if (!nrow(r)) return(0)
    anti <- paste0(r$w, comp[[substr(cod, 2, 2)]], comp[[substr(cod, 1, 1)]])
    cn <- gcn$gene_copy_number[match(anti, gcn$anticodon)]; cn[is.na(cn)] <- 0
    sum((1 - r$s) * cn)
  })
  names(W) <- SENSE_CODONS
  tw <- W / max(W)
  if (any(tw == 0)) tw[tw == 0] <- exp(mean(log(tw[tw > 0])))   # dos Reis substitution
  tw
})

.codon_w <- local({
  p <- file.path(OUT_DIR, "codon_weights_mdamb231.rds")
  stopifnot("codon_weights_mdamb231.rds not found - knit 01_build_feature_matrix.Rmd first" = file.exists(p))
  cw <- readRDS(p)
  stopifnot("the codon weights file is not the verified object" =
              all(c("w_int", "opt_int", "cai_dev", "fop_dev") %in% names(cw)),
            "the loaded codon weights were not verified against 01g_" =
              cw$cai_dev < 1e-9 && cw$fop_dev < 1e-9)
  cw
})

# ---- the feature functions ----------------------------------------------------------
# Initiator AUG and terminal stop removed; neither is decoded as a sense codon.
trim_internal <- function(x) {
  n <- nchar(x) %/% 3L
  a <- ifelse(substr(x, 1, 3) == "ATG", 4L, 1L)
  b <- ifelse(substr(x, (n - 1) * 3L + 1L, n * 3L) %in% STOPS, (n - 1) * 3L, n * 3L)
  substr(x, a, b)
}

cai_of <- function(cm, w) {
  use <- SENSE_CODONS[!is.na(w[SENSE_CODONS]) & w[SENSE_CODONS] > 0]
  apply(cm, 1, function(g) { cu <- g[use]; n <- sum(cu)
    if (n == 0) NA_real_ else exp(sum(cu * log(w[use])) / n) })
}
fop_of <- function(cm, opt) {
  sc <- setdiff(colnames(cm), STOPS)
  tot <- rowSums(cm[, sc, drop = FALSE]); o <- rowSums(cm[, opt, drop = FALSE])
  r <- o / tot; r[tot == 0] <- NA_real_; r
}
gc3_of <- function(s) {
  cd <- codons_of(s); if (!length(cd)) return(NA_real_)
  100 * mean(substr(cd, 3, 3) %in% c("G", "C"))
}
# K/R +1, D/E -1, maximised over any 30-aa window; short chains score their whole length.
charge_window_max <- function(aa, w = 30L) {
  ch <- strsplit(aa, "", fixed = TRUE)[[1]]
  if (!length(ch)) return(NA_real_)
  v <- integer(length(ch)); v[ch %in% c("K","R")] <- 1L; v[ch %in% c("D","E")] <- -1L
  if (length(v) <= w) return(as.numeric(sum(v)))
  cs <- cumsum(c(0L, v))
  as.numeric(max(cs[(w + 1L):length(cs)] - cs[seq_len(length(cs) - w)]))
}
max_codon_run <- function(s, codon) {
  cd <- codons_of(s); if (!length(cd)) return(0L)
  r <- rle(cd == codon)
  if (!any(r$values)) 0L else max(r$lengths[r$values])
}
next_inframe_stop <- function(u) {
  cd <- codons_of(u); if (!length(cd)) return(NA_integer_)
  w <- which(cd %in% STOPS)
  if (!length(w)) NA_integer_ else w[1]
}

#' Compute the 33 sequence-only intrinsic features.
#' @param d tibble with columns id, utr5, cds, utr3 (uppercase ACGT)
seq_features <- function(d) {
  stopifnot(all(c("id", "utr5", "cds", "utr3") %in% names(d)))
  d <- d %>% mutate(across(c(utr5, cds, utr3), toupper))
  n <- nchar(d$cds)

  # CDS eligibility, exactly as the builder defines it. The mask that matters is a MISSING
  # TERMINAL STOP, not a bad frame: an in-frame but truncated CDS passes every other filter.
  has_stop  <- substr(d$cds, (n %/% 3L - 1L) * 3L + 1L, (n %/% 3L) * 3L) %in% STOPS
  frame_ok  <- n %% 3L == 0L
  truncated <- !has_stop
  stopifnot("a CDS is not in frame - codon features would be computed on a broken frame" = all(frame_ok))

  s_int  <- trim_internal(d$cds)
  int_ss <- DNAStringSet(s_int); names(int_ss) <- d$id
  cm_int <- codonCounts(codonTable(int_ss)); rownames(cm_int) <- d$id
  pep_aa <- Biostrings::translate(DNAStringSet(s_int), if.fuzzy.codon = "solve")
  pep_len <- nchar(pep_aa); aa_chr <- as.character(pep_aa)

  # Start context. -6..-1 from the 5'UTR, +1..+4 (Kozak) or +1..+5 (Noderer) from the CDS.
  u5n   <- nchar(d$utr5)
  win10 <- paste0(str_sub(d$utr5, start = pmax(1L, u5n - 5L)), str_sub(d$cds, 1, 4))
  win11 <- paste0(str_sub(d$utr5, start = pmax(1L, u5n - 5L)), str_sub(d$cds, 1, 5))
  ok10  <- nchar(win10) == 10L & str_sub(win10, 7, 9) == "ATG"
  ok11  <- nchar(win11) == 11L & str_sub(win11, 7, 9) == "ATG"
  match_kozak <- function(s) {
    if (nchar(s) < 10) return(NA_integer_)
    if (grepl("^(GCCACC|GCCGCC)ATGG$", s)) return(3L)
    if (grepl("^...[AG]..ATGG$", s))       return(2L)
    if (grepl("^...[AG]..ATG[ACT]$", s))   return(1L)
    if (grepl("^...[CT]..ATGG$", s))       return(1L)
    0L
  }
  pwm <- .kozak_logodds
  pwm_score <- rep(NA_real_, nrow(d))
  if (any(ok10)) {
    mat <- str_split_fixed(win10[ok10], "", 10)
    pwm_score[ok10] <- rowSums(sapply(seq_along(pwm$win), function(k)
      pwm$lo[names(pwm$win)[k], mat[, pwm$win[k]]]))
  }
  kz <- vapply(win10, match_kozak, integer(1), USE.NAMES = FALSE)

  # Termination context: the stop codon plus the first bases of the 3'UTR.
  stop_cod <- substr(d$cds, (n %/% 3L - 1L) * 3L + 1L, (n %/% 3L) * 3L)
  ctx6 <- paste0(stop_cod, str_sub(d$utr3, 1, 3))
  ctx9 <- paste0(stop_cod, str_sub(d$utr3, 1, 6))
  p4   <- str_sub(d$utr3, 1, 1)
  rt   <- unname(.cridge[ctx6])

  mask <- function(x) if_else(truncated, NA_real_, as.numeric(x))
  out <- tibble(
    id = d$id,
    # +1 pseudocount, as 01_feature_extraction.Rmd defines these (not log2(n)).
    log2_utr5_length  = log2(nchar(d$utr5) + 1),
    log2_utr3_length  = log2(nchar(d$utr3) + 1),
    log2_cds_length   = log2(n + 1),
    log2_total_length = log2(nchar(d$utr5) + n + nchar(d$utr3) + 1),
    utr5_gc       = map_dbl(d$utr5, pct_gc),
    utr3_gc       = map_dbl(d$utr3, pct_gc),
    cds_gc        = map_dbl(d$cds, pct_gc),
    transcript_gc = map_dbl(paste0(d$utr5, d$cds, d$utr3), pct_gc),
    gc3_internal  = mask(map_dbl(s_int, gc3_of)),
    kozak_score   = as.numeric(kz),
    kozak_optimal = as.numeric(kz == 3L),
    kozak_pwm_score_v2    = pwm_score,
    noderer_tis_efficiency = if_else(ok11, unname(.noderer[win11]), NA_real_),
    tai_gtrnadb  = mask(map_dbl(d$cds, function(sq) {
      cd <- codons_of(sq); if (length(cd) < 3L) return(NA_real_)
      cd <- cd[-1]; if (tail(cd, 1) %in% STOPS) cd <- head(cd, -1)
      ww <- .tai_w[cd]; ww <- ww[!is.na(ww)]
      if (!length(ww)) NA_real_ else exp(mean(log(ww)))
    })),
    cai_internal = mask(cai_of(cm_int, .codon_w$w_int)),
    fop_internal = mask(fop_of(cm_int, .codon_w$opt_int)),
    csc_internal = mask(as.numeric(predict_stability("human")(s_int))),
    proline_fraction    = mask(as.numeric(letterFrequency(pep_aa, "P", as.prob = TRUE))),
    ppp_motif_density   = mask(if_else(pep_len > 0, 100 * vcountPattern("PPP", pep_aa) / pep_len, NA_real_)),
    max_net_charge_30aa = mask(vapply(aa_chr, charge_window_max, numeric(1), USE.NAMES = FALSE)),
    max_consecutive_aaa = mask(vapply(s_int, max_codon_run, integer(1), "AAA", USE.NAMES = FALSE)),
    max_consecutive_aag = mask(vapply(s_int, max_codon_run, integer(1), "AAG", USE.NAMES = FALSE)),
    term_stop_TAA = as.numeric(stop_cod == "TAA"),
    term_stop_TAG = as.numeric(stop_cod == "TAG"),
    term_stop_TGA = as.numeric(stop_cod == "TGA"),
    term_plus4_A  = as.numeric(p4 == "A"), term_plus4_C = as.numeric(p4 == "C"),
    term_plus4_G  = as.numeric(p4 == "G"), term_plus4_T = as.numeric(p4 == "T"),
    term_window_gc = if_else(nchar(ctx9) == 9L, 100 * str_count(ctx9, "[GC]") / 9, NA_real_),
    term_interval_codons   = as.numeric(map_int(d$utr3, next_inframe_stop)),
    term_readthrough_cridge = rt,
    term_log1p_readthrough  = log1p(rt))
  # The termination family is defined only where the CDS ends in a stop.
  tcols <- grep("^term_", names(out), value = TRUE)
  out[truncated, tcols] <- NA_real_

  # Range invariants that must hold by definition.
  stopifnot(
    "GC outside [0,100]"   = all(map_lgl(out[c("utr5_gc","utr3_gc","cds_gc","transcript_gc","gc3_internal","term_window_gc")],
                                         ~ all(.x >= 0 & .x <= 100, na.rm = TRUE))),
    "CAI/tAI/FOP outside [0,1]" = all(map_lgl(out[c("cai_internal","fop_internal","tai_gtrnadb")],
                                              ~ all(.x >= 0 & .x <= 1, na.rm = TRUE))),
    "kozak_score outside {0,1,2,3}" = all(out$kozak_score %in% c(0:3, NA)),
    "proline fraction outside [0,1]" = all(out$proline_fraction >= 0 & out$proline_fraction <= 1, na.rm = TRUE),
    "a stop-codon indicator set is not one-hot" =
      all(rowSums(out[c("term_stop_TAA","term_stop_TAG","term_stop_TGA")], na.rm = TRUE) %in% 0:1),
    "lengths are not finite" = all(map_lgl(out[grep("^log2_", names(out))], ~ all(is.finite(.x)))))
  out
}
