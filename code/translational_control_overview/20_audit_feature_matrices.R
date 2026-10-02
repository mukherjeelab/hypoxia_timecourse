# 20_audit_feature_matrices.R
# Independent audit of every INTRINSIC modelled feature in every feature matrix (Kate, 2026-10-01:
# "check all the transcripts to make sure the calculations are correct").
#
# Independence is the point. Each feature is recomputed here from PRIMARY sources - the GENCODE v49
# TxDb + hg38 genome (not the cached sequence files the builders read), the RNAplfold _lunp files
# (own parser), the G4mer summary files, and the published tables in
# accessories/collaborator_reference/ - with code written for this audit rather than copied from
# 01_ / 01l_. Definitions follow FEATURES.md. Shared constants that ARE the definition (the 231 codon
# weights, the 231 G4 residual coefficients, the CNOT3 P-site weights) are read from their saved
# files, which were themselves verified against their sources when written.
#
# Outputs, in output/translational_control_overview/:
#   audit_feature_matrices_summary.csv   one row per matrix x feature: n compared, value mismatches,
#                                         NA only in matrix, NA only in recompute, max |dev|
#   audit_feature_matrices_examples.csv  up to 5 offending transcripts per failing matrix x feature
#   audit_cross_matrix.csv               same transcript, different matrix: any disagreement
# csc_internal (iCodon) is recomputed on a random sample (n_csc per matrix), the only sampled check.
#   Rscript code/translational_control_overview/20_audit_feature_matrices.R [n_csc]
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(stringr); library(here)
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38); library(Biostrings)
  library(iCodon)   # attached, not just namespaced: predict_stability() looks up package data by name
})
args  <- commandArgs(TRUE); N_CSC <- if (length(args)) as.integer(args[1]) else 200L
O     <- here("output", "translational_control_overview")
PD    <- here("output", "predictive_modeling")
REF   <- here("accessories", "collaborator_reference")
PLF   <- here("accessories", "plfold_output", "rnaplfold_output")
TOL   <- 1e-9
MATS  <- c(mdamb231 = "feature_matrix_tco.rds", mcf7six1 = "feature_matrix_tco_mcf7six1.rds",
           hela = "feature_matrix_tco_hela.rds", dhx29 = "feature_matrix_tco_dhx29.rds",
           mane = "feature_matrix_tco_mane.rds")
MATS  <- MATS[file.exists(file.path(O, MATS))]
if (nzchar(Sys.getenv("AUDIT_CELLS"))) MATS <- MATS[strsplit(Sys.getenv("AUDIT_CELLS"), ",")[[1]]]
G4TAG <- c(mdamb231 = "_mdamb231", mcf7six1 = "_mcf7six1", hela = "_hela", dhx29 = "_dhx29", mane = "_mane")
intr  <- read_csv(file.path(O, "rfreg_importance_eif33d_promotes_hypoxia_1hr_lfc0.5_tcoreg_intrinsic_corrected_both_clipaggregate_yraw.csv"),
                  show_col_types = FALSE)$feature
mats <- map(MATS, ~ readRDS(file.path(O, .x)))
cat("matrices:", paste(sprintf("%s (%d)", names(mats), map_int(mats, nrow)), collapse = ", "), "\n")
tx_all <- unique(unlist(map(mats, "transcript_id_clean")))

# ---- 1. sequences straight from GENCODE v49 + hg38 ----------------------------------------------
txdb <- loadDb(here("accessories", "human", "txdb.gencode49.sqlite"))
ext <- function(grl) {
  names(grl) <- sub("[.].*", "", names(grl)); grl <- grl[names(grl) %in% tx_all]
  s <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, grl); setNames(toupper(as.character(s)), names(s))
}
U5 <- ext(fiveUTRsByTranscript(txdb, use.names = TRUE))
CD <- ext(cdsBy(txdb, by = "tx", use.names = TRUE))
U3 <- ext(threeUTRsByTranscript(txdb, use.names = TRUE))
cat("sequences from v49: 5'UTR", length(U5), "| CDS", length(CD), "| 3'UTR", length(U3), "of", length(tx_all), "\n")

gcp  <- function(s) ifelse(nchar(s) > 0, 100 * str_count(s, "[GC]") / nchar(s), NA_real_)
cods <- function(s) { n <- nchar(s) %/% 3L; if (n < 1) character(0) else substring(s, seq(1, 3 * n, 3), seq(3, 3 * n, 3)) }
STOPS <- c("TAA", "TAG", "TGA")
CODE <- GENETIC_CODE   # Biostrings standard code, used only as a lookup table

# ---- 2. constants that are part of the definitions ----------------------------------------------
cw  <- readRDS(file.path(O, "codon_weights_mdamb231.rds")); w_int <- cw$w_int; opt_int <- cw$opt_int
g4c <- readRDS(file.path(O, "g4_resid_coefs_mdamb231.rds"))
cn3 <- readRDS(file.path(O, "cnot3_psite_weights.rds"))
kz  <- read_tsv(file.path(REF, "kozak_1987_table1_frequencies.tsv"), show_col_types = FALSE)
pos <- c(`1` = 1, `2` = 2, `3` = 3, `4` = 4, `5` = 5, `6` = 6, `10` = 10)   # window index per row
kz_lo <- kz %>% transmute(idx = context_index, A = log2((A_frequency + 0.001) / 0.25),
                          C = log2((C_frequency + 0.001) / 0.25), G = log2((G_frequency + 0.001) / 0.25),
                          T = log2((T_frequency + 0.001) / 0.25))
nod <- read_tsv(file.path(REF, "noderer_2014_tis_efficiency.tsv"), skip = 1, show_col_types = FALSE) %>%
  transmute(ctx = chartr("U", "T", sequence), eff = efficiency)
nod_v <- setNames(nod$eff, nod$ctx)
crd <- read_tsv(file.path(REF, "cridge2018_termination_contexts.tsv"), show_col_types = FALSE)
crd_v <- setNames(crd$readthrough_measure, crd$termination_context_plus1_to_plus6)
# tAI weights from the GtRNAdb copy-number table and the collaborator wobble scheme (dos Reis)
gcn <- read_tsv(file.path(REF, "human_trna_gene_copy_hg38.tsv"), comment = "#", show_col_types = FALSE)
wob <- read_tsv(file.path(REF, "trna_wobble_scheme_v1.tsv"), show_col_types = FALSE)
cmp <- c(A = "T", C = "G", G = "C", T = "A")
sense <- setdiff(names(CODE), STOPS)
tW <- vapply(sense, function(cod) {
  r <- wob[chartr("U", "T", wob$codon_third_base_rna) == substr(cod, 3, 3), ]
  if (!nrow(r)) return(0)
  ac <- paste0(chartr("U", "T", r$anticodon_wobble_base_rna), cmp[substr(cod, 2, 2)], cmp[substr(cod, 1, 1)])
  cn <- gcn$gene_copy_number[match(ac, gcn$anticodon)]; cn[is.na(cn)] <- 0
  sum((1 - r$penalty_s) * cn)
}, numeric(1))
tW <- tW / max(tW); tW[tW == 0] <- exp(mean(log(tW[tW > 0])))

# ---- 3. per-transcript recomputation (once per transcript, shared by every matrix) --------------
cat("recomputing sequence features...\n")
one <- function(t) {
  u5 <- U5[t]; cd <- CD[t]; u3 <- U3[t]
  u5 <- if (is.na(u5)) "" else u5; u3 <- if (is.na(u3)) "" else u3; has_cd <- !is.na(cd)
  r <- list(transcript_id_clean = t,
            log2_utr5_length = if (nchar(u5)) log2(nchar(u5) + 1) else NA_real_,
            log2_cds_length  = if (has_cd) log2(nchar(cd) + 1) else NA_real_,
            log2_utr3_length = if (nchar(u3)) log2(nchar(u3) + 1) else NA_real_,
            utr5_gc = if (nchar(u5)) gcp(u5) else NA_real_, cds_gc = if (has_cd) gcp(cd) else NA_real_,
            utr3_gc = if (nchar(u3)) gcp(u3) else NA_real_)
  r$log2_total_length <- if (nchar(u5) && has_cd && nchar(u3)) log2(nchar(u5) + nchar(cd) + nchar(u3) + 1) else NA_real_
  r$transcript_gc <- if (nchar(u5) && has_cd && nchar(u3))
    100 * (str_count(u5, "[GC]") + str_count(cd, "[GC]") + str_count(u3, "[GC]")) / (nchar(u5) + nchar(cd) + nchar(u3)) else NA_real_
  # Kozak: last 6 nt of 5'UTR + first 4 of CDS (tiers, PWM), -6..+5 (Noderer)
  w10 <- paste0(str_sub(u5, -6), str_sub(cd, 1, 4)); w11 <- paste0(str_sub(u5, -6), str_sub(cd, 1, 5))
  r$kozak_score <- if (!has_cd || nchar(w10) < 10) NA_integer_ else
    if (grepl("^(GCCACC|GCCGCC)ATGG$", w10)) 3L else if (grepl("^...[AG]..ATGG$", w10)) 2L else
    if (grepl("^...[AG]..ATG[ACT]$", w10) || grepl("^...[CT]..ATGG$", w10)) 1L else 0L
  r$kozak_optimal <- as.integer(r$kozak_score == 3L)
  r$kozak_pwm_score_v2 <- if (has_cd && nchar(w10) == 10 && substr(w10, 7, 9) == "ATG" && !grepl("[^ACGT]", w10))
    sum(vapply(seq_len(nrow(kz_lo)), function(k) kz_lo[[substr(w10, kz_lo$idx[k], kz_lo$idx[k])]][k], numeric(1))) else NA_real_
  r$noderer_tis_efficiency <- if (has_cd && nchar(w11) == 11 && substr(w11, 7, 9) == "ATG") unname(nod_v[w11]) else NA_real_
  # codon-level: complete CDS = frame ok and a terminal stop; internal = drop ATG initiator + stop
  cc <- if (has_cd) cods(cd) else character(0)
  frame_ok <- has_cd && nchar(cd) %% 3 == 0
  has_stop <- length(cc) > 0 && tail(cc, 1) %in% STOPS
  na_codon <- c("gc3_internal", "cai_internal", "fop_internal", "tai_gtrnadb", "proline_fraction",
                "ppp_motif_density", "max_net_charge_30aa", "max_consecutive_aaa", "max_consecutive_aag")
  if (!(frame_ok && has_stop)) { for (f in na_codon) r[[f]] <- NA_real_ } else {
    ic <- cc; if (ic[1] == "ATG") ic <- ic[-1]; ic <- head(ic, -1)
    r$gc3_internal <- if (length(ic)) 100 * mean(substr(ic, 3, 3) %in% c("G", "C")) else NA_real_
    ws <- w_int[ic]; ws <- ws[!is.na(ws) & ws > 0]
    r$cai_internal <- if (length(ws)) exp(mean(log(ws))) else NA_real_
    # FOP's denominator is sense codons only: selenoproteins carry internal UGA (selenocysteine),
    # which 01g_/01_ exclude as a stop column (found by this audit: 18 SELENO* transcripts).
    sc_ic <- ic[!ic %in% STOPS]
    r$fop_internal <- if (length(sc_ic)) mean(sc_ic %in% opt_int) else NA_real_
    tc <- head(cc[-1], -1); tw <- tW[tc]; tw <- tw[!is.na(tw)]
    r$tai_gtrnadb <- if (length(cc) >= 3 && length(tw)) exp(mean(log(tw))) else NA_real_
    aa <- unname(ifelse(ic %in% names(CODE), CODE[ic], "X"))
    r$proline_fraction <- if (length(aa)) mean(aa == "P") else NA_real_
    s <- paste(aa, collapse = "")
    r$ppp_motif_density <- if (length(aa)) 100 * lengths(gregexpr("(?=PPP)", s, perl = TRUE)) * grepl("PPP", s) / length(aa) else NA_real_
    v <- ifelse(aa %in% c("K", "R"), 1, ifelse(aa %in% c("D", "E"), -1, 0))
    r$max_net_charge_30aa <- if (!length(v)) NA_real_ else if (length(v) <= 30) sum(v) else
      max(stats::filter(v, rep(1, 30), sides = 1)[30:length(v)])
    run <- function(cod) { x <- rle(ic == cod); if (any(x$values)) max(x$lengths[x$values]) else 0 }
    r$max_consecutive_aaa <- run("AAA"); r$max_consecutive_aag <- run("AAG")
  }
  # termination: needs a terminal stop and a 3'UTR
  tcols <- c("term_stop_TAA", "term_stop_TAG", "term_stop_TGA", "term_plus4_A", "term_plus4_C", "term_plus4_G",
             "term_plus4_T", "term_window_gc", "term_interval_codons", "term_readthrough_cridge", "term_log1p_readthrough")
  if (has_cd && has_stop && nchar(u3)) {
    sc <- tail(cc, 1); p4 <- substr(u3, 1, 1)
    for (k in STOPS) r[[paste0("term_stop_", k)]] <- as.numeric(sc == k)
    for (b in c("A", "C", "G", "T")) r[[paste0("term_plus4_", b)]] <- as.numeric(p4 == b)
    c9 <- paste0(sc, substr(u3, 1, 6))
    r$term_window_gc <- if (nchar(c9) == 9) 100 * str_count(c9, "[GC]") / 9 else NA_real_
    uc <- cods(u3); hit <- which(uc %in% STOPS)
    r$term_interval_codons <- if (length(hit)) hit[1] else NA_real_
    rt <- unname(crd_v[paste0(sc, substr(u3, 1, 3))])
    r$term_readthrough_cridge <- if (length(rt)) rt else NA_real_
    r$term_log1p_readthrough <- log1p(r$term_readthrough_cridge)
  } else for (f in tcols) r[[f]] <- NA_real_
  r$cnot3_weighted_codon_score_tx <- if (has_cd && frame_ok && length(cc) >= cn3$min_codons)
    sum(cn3$weights[cc], na.rm = TRUE) / length(cc) else NA_real_
  as_tibble(map(r, ~ if (is.null(.x) || !length(.x)) NA else .x))
}
seqf <- map_dfr(tx_all, one)

# ---- 4. structure: own parser of the RNAplfold _lunp files --------------------------------------
cat("parsing _lunp files...\n")
read_lunp <- function(f) {
  if (!file.exists(f)) return(NULL)
  l <- readLines(f, warn = FALSE); l <- l[!grepl("^\\s*#", l) & nzchar(trimws(l))]
  if (!length(l)) return(NULL)
  v <- suppressWarnings(as.numeric(vapply(strsplit(l, "\t", fixed = TRUE), `[`, "", 2)))
  v[!is.na(v)]
}
SDIR <- c(utr5 = "rnaplfold_output", cds = "cds_rnaplfold_output", utr3 = "utr3_rnaplfold_output")
SNAM <- list(utr5 = c("cap_proximal", "aug_context", "utr5_mean", "utr5_min", "utr5_sd"),
             cds  = c("start_proximal_cds", "stop_proximal_cds", "cds_mean", "cds_min", "cds_sd"),
             utr3 = c("stop_proximal_utr3", "distal_utr3", "utr3_mean", "utr3_min", "utr3_sd"))
struc <- map_dfr(tx_all, function(t) {
  out <- list(transcript_id_clean = t)
  for (rg in names(SDIR)) {
    a <- read_lunp(file.path(PLF, SDIR[[rg]], paste0(t, "_lunp"))); nm <- paste0("tco_struct_accessibility_", SNAM[[rg]])
    vals <- if (is.null(a) || !length(a)) rep(NA_real_, 5) else
      c(mean(a[1:min(30, length(a))]), mean(a[max(1, length(a) - 29):length(a)]), mean(a), min(a), sd(a))
    out[nm] <- as.list(vals)
  }
  as_tibble(out)
})

# ---- 5. compare --------------------------------------------------------------------------------
cmp_one <- function(cell, m) {
  ref <- seqf %>% inner_join(struc, by = "transcript_id_clean")
  # G4: this cell's own summaries; residual with the shared 231 coefficients
  for (rg in c("utr5", "cds", "utr3")) {
    g <- read_tsv(file.path(here("output", "g4mer"), paste0("g4mer_summary_", rg, G4TAG[[cell]], ".tsv")),
                  show_col_types = FALSE)
    i <- match(ref$transcript_id_clean, g$transcript_id_clean)
    ref[[paste0("g4mer_mean_", rg)]] <- g$g4mer_mean[i]
    ref[[paste0("g4mer_frac_above_", rg)]] <- g$g4mer_frac_above[i]
    cf <- g4c[[rg]]; gcv <- m[[paste0(rg, "_gc")]][match(ref$transcript_id_clean, m$transcript_id_clean)]
    ref[[paste0("g4mer_max_resid_", rg)]] <- g$g4mer_max[i] - (cf[1] + cf[2] * log2(g$region_length[i] + 1) + cf[3] * gcv)
  }
  feats <- intersect(c(intr, "cnot3_weighted_codon_score_tx"), names(ref))
  feats <- setdiff(feats, "csc_internal")
  x <- m[match(ref$transcript_id_clean, m$transcript_id_clean), ]; keep <- !is.na(x$transcript_id_clean)
  x <- x[keep, ]; ref <- ref[keep, ]
  res <- map_dfr(feats, function(f) {
    u <- as.numeric(unlist(x[[f]])); v <- as.numeric(ref[[f]]); k <- !is.na(u) & !is.na(v)
    bad <- which((k & abs(u - v) > TOL) | (is.na(u) != is.na(v)))
    tibble(matrix = cell, feature = f, n_compared = sum(k), value_mismatch = sum(k & abs(u - v) > TOL),
           na_only_in_matrix = sum(is.na(u) & !is.na(v)), na_only_in_recompute = sum(!is.na(u) & is.na(v)),
           max_abs_dev = if (any(k)) max(abs(u[k] - v[k])) else NA_real_,
           examples = paste(head(x$transcript_id_clean[bad], 5), collapse = " "))
  })
  # csc_internal: iCodon on a random sample of complete CDS
  set.seed(9)
  cand <- x$transcript_id_clean[!is.na(x$csc_internal)]; smp <- sample(cand, min(N_CSC, length(cand)))
  ints <- vapply(smp, function(t) { cc <- cods(CD[t]); if (cc[1] == "ATG") cc <- cc[-1]; paste(head(cc, -1), collapse = "") }, "")
  pr <- tryCatch(suppressWarnings(as.numeric(predict_stability("human")(unname(ints)))),
                 error = function(e) { cat("CSC sample FAILED:", conditionMessage(e), "\n"); rep(NA_real_, length(ints)) })
  st <- x$csc_internal[match(smp, x$transcript_id_clean)]
  bind_rows(res, tibble(matrix = cell, feature = "csc_internal (sample)", n_compared = length(smp),
                        value_mismatch = sum(abs(pr - st) > TOL, na.rm = TRUE), na_only_in_matrix = 0L,
                        na_only_in_recompute = sum(is.na(pr)), max_abs_dev = max(abs(pr - st), na.rm = TRUE),
                        examples = paste(head(smp[abs(pr - st) > TOL], 5), collapse = " ")))
}
summ <- imap_dfr(mats, ~ { cat("comparing", .y, "\n"); cmp_one(.y, .x) })
write_csv(summ, file.path(O, "audit_feature_matrices_summary.csv"))   # saved before the cross step

# ---- 6. cross-matrix: the same transcript must carry the same value everywhere ------------------
pairs <- if (length(mats) >= 2) combn(names(mats), 2, simplify = FALSE) else list()
cross <- map_dfr(pairs, function(p) {
  a <- mats[[p[1]]]; b <- mats[[p[2]]]; sh <- intersect(a$transcript_id_clean, b$transcript_id_clean)
  xa <- a[match(sh, a$transcript_id_clean), ]; xb <- b[match(sh, b$transcript_id_clean), ]
  map_dfr(intersect(intr, intersect(names(a), names(b))), function(f) {
    u <- as.numeric(unlist(xa[[f]])); v <- as.numeric(unlist(xb[[f]])); k <- !is.na(u) & !is.na(v)
    tibble(pair = paste(p, collapse = " vs "), feature = f, shared = length(sh),
           value_mismatch = sum(k & abs(u - v) > TOL), na_pattern_mismatch = sum(is.na(u) != is.na(v)))
  })
})
cross <- if (nrow(cross)) cross %>% filter(value_mismatch > 0 | na_pattern_mismatch > 0) else
  tibble(pair = character(), feature = character(), shared = integer(), value_mismatch = integer(), na_pattern_mismatch = integer())

write_csv(summ, file.path(O, "audit_feature_matrices_summary.csv"))
write_csv(cross, file.path(O, "audit_cross_matrix.csv"))
fail <- summ %>% filter(value_mismatch > 0 | na_only_in_matrix > 0 | na_only_in_recompute > 0)
cat(sprintf("\nAUDIT: %d matrix x feature checks | %d clean | %d with a discrepancy | cross-matrix discrepancies: %d\n",
            nrow(summ), nrow(summ) - nrow(fail), nrow(fail), nrow(cross)))
if (nrow(fail)) print(as.data.frame(fail %>% dplyr::select(-examples)), row.names = FALSE)
if (nrow(cross)) print(as.data.frame(cross), row.names = FALSE)
