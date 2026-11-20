# 01r_positional_gc_matrix.R - GC1 and GC2 on the same internal codons as gc3_internal.
# Writes a SEPARATE matrix + annotation so feature_matrix_tco.rds and feature_annotation.csv
# stay untouched (Kate, 2026-11-18):
#   output/translational_control_overview/feature_matrix_tco_posgc.rds
#   output/translational_control_overview/feature_annotation_tco_posgc.csv
# Use with 02_rf_model.Rmd via matrix_rds / annotation_csv, and drop the aggregate GC columns
# with exclude_features = "cds_gc,transcript_gc" - otherwise the forest can fall back on the
# aggregate and the per-position split is never tested (the 01v_ rationale).
#
# Definition: "internal" = initiator ATG and terminal stop removed (01_'s trim_internal), and
# the column is NA wherever gc3_internal is NA (out-of-frame CDS, missing terminal stop) - so
# the three positions are always computed on exactly the same codons. The script proves that
# by recomputing gc3 with the same function and asserting it reproduces gc3_internal exactly.
# Other cell lines: pass "mcf7six1" to build on that line's own-isoform matrix and CDS
# sequences (-> feature_matrix_tco_mcf7six1_posgc.rds + feature_annotation_mcf7six1_posgc.csv).
#   Rscript code/translational_control_overview/01r_positional_gc_matrix.R [mda231|mcf7six1]
suppressPackageStartupMessages({ library(tidyverse); library(here) })

cell <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "mda231" }
cfg <- switch(cell,
  mda231   = list(mat = "feature_matrix_tco.rds", ann = "feature_annotation.csv",
                  seq = "precomputed_kozak_sequences.rds",
                  out_mat = "feature_matrix_tco_posgc.rds", out_ann = "feature_annotation_tco_posgc.csv"),
  mcf7six1 = list(mat = "feature_matrix_tco_mcf7six1.rds", ann = "feature_annotation_mcf7six1.csv",
                  seq = "precomputed_mcf7six1_kozak_sequences.rds",
                  out_mat = "feature_matrix_tco_mcf7six1_posgc.rds",
                  out_ann = "feature_annotation_mcf7six1_posgc.csv"),
  stop("cell must be mda231 or mcf7six1"))
cat("cell line:", cell, "\n")
O    <- here("output", "translational_control_overview")
feat <- readRDS(file.path(O, cfg$mat))
ann  <- read_csv(file.path(O, cfg$ann), show_col_types = FALSE)
seqs <- readRDS(here("output", "predictive_modeling", cfg$seq))
stopifnot("annotation does not match the base matrix" = setequal(ann$feature, colnames(feat)),
          "positional GC columns already present" =
            !any(c("gc1_internal", "gc2_internal") %in% colnames(feat)))

cds <- seqs$cds_seq %>%
  filter(transcript_id_clean %in% feat$transcript_id_clean) %>%
  transmute(transcript_id_clean, sq = str_to_upper(cds_sequence))
stopifnot("duplicate CDS rows" = !anyDuplicated(cds$transcript_id_clean))

stop_codons <- c("TAA", "TAG", "TGA")
trim_internal <- function(x) {                     # verbatim from 01_
  n <- nchar(x) %/% 3L
  a <- ifelse(substr(x, 1, 3) == "ATG", 4L, 1L)
  b <- ifelse(substr(x, (n - 1) * 3L + 1L, n * 3L) %in% stop_codons, (n - 1) * 3L, n * 3L)
  substr(x, a, b)
}
gc_pos <- function(x, pos) vapply(x, function(z) {
  n <- nchar(z) %/% 3L; if (n < 1L) return(NA_real_)
  k <- seq.int(pos, n * 3L, 3L)
  sum(substr(rep.int(z, n), k, k) %in% c("G", "C")) * 100 / n
}, numeric(1L), USE.NAMES = FALSE)

cds <- cds %>%
  filter(nchar(sq) %% 3L == 0L) %>%                # 01_'s codon_tbl: in-frame only
  mutate(s_int = trim_internal(sq))
stopifnot("internal sequences lost frame" = all(nchar(cds$s_int) %% 3L == 0L))
pos <- cds %>% transmute(transcript_id_clean,
                         gc1_internal = gc_pos(s_int, 1L),
                         gc2_internal = gc_pos(s_int, 2L),
                         gc3_check    = gc_pos(s_int, 3L))

out <- feat %>% left_join(pos, by = "transcript_id_clean")
stopifnot("join changed the row count" = nrow(out) == nrow(feat))
# same codons as gc3_internal: identical NA pattern after masking, and gc3 reproduces exactly
mask <- is.na(out$gc3_internal)
out$gc1_internal[mask] <- NA_real_
out$gc2_internal[mask] <- NA_real_
dev <- max(abs(out$gc3_check - out$gc3_internal), na.rm = TRUE)
stopifnot("recomputed GC3 does not reproduce gc3_internal - not the same codons" =
            dev < 1e-9 && identical(is.na(out$gc1_internal), mask))
cat(sprintf("GC3 reproduction: max deviation %.1e over %d transcripts\n",
            dev, sum(!mask)))
out$gc3_check <- NULL

# --- computational checks ---
stopifnot(
  "gc1_internal outside [0, 100]" = all(out$gc1_internal >= 0 & out$gc1_internal <= 100, na.rm = TRUE),
  "gc2_internal outside [0, 100]" = all(out$gc2_internal >= 0 & out$gc2_internal <= 100, na.rm = TRUE),
  "no target leakage" = !any(c("eif3d_target") %in% colnames(out))
)

# --- biological spot-checks ---
# Human CDS: GC1 ~ 55-60%, GC2 ~ 40% (second position is constrained by amino-acid
# hydrophobicity), GC3 broadest and highest-variance. GC1 > GC2 is the invariant.
med <- sapply(out[c("gc1_internal", "gc2_internal", "gc3_internal")], median, na.rm = TRUE)
sdv <- sapply(out[c("gc1_internal", "gc2_internal", "gc3_internal")], sd, na.rm = TRUE)
print(round(rbind(median = med, sd = sdv), 2))
stopifnot("GC1 median not above GC2 median" = med[["gc1_internal"]] > med[["gc2_internal"]],
          "GC1 median outside [48, 65]%" = between(med[["gc1_internal"]], 48, 65),
          "GC2 median outside [33, 48]%" = between(med[["gc2_internal"]], 33, 48),
          "GC3 is not the most variable position" = sdv[["gc3_internal"]] > max(sdv[1:2]))
# mean of the three positions on internal codons should track cds_gc closely
r_mean <- cor(rowMeans(out[c("gc1_internal", "gc2_internal", "gc3_internal")]), out$cds_gc,
              use = "complete.obs")
cat(sprintf("cor(mean of GC1-3, cds_gc) = %.4f\n", r_mean))
stopifnot("mean of GC1-3 does not track cds_gc" = r_mean > 0.98)
actb <- out %>% filter(symbol == "ACTB")
stopifnot("ACTB missing" = nrow(actb) >= 1,
          "ACTB GC1 outside [50, 75]%" = all(between(actb$gc1_internal, 50, 75)))
cat(sprintf("ACTB GC1/GC2/GC3: %.1f / %.1f / %.1f\n",
            actb$gc1_internal[1], actb$gc2_internal[1], actb$gc3_internal[1]))

ann_new <- bind_rows(ann, tibble(
  feature = c("gc1_internal", "gc2_internal"), block = "intrinsic", family = "gc",
  variant = "single", source = "01r_positional_gc_matrix.R",
  definition = c("percent (0-100); codon position 1, internal codons: initiator and stop removed",
                 "percent (0-100); codon position 2, internal codons: initiator and stop removed"),
  deprecated = FALSE, added_in = "tco"))
stopifnot(setequal(ann_new$feature, colnames(out)))

saveRDS(out, file.path(O, cfg$out_mat))
write_csv(ann_new, file.path(O, cfg$out_ann))
cat("wrote", cfg$out_mat, "(", nrow(out), "x", ncol(out), ") and", cfg$out_ann, "\n")
