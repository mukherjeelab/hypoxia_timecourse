# Shared machinery for the CDS reporter design notebooks.
#
# Sourced by 01_ (MDA-MB-231) and 02_ (MCF7-SIX1). Everything here is
# cell-line agnostic: the reporter sequences, the CDS metrics, and the density
# plotter. Anything cell-line specific (which TE table, which feature matrix,
# which direction check) stays in the notebook.
#
# Expects tidyverse, Biostrings, here and GenomicFeatures to be attached by the
# caller.

STOP_CODONS <- c("TAA", "TAG", "TGA")

# The mouse ODC PEST degron carried by all three lab reporters.
PEST_AA <- "SHGFPPEVEEQAAGTLPMSCAQESGMDRHPAACASARINV"

# The C-terminal extension is 41 aa, not 40: one asparagine sits between the
# reporter ORF and the degron, encoded by the AAT inside the GAATTC EcoRI site
# at the cloning junction (...CTGGCG|AATTCTCACGGC...). Splitting a construct on
# the 40-mer degron alone puts the boundary 3 nt late.
TAIL_AA <- paste0("N", PEST_AA)

codon_split <- function(s) {
  n <- nchar(s) %/% 3L
  substring(s,
            seq.int(1L, by = 3L, length.out = n),
            seq.int(3L, by = 3L, length.out = n))
}

compute_gc3_pct <- function(s) {
  s <- toupper(gsub("\\s", "", s))
  if (is.na(s) || nchar(s) < 3L) return(NA_real_)
  n_codons <- nchar(s) %/% 3L
  starts   <- seq.int(3L, n_codons * 3L, 3L)
  thirds   <- substr(rep.int(s, n_codons), starts, starts)
  sum(thirds %in% c("G", "C")) * 100 / n_codons
}

gc_pct_of <- function(s) {
  as.numeric(Biostrings::letterFrequency(Biostrings::DNAString(s),
                                         letters = "GC", as.prob = TRUE)) * 100
}

# Validate a reporter ORF: in frame, ATG start, no internal stop, and already
# stop-free so it matches the stop-stripped gene CDS sequences.
check_orf <- function(s, name) {
  s   <- toupper(gsub("\\s", "", s))
  cdn <- codon_split(s)
  stopifnot(
    "Reporter ORF length not divisible by 3"     = nchar(s) %% 3 == 0,
    "Reporter ORF does not start with ATG"       = cdn[1] == "ATG",
    "Reporter ORF has an internal stop codon"    =
      !any(cdn[-c(1L, length(cdn))] %in% STOP_CODONS),
    "Reporter ORF still carries a terminal stop" =
      !cdn[length(cdn)] %in% STOP_CODONS
  )
  cat(sprintf("  %-22s OK: %d nt, %d codons\n", name, nchar(s), length(cdn)))
  invisible(TRUE)
}

fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)


# ---------------------------------------------------------------------------
# NanoLuc construct: import and verify
# ---------------------------------------------------------------------------

# Read the construct from the Benchling export and assert it is identical to the
# string notebook 62 hard-codes, by parsing notebook 62's own source. Importing
# and cross-checking rather than retyping is the project rule; parsing the other
# notebook means a change to either file fires this check instead of the two
# drifting apart silently.
load_nanoluc_construct <- function() {
  nl_xlsx <- here::here("accessories", "cds seq.xlsx")
  stopifnot("accessories/cds seq.xlsx not found" = file.exists(nl_xlsx))

  nl_tbl <- readxl::read_excel(nl_xlsx)
  stopifnot("Expected a 'cds' column in cds seq.xlsx" = "cds" %in% names(nl_tbl))
  names(nl_tbl)[1] <- "construct"

  nl_row <- nl_tbl %>% dplyr::filter(stringr::str_detect(tolower(construct), "nanoluc"))
  stopifnot("No nanoluc row found in cds seq.xlsx" = nrow(nl_row) == 1)

  nl_raw <- toupper(gsub("\\s", "", nl_row$cds[1]))
  cat("Imported construct:", nl_row$construct[1], "|", nchar(nl_raw), "nt\n")

  stopifnot(
    "NanoLuc CDS length not divisible by 3" = nchar(nl_raw) %% 3 == 0,
    "NanoLuc CDS has non-ACGT characters"   = !grepl("[^ACGT]", nl_raw),
    "NanoLuc CDS does not start with ATG"   = substr(nl_raw, 1, 3) == "ATG"
  )

  nb62_path <- here::here("code", "predictive_modeling",
                          "62_cds_length_gc_hypoxia_promotes.Rmd")
  stopifnot("Notebook 62 source not found" = file.exists(nb62_path))
  nb62_src <- paste(readLines(nb62_path, warn = FALSE), collapse = "\n")
  nb62_nl  <- stringr::str_match(nb62_src,
                                 'nl_seq <- toupper\\("([ACGTacgt]+)"\\)')[, 2]
  stopifnot("Could not parse nl_seq out of notebook 62" = !is.na(nb62_nl))
  stopifnot("Imported NanoLuc CDS does not match notebook 62's nl_seq" =
              identical(nl_raw, toupper(nb62_nl)))
  cat("Imported sequence is identical to notebook 62 nl_seq (", nchar(nl_raw),
      "nt )\n")

  nl_codons <- codon_split(nl_raw)
  has_stop  <- nl_codons[length(nl_codons)] %in% STOP_CODONS
  nl_full   <- if (has_stop) substr(nl_raw, 1, nchar(nl_raw) - 3) else nl_raw
  cat("Terminal stop codon present:", has_stop,
      "| stop-free construct:", nchar(nl_full), "nt\n")

  nl_aa <- as.character(Biostrings::translate(Biostrings::DNAString(nl_full)))
  stopifnot(
    "Internal stop codon in NanoLuc construct" = !grepl("\\*", nl_aa),
    "NanoLuc construct is missing its C-terminal PEST degron" =
      endsWith(nl_aa, PEST_AA),
    "Construct does not end with the 41-aa junction + degron extension" =
      endsWith(nl_aa, TAIL_AA)
  )

  nl_core <- substr(nl_full, 1, nchar(nl_full) - 3 * nchar(TAIL_AA))
  nl_pest <- substr(nl_full, nchar(nl_core) + 1, nchar(nl_full))
  stopifnot(
    "Tag split did not partition the construct" =
      paste0(nl_core, nl_pest) == nl_full,
    "NanoLuc catalytic domain is not 171 codons" = nchar(nl_core) == 513,
    "NanoLuc does not start with the expected MVFTLEDFVGDWRQTAGYNL" =
      startsWith(as.character(
        Biostrings::translate(Biostrings::DNAString(nl_core))),
        "MVFTLEDFVGDWRQTAGYNL")
  )

  cat("Reporter ORF validation:\n")
  check_orf(nl_full, "construct")
  check_orf(nl_core, "NanoLuc only")

  list(full = nl_full, core = nl_core, pest = nl_pest, aa = nl_aa)
}


# Load the codon-adapted / anti-adapted / GC3-matched control NanoLuc sequences
# designed in 64_nanoluc_codon_adaptation.Rmd, and attach the unmodified
# junction + PEST tail each was cloned with. Notebook 64 exports 513 nt ORFs,
# but the tail is part of the CDS and was never recoded, so the sequence that
# actually exists is ORF + tail.
load_nanoluc_variants <- function(nl) {
  des_path <- here::here("output", "predictive_modeling",
                         "nanoluc_codon_designed_sequences.csv")
  stopifnot("Notebook 64 designed-sequence CSV not found; knit 64 first" =
              file.exists(des_path))

  des <- readr::read_csv(des_path, show_col_types = FALSE)
  stopifnot("Unexpected columns in designed-sequence CSV" =
              all(c("name", "bases", "length_nt", "gc3_pct") %in% names(des)))

  var_order <- c(nanoluc_original       = "NanoLuc parent",
                 nanoluc_3d_adapted     = "3d-adapted",
                 nanoluc_3d_antiadapted = "3d-anti-adapted",
                 nanoluc_3d_negctrl     = "3d-negative control")
  stopifnot("Designed-sequence CSV does not hold the four expected variants" =
              setequal(des$name, names(var_order)))

  des <- des %>%
    dplyr::mutate(bases = toupper(gsub("\\s", "", bases)),
                  label = unname(var_order[name])) %>%
    dplyr::arrange(match(name, names(var_order)))

  stopifnot("CSV nanoluc_original is not identical to the construct's ORF" =
              identical(des$bases[des$name == "nanoluc_original"], nl$core))
  cat("CSV parent sequence is identical to the imported construct's 513 nt ORF\n")

  cat("Variant ORF validation:\n")
  purrr::walk2(des$bases, des$label, check_orf)

  # Biological spot-check: every variant is a synonymous recoding, so each must
  # encode a protein identical to the parent. A design error that changed an
  # amino acid would still pass every range and GC check.
  ref_aa <- as.character(Biostrings::translate(Biostrings::DNAString(nl$core)))
  var_aa <- purrr::map_chr(des$bases, ~ as.character(
    Biostrings::translate(Biostrings::DNAString(.x))))
  stopifnot(
    "A designed variant does not encode the parent NanoLuc protein" =
      all(var_aa == ref_aa),
    "Designed variants are not all the same length" =
      dplyr::n_distinct(nchar(des$bases)) == 1
  )
  cat("All", nrow(des), "variants encode the identical 171 aa NanoLuc protein\n")

  des <- des %>% dplyr::mutate(bases_full = paste0(bases, nl$pest))

  cat("Variant construct validation (ORF + unmodified tail):\n")
  purrr::walk2(des$bases_full, des$label, check_orf)

  stopifnot(
    "A variant construct is not the parent construct's length" =
      all(nchar(des$bases_full) == nchar(nl$full)),
    "A variant construct does not end with the PEST degron" =
      all(purrr::map_lgl(des$bases_full, ~ endsWith(as.character(
        Biostrings::translate(Biostrings::DNAString(.x))), PEST_AA))),
    "The parent variant's construct is not identical to the imported construct" =
      des$bases_full[des$name == "nanoluc_original"] == nl$full
  )
  cat("All variant constructs are", nchar(nl$full),
      "nt and carry the same unmodified tail\n")

  var_df <- des %>%
    dplyr::transmute(
      reporter   = label,
      basis      = "full construct (ORF + unmodified junction/PEST tail)",
      # gc3_orf_csv must be bound BEFORE anything recomputes GC3: transmute
      # evaluates in order, so reading gc3_pct afterwards would return the
      # recomputed value and the assertion below could never fail.
      gc3_orf_csv = gc3_pct,
      orf_length  = nchar(bases),
      gc_orf      = purrr::map_dbl(bases, gc_pct_of),
      gc3_orf     = purrr::map_dbl(bases, compute_gc3_pct),
      cds_length  = nchar(bases_full),
      gc_pct      = purrr::map_dbl(bases_full, gc_pct_of),
      gc3_pct     = purrr::map_dbl(bases_full, compute_gc3_pct),
      n_codons_changed = n_codons_changed
    )

  stopifnot(
    "Recomputed ORF GC3 does not reproduce notebook 64's value" =
      max(abs(var_df$gc3_orf - var_df$gc3_orf_csv)) < 1e-8,
    "Variant GC% out of range"  =
      all(var_df$gc_pct  >= 0 & var_df$gc_pct  <= 100),
    "Variant GC3% out of range" =
      all(var_df$gc3_pct >= 0 & var_df$gc3_pct <= 100),
    "Parent construct GC3 disagrees with the imported construct" =
      abs(var_df$gc3_pct[var_df$reporter == "NanoLuc parent"] -
            compute_gc3_pct(nl$full)) < 1e-8
  )
  cat("Recomputed ORF GC3 reproduces notebook 64's values exactly\n")

  var_df
}


# ---------------------------------------------------------------------------
# CDS metrics from the genome
# ---------------------------------------------------------------------------

# One row per transcript with all three metrics, so every panel is guaranteed to
# describe the same transcript set and the printed n cannot drift between
# panels. GENCODE txdb CDS ranges INCLUDE the terminal stop codon (ACTB CDS =
# 1128 nt = 375 aa + stop), so it is stripped to share the stop-free basis used
# for the reporters.
# Internal codons: the basis every composition metric in this folder uses.
# Drop the initiator, whatever it is, and the terminal stop, whichever it is.
#
# *** THIS DIVERGES FROM THE MODELLING, DELIBERATELY. ***
#
# code/translational_control_overview/01_build_feature_matrix.Rmd defines
# `gc3_internal` with a CONDITIONAL start trim - `ifelse(substr(x,1,3)=="ATG",
# 4L, 1L)` - so a transcript that initiates on a near-cognate codon keeps its
# first codon. That guard exists to avoid over-trimming anything non-canonical.
#
# The trim here is unconditional at position 1, because the rationale for
# removing the initiator does not depend on its identity: a CUG that genuinely
# initiates a transcript is as forced, and as unavailable to synonymous choice,
# as an AUG. Leaving it in gives those transcripts one extra codon in a
# denominator that is supposed to count only choosable positions. It affects 25
# of the 10,054 expressed transcripts - CTG 14, GTG 6, ATT 2, ACG 1, ATA 1,
# TTG 1 - and moves no reported median.
#
# The STOP trim stays conditional, and that asymmetry is load-bearing: the
# reporter constructs are built stop-free by design, so an unconditional trim
# would remove a real codon from every one of them. On the gene side it makes no
# difference, since build_cds_seqs has already dropped anything without a stop.
#
# Composition only. CDS LENGTH and total GC keep the stop-stripped full CDS.
trim_internal <- function(x) {
  stopifnot(
    "trim_internal needs in-frame sequences" = all(nchar(x) %% 3 == 0),
    # Position 1 is removed unconditionally, so anything passed here must be a
    # whole ORF. Two codons is the minimum that leaves something behind.
    "trim_internal needs at least two codons" = all(nchar(x) >= 6)
  )
  n <- nchar(x) %/% 3L
  b <- ifelse(substr(x, (n - 1) * 3L + 1L, n * 3L) %in% STOP_CODONS,
              (n - 1) * 3L, n * 3L)
  out <- substr(x, 4L, b)
  # A stop-bearing CDS needs 3 codons to leave one behind, a stop-free
  # construct only 2. Rather than encode that split, assert the result.
  stopifnot("trim_internal left nothing: sequence too short to have internal codons" =
              all(nchar(out) > 0))
  out
}


# Composition of a whole ORF, on the same bases the gene distributions use. Use
# these for anything that becomes a reference line against those distributions;
# part-composition tables describe fragments as ordered and keep the plain
# versions.
#
# The two metrics sit on DIFFERENT bases, and that is deliberate rather than an
# oversight: it is the modelling's pairing. GC3 is internal-codon, matching
# `gc3_internal`, because it exists to measure synonymous choice and the
# initiator is not a choice. Total GC is the full stop-stripped CDS, matching
# `cds_gc`, which the modelling has no internal variant of - so this folder's GC
# column stays comparable with every `cds_gc` elsewhere in the project. The
# difference between the two bases is one codon, under 0.2 points on a median
# CDS.
orf_gc_pct  <- function(s) gc_pct_of(s)
orf_gc3_pct <- function(s) compute_gc3_pct(trim_internal(s))


# Complete CDS: in frame, and carrying a terminal stop. Transcripts without one
# are GENCODE cds_end_NF - the CDS starts at ATG but the annotation runs out -
# and the modelling masks every codon-derived feature for them rather than
# trimming them, because the sequence is incomplete and compositionally skewed
# (GC3 median 48.9 against 57.2). They are dropped here for the same reason.
# The terminal stop is still present in what this returns; trim_internal takes
# it off.
build_cds_seqs <- function(txdb_path) {
  txdb <- AnnotationDbi::loadDb(txdb_path)
  cds  <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)

  seqs <- GenomicFeatures::extractTranscriptSeqs(
    BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, cds)
  names(seqs) <- sub("\\..*", "", names(seqs))
  stopifnot("Duplicate transcript IDs after version stripping" =
              !any(duplicated(names(seqs))))

  # At least 3 codons, because the initiator and the terminal stop both come off
  # and a composition metric needs something left. A 1- or 2-codon CDS has no
  # internal codons at all, so it cannot contribute to GC3 under any definition.
  seqs <- seqs[BiocGenerics::width(seqs) >= 9 &
                 BiocGenerics::width(seqs) %% 3 == 0]

  last_codon <- as.character(
    Biostrings::subseq(seqs, start = BiocGenerics::width(seqs) - 2L))
  tx_has_stop <- last_codon %in% STOP_CODONS
  cat("Transcripts in frame:", length(seqs),
      "| with a terminal stop:", sum(tx_has_stop),
      paste0("(", round(100 * mean(tx_has_stop), 1), "%)"),
      "| dropped as cds_end_NF:", sum(!tx_has_stop), "\n")

  seqs[tx_has_stop]
}


build_cds_metrics <- function(txdb_path) {
  complete <- build_cds_seqs(txdb_path)

  # Three quantities, two bases, matching the modelling's own pairing.
  # LENGTH and TOTAL GC are the stop-stripped CDS - what this folder has always
  # plotted, and what `cds_gc` is elsewhere in the project. GC3 is internal
  # codons, matching `gc3_internal`, because it exists to measure synonymous
  # choice and the forced initiator is not one. The initiator biases short
  # sequences hardest: one codon is 0.2% of a median CDS but 2.6% of the 117 nt
  # reporter, which is why the basis matters here at all.
  nostop   <- Biostrings::subseq(complete, start = 1L,
                                 end = BiocGenerics::width(complete) - 3L)
  internal <- trim_internal(as.character(complete))

  stopifnot(
    "Internal sequences lost frame"   = all(nchar(internal) %% 3 == 0),
    "An internal sequence is empty"   = all(nchar(internal) > 0),
    # Every sequence here carries a stop, so trimming must always remove it,
    # and must remove the initiator from all but the non-AUG starters.
    "Internal trim did not remove the terminal stop" =
      all(nchar(internal) <= BiocGenerics::width(complete) - 3L)
  )

  out <- tibble::tibble(
    transcript_id_clean = names(complete),
    cds_length = BiocGenerics::width(nostop),
    gc_pct  = as.numeric(Biostrings::letterFrequency(
      nostop, letters = "GC", as.prob = TRUE)) * 100,
    gc3_pct = sapply(internal, compute_gc3_pct, USE.NAMES = FALSE)
  )
  stopifnot("Stop-stripped CDS lengths must remain divisible by 3" =
              all(out$cds_length %% 3 == 0))
  out
}


# ---------------------------------------------------------------------------
# Density panels
# ---------------------------------------------------------------------------

# expressed / promoted: data frames carrying transcript_id_clean,
# gene_id_clean and the three metric columns.
plot_density_luc <- function(expressed, promoted, xvar, xlab, title, outfile,
                             luc_tbl = NULL, group_cols,
                             xscale = 1, xmax = NULL) {
  # luc_tbl = NULL draws the gene distributions alone, with no reference lines.
  # An empty table is built with the columns the code below expects rather than
  # branching around every use of it, so the no-reference case walks the same
  # path as the normal one instead of being a second implementation.
  if (is.null(luc_tbl) || nrow(luc_tbl) == 0) {
    luc_tbl <- tibble::tibble(reporter = character(0),
                              color    = character(0),
                              linetype = character(0))
    luc_tbl[[xvar]] <- numeric(0)
  }
  df <- dplyr::bind_rows(
    expressed %>% dplyr::mutate(group = names(group_cols)[1]),
    promoted  %>% dplyr::mutate(group = names(group_cols)[2])
  ) %>%
    dplyr::mutate(group = factor(group, levels = names(group_cols)),
                  .x    = .data[[xvar]] / xscale)

  # Every legend entry carries the number it stands for: the group median for
  # each density curve, the construct's own value for each reference line. Units
  # follow the metric, not the axis - the length axis is drawn in kb because kb
  # ticks are readable, but a construct is specified in nt, so the label says nt.
  fmt_metric <- function(x) {
    if (xvar == "cds_length") sprintf("%s nt", fmt_n(round(x)))
    else sprintf("%.1f%%", x)
  }

  lev    <- levels(df$group)
  n_tab  <- df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(n = dplyr::n(), med = stats::median(.data[[xvar]]),
                     .groups = "drop")
  labs_v <- setNames(sprintf("%s (n = %s, median %s)",
                             as.character(n_tab$group), fmt_n(n_tab$n),
                             fmt_metric(n_tab$med)),
                     as.character(n_tab$group))
  df <- df %>% dplyr::mutate(
    group_lab = factor(labs_v[as.character(group)],
                       levels = unname(labs_v[lev])))
  cols_lab <- setNames(unname(group_cols[lev]), unname(labs_v[lev]))

  # Disjoint comparison: promoted vs every other expressed transcript. Reported
  # to the console, not annotated on the panel; the per-group n is already
  # carried by the legend labels.
  x_pro  <- promoted[[xvar]] / xscale
  x_rest <- expressed[[xvar]][
    !expressed$gene_id_clean %in% promoted$gene_id_clean] / xscale
  wt <- stats::wilcox.test(x_pro, x_rest)

  luc_plot     <- luc_tbl %>%
    dplyr::mutate(.x = .data[[xvar]] / xscale,
                  reporter_lab = sprintf("%s (%s)", reporter,
                                         fmt_metric(.data[[xvar]])))
  stopifnot("Two reference lines share a legend label" =
              !any(duplicated(luc_plot$reporter_lab)))
  luc_cols_use <- setNames(luc_plot$color, luc_plot$reporter_lab)

  # Colour and linetype are both keyed on the same set of names (the group
  # labels, then the reporters), so the two scales merge into one legend rather
  # than drawing two. Linetype matters because synonymous variants can land on
  # identical GC/GC3, where a solid line would be completely hidden by the one
  # drawn over it.
  luc_lty_use <- setNames(luc_plot$linetype, luc_plot$reporter_lab)
  lty_lab     <- setNames(rep("solid", length(cols_lab)), names(cols_lab))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .x, color = group_lab,
                                        linetype = group_lab)) +
    ggplot2::geom_density(fill = NA, linewidth = 1, trim = TRUE,
                          key_glyph = "path") +
    # With no reference lines the vline layer draws nothing, but it still
    # contributes a vertical stroke to every legend key, turning each group's
    # line glyph into a cross. Suppressing the layer's legend in that case
    # leaves the keys as plain horizontal lines. NA is ggplot's default, so the
    # annotated panels are unaffected.
    ggplot2::geom_vline(
      data = luc_plot,
      ggplot2::aes(xintercept = .x, color = reporter_lab,
                   linetype = reporter_lab),
      linewidth = 0.9,
      show.legend = if (nrow(luc_plot) == 0) FALSE else NA) +
    ggplot2::scale_color_manual(
      values = c(cols_lab, luc_cols_use),
      breaks = c(names(cols_lab), names(luc_cols_use)), name = NULL) +
    ggplot2::scale_linetype_manual(
      values = c(lty_lab, luc_lty_use),
      breaks = c(names(lty_lab), names(luc_lty_use)), name = NULL) +
    ggplot2::labs(x = xlab, y = "Density", title = title) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.title  = ggplot2::element_text(size = 17),
      axis.text   = ggplot2::element_text(size = 15, colour = "black"),
      legend.text = ggplot2::element_text(size = 14),
      plot.title  = ggplot2::element_text(size = 18))

  # Density is bounded below by 0, so the y axis starts there rather than
  # leaving ggplot's default gap under the curves.
  #
  # The x axis is anchored at 0 only where xmax is given, i.e. the length panel,
  # where 0 is a meaningful origin and the long tail needs clipping. GC and GC3
  # autoscale: forcing those to start at 0 would push the whole distribution
  # into the right-hand third and squash the gap between the reporter line and
  # the gene curves, which is the thing to read.
  #
  # coord_cartesian zooms without dropping rows, so the densities are still
  # estimated from the full data and only the visible range is clipped.
  if (!is.null(xmax)) {
    n_off <- sum(df$.x > xmax)
    cat(sprintf("%s: x axis clipped at %g; %s of %s transcripts (%.1f%%) beyond it\n",
                title, xmax, fmt_n(n_off), fmt_n(nrow(df)),
                100 * n_off / nrow(df)))
  }
  p <- p +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::coord_cartesian(
      xlim = if (is.null(xmax)) NULL else c(0, xmax), ylim = c(0, NA))

  # coord_cartesian's xlim still carries ggplot's default 5% expansion, so the
  # drawn axis would start below zero without this.
  if (!is.null(xmax)) {
    p <- p + ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.02)))
  }

  print(p)
  ggplot2::ggsave(here::here("plots", outfile), p, width = 9.5, height = 5.5)
  cat(title, "- Wilcoxon promotes vs other expressed: p =",
      format.pval(wt$p.value, digits = 3), "\n")
  invisible(p)
}


# Percentile of each reference line within each gene distribution. A reference
# line is only interpretable with the percentile attached.
reporter_percentiles <- function(expressed, promoted, ref_all) {
  pct_of <- function(x, v) 100 * mean(x <= v)

  tidyr::expand_grid(
    reporter = ref_all$reporter,
    metric   = c("cds_length", "gc_pct", "gc3_pct")
  ) %>%
    dplyr::left_join(
      ref_all %>%
        dplyr::select(reporter, cds_length, gc_pct, gc3_pct) %>%
        tidyr::pivot_longer(-reporter, names_to = "metric",
                            values_to = "reporter_value"),
      by = c("reporter", "metric")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pct_expressed    = pct_of(expressed[[metric]], reporter_value),
      pct_promoted     = pct_of(promoted[[metric]],  reporter_value),
      median_expressed = median(expressed[[metric]]),
      median_promoted  = median(promoted[[metric]])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 1)))
}
