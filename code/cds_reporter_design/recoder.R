# ---------------------------------------------------------------------------
# The synonymous recoder, as a shared implementation.
#
# This is a faithful port of section 5 of
# `08_gc_matched_recode_egfp_ladder.Rmd`, with the quantities that were globals
# there - the GC3 target, the transcriptome's within-class third-base splits,
# the codon menu and usage weights - passed in instead. `08_` is the design of
# record and is deliberately NOT changed to source this file: it keeps its
# inline copy, and `12_` gates itself on reproducing `08_`'s exported sequences
# bit-for-bit. That gate is what makes two copies safe. If this port ever drifts
# from the original, the gate fires rather than a wrong construct being drawn.
#
# The one functional difference is that `target_gc3` is an argument rather than
# the expressed-transcriptome median. Everything else - the two-stage third-base
# assignment, the GC3-before-balance priority, the clamping in `block_needs`,
# the motif redraw loop - is unchanged.
# ---------------------------------------------------------------------------

BASES_R <- c("A", "C", "G", "T")

# Count third bases in a codon vector.
base3_counts <- function(cv)
  table(factor(substr(cv, 3, 3), levels = BASES_R)) |>
    as.integer() |> setNames(BASES_R)

# Split one class of positions between its two bases, putting `n_first` of them
# on the first base. Positions that can only take one of the two are settled
# first; the rest are filled at random, which is where copy-to-copy divergence
# comes from.
split_class <- function(av, pair, n_first) {
  can1 <- vapply(av, function(x) pair[1] %in% x, logical(1))
  can2 <- vapply(av, function(x) pair[2] %in% x, logical(1))
  res  <- ifelse(can1 & !can2, pair[1], pair[2])
  free <- which(can1 & can2)
  n_more <- min(max(n_first - sum(res == pair[1]), 0), length(free))
  if (n_more > 0)
    res[if (n_more == length(free)) free else sample(free, n_more)] <- pair[1]
  res
}

# Give each position a third base, in two stages that encode the priority
# `block_needs` states: GC3 first, the within-class balance second.
#
#   Stage 1  decide which positions end G/C at all, hitting the GC3 count.
#   Stage 2  split those between G and C, and the rest between A and T.
#
# `avail` is the third bases allowed AT EACH POSITION, not per amino acid, so a
# block can narrow the menu at particular positions.
assign_third_bases <- function(avail, need, max_try = 50) {
  gc_ok <- vapply(avail, function(x) any(x %in% c("G", "C")), logical(1))
  at_ok <- vapply(avail, function(x) any(x %in% c("A", "T")), logical(1))
  n_gc  <- need[["G"]] + need[["C"]]

  forced_gc <- gc_ok & !at_ok
  free      <- which(gc_ok & at_ok)
  n_more    <- min(max(n_gc - sum(forced_gc), 0), length(free))

  best <- NULL; best_cost <- Inf
  for (i in seq_len(max_try)) {
    is_gc <- forced_gc
    if (n_more > 0)
      is_gc[if (n_more == length(free)) free else sample(free, n_more)] <- TRUE

    out <- character(length(avail))
    out[is_gc]  <- split_class(avail[is_gc],  c("G", "C"), need[["G"]])
    out[!is_gc] <- split_class(avail[!is_gc], c("A", "T"), need[["A"]])

    cnt <- base3_counts(paste0("..", out))
    if (all(cnt == need)) return(out)
    # Rank near-misses by GC3 error first, then by the within-class error.
    cost <- 1000 * abs(sum(cnt[c("G", "C")]) - n_gc) + sum(abs(cnt - need))
    if (cost < best_cost) { best_cost <- cost; best <- out }
  }
  best
}

# Build the recoder for one GC3 target.
#
# target_gc3 : construct-wide GC3 percentage to aim at
# R_G, R_A   : the expressed transcriptome's within-class splits - the share of
#              GC-enders that are G, and of AT-enders that are A. These are held
#              at their transcriptome values across every target, so a variant
#              differs from its siblings in how many third bases are G/C, not in
#              how the G/C and A/T classes are divided internally.
# MENU       : per amino acid, usable codons grouped by third base
# USE        : named vector of synonymous usage proportions, for weighting
# AA         : codon-to-amino-acid lookup
# max_homopolymer : longest run of one base allowed in the assembled construct
make_recoder <- function(target_gc3, R_G, R_A, MENU, USE, AA,
                         max_homopolymer) {
  stopifnot(
    "target_gc3 must be a single percentage in (0, 100)" =
      length(target_gc3) == 1 && target_gc3 > 0 && target_gc3 < 100,
    "R_G and R_A must be proportions" =
      R_G >= 0 && R_G <= 1 && R_A >= 0 && R_A <= 1,
    "USE must be a named numeric vector of codon proportions" =
      is.numeric(USE) && !is.null(names(USE)),
    "AA lookup lost its names" = !is.null(names(AA))
  )

  # How many codons of the WHOLE construct should end in each base, given the
  # GC3 target and the transcriptome's within-class splits.
  target_counts <- function(n_total) {
    n_gc <- round(target_gc3 / 100 * n_total)
    n_at <- n_total - n_gc
    out  <- c(A = round(n_at * R_A), C = 0, G = round(n_gc * R_G), T = 0)
    out["C"] <- n_gc - out[["G"]]
    out["T"] <- n_at - out[["A"]]
    out[BASES_R]
  }

  # What this block must contribute for the whole construct to land on target.
  #
  # GC3 is the primary target and the within-class balance is secondary, so they
  # are resolved in that order. It matters whenever the fixed prefix already
  # overshoots or undershoots - as it does here, where HiBiT and P2A are 30 of
  # the bare rung's 54 codons at GC3 74%. Settling the GC total first spends
  # every free codon on pulling GC3 toward the target; splitting the four counts
  # at once would instead spend some of them fixing A against T.
  #
  # The clamps are what make an unreachable target degrade rather than error: a
  # block that cannot supply the whole shortfall supplies as much as it can, and
  # the next block's `need` is recomputed against the prefix actually built.
  block_needs <- function(n_block, prefix_cv) {
    clamp <- function(x, hi) min(max(x, 0), hi)
    # Internal basis: the construct's first codon is its fixed initiator ATG and
    # is no part of a target about synonymous choice, so it is excluded from
    # both the prefix counts and the total the target is computed over.
    pre <- base3_counts(prefix_cv[-1])
    tot <- target_counts(length(prefix_cv) - 1 + n_block)

    n_gc <- clamp(tot[["G"]] + tot[["C"]] - pre[["G"]] - pre[["C"]], n_block)
    n_at <- n_block - n_gc
    g    <- clamp(tot[["G"]] - pre[["G"]], n_gc)
    a    <- clamp(tot[["A"]] - pre[["A"]], n_at)

    c(A = a, C = n_gc - g, G = g, T = n_at - a)[BASES_R]
  }

  # cv    : the codons to recode (the protein is read off these)
  # need  : how many of them must end in each of A, C, G, T
  # menus : one entry per position, each the codon options grouped by third
  #         base. Defaults to the construct-wide MENU for that position's amino
  #         acid; a block passes its own to drop a codon at some positions and
  #         not others, which a menu keyed by amino acid cannot express.
  recode_codons <- function(cv, need, menus = NULL) {
    if (is.null(menus)) menus <- MENU[AA[cv]]
    stopifnot("menus must have one entry per codon" = length(menus) == length(cv))
    b3 <- assign_third_bases(lapply(menus, names), need)
    # The third base is now fixed, so the remaining choice is between codons
    # that share it and differ only in positions 1-2. Weighted by 231 usage.
    mapply(function(m, b) {
      opts <- m[[b]]
      if (length(opts) == 1) opts else sample(opts, 1, prob = USE[opts])
    }, menus, b3, USE.NAMES = FALSE)
  }

  # Two motif constraints have to hold: the Waldo linker's GAATTC must stay the
  # ONLY EcoRI site in the ORF, and long homopolymers are a synthesis problem.
  # Rather than repair a sequence in place, just draw again.
  # stringr rather than base grepl: `08_` tests the homopolymer bound with
  # str_detect, and TRE and ICU do not have to agree on a backreference. The
  # bit-identity gate in `12_` would catch a disagreement, but matching the
  # original's regex engine means there is nothing to catch.
  #
  # `n_linkers = NA` means the construct has no EcoRI requirement at all, and
  # the site count is simply not constrained. The eGFP ladder needs the count
  # pinned because its Waldo linker carries the only intended site and a second
  # one would break the cloning; a construct built without that linker has no
  # such site to protect, and demanding zero would be an invented constraint
  # that silently rejects good draws. The homopolymer bound always applies.
  is_clean <- function(seq_str, n_linkers) {
    ecori_ok <- is.na(n_linkers) ||
      stringr::str_count(seq_str, stringr::fixed("GAATTC")) == n_linkers
    ecori_ok &&
      !stringr::str_detect(
        seq_str, paste0("(A|C|G|T)\\1{", max_homopolymer, ",}"))
  }

  # Redraw the block until the assembled construct passes both motif checks, and
  # any block-specific constraint `extra_ok` adds.
  #
  # `assemble` places the drawn block into the construct instead of a prefix to
  # paste it after, because the GS linker sits UPSTREAM of the fixed P2A rather
  # than at the end of what is already built. Checking motifs on the wrong
  # assembly order would test junctions the construct does not have.
  draw_clean_block <- function(cv, need, assemble, n_linkers,
                               extra_ok = function(block) TRUE, menus = NULL,
                               max_try = 200) {
    for (i in seq_len(max_try)) {
      block <- recode_codons(cv, need, menus = menus)
      if (!extra_ok(block)) next
      full <- paste0(assemble(block), collapse = "")
      if (is_clean(full, n_linkers)) return(list(block = block, tries = i))
    }
    stop("No clean draw after ", max_try, " attempts at GC3 target ", target_gc3)
  }

  list(target_gc3       = target_gc3,
       target_counts    = target_counts,
       block_needs      = block_needs,
       recode_codons    = recode_codons,
       is_clean         = is_clean,
       draw_clean_block = draw_clean_block)
}
