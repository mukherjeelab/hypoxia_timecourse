# 34a_synonymous_variants.R - generate synonymous recodings of a CDS across a GC3 gradient.
#
# For the question "what GC3 do the models actually prefer?", the sweep has to move along the
# axis a codon design can actually move along: synonymous substitution. The protein is fixed, so
# every variant encodes the same peptide and differs only in wobble choice.
#
# The knob is p, the probability of taking a G/C-ending synonymous codon where one exists.
# p = 0 gives the lowest reachable GC3, p = 1 the highest. Amino acids with a single codon (Met,
# Trp) are untouched, so neither extreme reaches 0 or 100.
#
# This deliberately does NOT hold CAI constant. Across the transcriptome GC3 and CAI correlate
# at rho = 0.85, so they cannot be moved independently by synonymous choice; the sweep therefore
# measures the NET effect of recoding, which is what a construct design actually buys.
suppressPackageStartupMessages({ library(tidyverse); library(Biostrings) })

.gc_code <- Biostrings::getGeneticCode("SGC0")
.syn_map <- split(names(.gc_code), unname(.gc_code))     # amino acid (and *) -> its codons

#' One synonymous variant of `cds` at GC-preference `p`.
#' @param cds character, in frame, including the terminal stop
#' @param p numeric in [0,1]
#' @param keep_ends keep the initiator and the terminal stop codon untouched
synonymous_variant <- function(cds, p, keep_ends = TRUE) {
  n <- nchar(cds) %/% 3L
  cod <- substring(cds, seq.int(1L, n * 3L, 3L), seq.int(3L, n * 3L, 3L))
  idx <- seq_along(cod)
  if (keep_ends) idx <- idx[-c(1L, length(cod))]
  for (i in idx) {
    aa  <- .gc_code[[cod[i]]]
    alt <- .syn_map[[aa]]
    if (length(alt) < 2L) next
    gc_end <- alt[substr(alt, 3, 3) %in% c("G", "C")]
    at_end <- alt[substr(alt, 3, 3) %in% c("A", "T")]
    pool <- if (!length(gc_end)) at_end else if (!length(at_end)) gc_end else
              if (runif(1) < p) gc_end else at_end
    cod[i] <- if (length(pool) == 1L) pool else sample(pool, 1L)
  }
  out <- paste0(cod, collapse = "")
  stopifnot("variant changed the protein" =
              identical(as.character(Biostrings::translate(DNAString(out))),
                        as.character(Biostrings::translate(DNAString(cds)))),
            "variant changed length" = nchar(out) == nchar(cds))
  out
}

#' A sweep: `reps` variants at each p, labelled with the realised GC3.
synonymous_sweep <- function(cds, ps = seq(0, 1, by = 0.05), reps = 4L, seed = 9L) {
  set.seed(seed)
  # substring(), not substr(): substr() is not vectorised over start/stop and silently returns
  # a single character, which makes every variant score 0 or 100.
  gc3_of <- function(s) { k <- seq.int(3L, nchar(s), 3L); 100 * mean(substring(s, k, k) %in% c("G", "C")) }
  map_dfr(ps, function(p) map_dfr(seq_len(reps), function(r) {
    v <- synonymous_variant(cds, p)
    tibble(p = p, rep = r, id = sprintf("p%03d_r%d", round(p * 100), r), cds = v, gc3 = gc3_of(v))
  }))
}
