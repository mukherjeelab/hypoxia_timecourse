# 33a_parse_plasmid_sequences.R - parse accessories/transcribed_plasmids.xlsx into the canonical
# per-construct 5'UTR / CDS / 3'UTR table that 33_ turns into model features.
#
# The spreadsheet is Kate's record and is never modified; this writes a derived CSV.
#
# SCOPE: the four synonymous NanoLuc constructs AND the firefly normaliser - five rows. The
# four are a controlled synonymous series and carry assertions firefly cannot satisfy (one CDS
# length, one protein, one 5'UTR, one 3'UTR); those stay scoped to the four. See "firefly".
#
# SPLICING. These transcripts carry a 229 nt intron in the 5'UTR (sequence supplied by Kate,
# 2026-10-06; verified here as an exact single-copy match at 236-464 in all four, obeying GT..AG).
# The models describe mature mRNA, so the intron is removed: the 5'UTR goes 554 -> 325 nt. The
# intron lies wholly inside the 5'UTR, which is asserted, so the CDS and 3'UTR are untouched.
#
# FIREFLY (included 2026-11-18; previously excluded). Two issues, now resolved:
#   1. TSS. Its transcript as entered begins 644 nt upstream of the transcription start site,
#      inside the CMV promoter. The start is located from the CMV landmark "GTTTAGTGAACCGTCAGAT"
#      - the "GTCAGAT" inside it is exactly where the NanoLuc transcripts begin - giving a 433 nt
#      unspliced 5'UTR. Corroborated independently: the 90 nt by which its intron's upstream flank
#      is shorter than the NanoLuc's equals the 90 nt difference in their TSS-to-shared-element
#      linkers. Asserted below, and the 5'UTR must end on the gccacc Kozak.
#   2. Intron. Its copy carries a CAGG deletion at the 3' splice site relative to the NanoLuc
#      backbone, so the annotated intron's last two nt are CA, not the AG a canonical acceptor
#      needs. PER KATE (2026-11-18): trust the Promega annotation for their own vector. The
#      annotated intron is FF_INTRON below (226 nt), giving a 207 nt mature 5'UTR.
#      The deletion itself is ASSERTED rather than excepted: firefly's intron tail + "CAGG" must
#      reproduce the NanoLuc tail plus its first exon base. The alternative reading - splicing to
#      the first downstream AG, a 239 nt intron and a 194 nt 5'UTR - is 0.09 log2 units away and
#      is recorded in code/cds_reporter_design/19_firefly_vs_nanoluc_target_features.Rmd, which
#      derives all of this from sequence. Deciding between them needs RT-PCR, RNA-seq or a
#      sequenced transcript.
#   Firefly's intron length CANNOT be derived from where the two copies diverge: the divergence is
#   internal (AGC -> AAG at nt 152 of 229), not at the 3' end. Hence an explicit constant, located
#   and checked, exactly as INTRON is for the NanoLuc set.
#   Rscript code/translational_control_overview/33a_parse_plasmid_sequences.R
suppressPackageStartupMessages({ library(tidyverse); library(readxl); library(here); library(Biostrings) })

XLSX <- here("accessories", "transcribed_plasmids.xlsx")
OUT  <- here("accessories", "transcribed_plasmids_parsed.csv")
CONSTRUCTS <- c("parent_bidir", "negative_control_3d", "3d-adapted", "3d-anti_adapted")
INTRON <- paste0("GTGAGTACTCCCTCTCAAAAGCGGGCATGACTTCTGCGCTAAGATTGTCAGTTTCCAAAAACGAGGAGGATTTGATAT",
                 "TCACCTGGCCCGCGGTGATGCCTTTGAGGGTGGCCGCGTCCATCTGGTCAGAAAAGACAATCTTTTTGTTGTCAGCTT",
                 "GAGGTGTGGCAGGCTTGAGATCTGGCCATACACTTGAGTGACAATGACATCCACTTTGCCTTTCTCTCCACAG")
PAS <- c("AATAAA", "ATTAAA")

# --- firefly ------------------------------------------------------------------------
# The spreadsheet row is "firefly luc"; the EMITTED id drops the space. Construct ids become
# FASTA names, data-frame rownames and plot labels downstream, and RNAplfold writes
# <name>_lunp, so an id with a space produces a file the reader can never find.
FIREFLY    <- "firefly luc"   # spreadsheet key
FIREFLY_ID <- "firefly"       # emitted id
# Promega's annotated 5'UTR intron for this vector. Located and checked below, never trusted.
FF_INTRON <- paste0("GTGAGTACTCCCTCTCAAAAGCGGGCATGACTTCTGCGCTAAGATTGTCAGTTTCCAAAAACGAGGAGGATTTGATAT",
                    "TCACCTGGCCCGCGGTGATGCCTTTGAGGGTGGCCGCGTCCATCTGGTCAGAAAAGACAATCTTTTTGTTGTCAAGTT",
                    "GAGGTGTGGCAGGCTTGAGATCTGGCCATACACTTGAGTGACAATGACATCCACTTTGCCTTTCTCTCCA")
# The CMV landmark that spans the transcription start site. The TSS is the G of TSS_MOTIF.
LANDMARK  <- "GTTTAGTGAACCGTCAGAT"
TSS_MOTIF <- "GTCAGAT"
# What firefly lost at the 3' splice site, relative to the NanoLuc backbone.
CAGG      <- "CAGG"

clean <- function(x) toupper(gsub("[^A-Za-z]", "", x))
codon_vec <- function(s) substring(s, seq(1, nchar(s) - 2, 3), seq(3, nchar(s), 3))
n_hits <- function(pat, s) sum(gregexpr(pat, s, fixed = TRUE)[[1]] > 0)

raw <- read_excel(XLSX) %>% transmute(construct, tx = clean(transcribed_seq), cds = clean(CDS))
stopifnot("construct names changed" = all(CONSTRUCTS %in% raw$construct),
          "sequences contain non-ACGT characters" =
            all(strsplit(paste0(raw$tx, raw$cds, collapse = ""), "")[[1]] %in% c("A","C","G","T")))
pl <- raw %>% filter(construct %in% CONSTRUCTS)

# --- splice, then split into regions -------------------------------------------------
stopifnot("the intron does not obey the GT..AG rule" =
            substr(INTRON, 1, 2) == "GT" && substr(INTRON, nchar(INTRON) - 1, nchar(INTRON)) == "AG")
pl <- pl %>%
  rowwise() %>%
  mutate(intron_start = { h <- gregexpr(INTRON, tx, fixed = TRUE)[[1]]
                          stopifnot("the intron does not occur exactly once in this transcript" =
                                      length(h) == 1 && h[1] > 0); h[1] },
         pre_cds_start = { h <- gregexpr(cds, tx, fixed = TRUE)[[1]]
                           stopifnot("CDS does not occur exactly once in this transcript" =
                                       length(h) == 1 && h[1] > 0); h[1] },
         # The intron must lie wholly within the 5'UTR, or splicing would alter the CDS.
         ok_in_utr5 = intron_start + nchar(INTRON) - 1L < pre_cds_start,
         tx = paste0(substr(tx, 1, intron_start - 1L),
                     substr(tx, intron_start + nchar(INTRON), nchar(tx))),
         cds_start = { h <- gregexpr(cds, tx, fixed = TRUE)[[1]]
                       stopifnot("CDS is not unique after splicing" = length(h) == 1 && h[1] > 0); h[1] },
         utr5 = substr(tx, 1, cds_start - 1L),
         utr3 = substr(tx, cds_start + nchar(cds), nchar(tx)),
         protein = as.character(translate(DNAString(cds)))) %>%
  ungroup()
stopifnot("the intron is not contained in the 5'UTR - splicing would change the CDS" = all(pl$ok_in_utr5),
          "splicing did not remove exactly the intron length" =
            all(pl$pre_cds_start - pl$cds_start == nchar(INTRON)),
          "the intron sequence survives in a spliced transcript" = all(map_int(pl$tx, ~ n_hits(INTRON, .x)) == 0))
cat(sprintf("spliced out a %d nt 5'UTR intron from %d constructs (5'UTR %d -> %d nt)\n",
            nchar(INTRON), nrow(pl), pl$pre_cds_start[1] - 1L, nchar(pl$utr5[1])))

# --- structural checks ---------------------------------------------------------------
# These must hold for ANY construct - the four or firefly - so they are a function called on
# both rather than a block that silently only ever sees the NanoLuc set.
universal_checks <- function(d, who) {
  stopifnot(
    "regions do not reconstitute the transcript" = all(paste0(d$utr5, d$cds, d$utr3) == d$tx),
    "CDS length is not a multiple of 3"          = all(nchar(d$cds) %% 3 == 0),
    "CDS does not start with ATG"                = all(substr(d$cds, 1, 3) == "ATG"),
    "CDS does not end in a stop codon"           = all(substr(d$cds, nchar(d$cds) - 2, nchar(d$cds)) %in% c("TAA","TAG","TGA")),
    "CDS contains an internal stop codon"        = all(map_int(d$cds, ~ sum(head(codon_vec(.x), -1) %in% c("TAA","TAG","TGA"))) == 0),
    "a 5'UTR or 3'UTR is empty"                  = all(nchar(d$utr5) > 0), all(nchar(d$utr3) > 0),
    "a 3'UTR has no polyadenylation signal"      = all(map_lgl(d$utr3, ~ any(map_int(PAS, n_hits, s = .x) > 0))))
  cat(sprintf("  universal structural checks passed: %s (%d)\n", who, nrow(d)))
  invisible(TRUE)
}
universal_checks(pl, "NanoLuc four")

stopifnot(
  # These four are synonymous variants of one reporter in one backbone: everything but the
  # codons must be shared, or they are not the controlled comparison they are taken to be.
  # Scoped to the four ON PURPOSE - firefly satisfies none of them.
  "CDS lengths differ"                         = n_distinct(nchar(pl$cds)) == 1,
  "variants do not encode the same protein"    = n_distinct(pl$protein) == 1,
  "variants do not share one 5'UTR"            = n_distinct(pl$utr5) == 1,
  "variants do not share one 3'UTR"            = n_distinct(pl$utr3) == 1,
  "the negative control is not a synonymous shuffle of the parent" =
    identical(sort(codon_vec(pl$cds[pl$construct == "negative_control_3d"])),
              sort(codon_vec(pl$cds[pl$construct == "parent_bidir"]))))

# --- biological spot-check: the designed GC3 ordering, and the designed CDS themselves -
gc3 <- function(s) { v <- codon_vec(s); 100 * mean(substr(v, 3, 3) %in% c("G", "C")) }
pl$gc3_pct <- map_dbl(pl$cds, gc3)
g <- set_names(pl$gc3_pct, pl$construct)
stopifnot("designed GC3 ordering is wrong: adapted < parent < anti-adapted expected" =
            g[["3d-adapted"]] < g[["parent_bidir"]] && g[["parent_bidir"]] < g[["3d-anti_adapted"]],
          "the negative control should match the parent's GC3 (same codon composition)" =
            abs(g[["negative_control_3d"]] - g[["parent_bidir"]]) < 1e-9)
# Independent check against the file 64_nanoluc_codon_adaptation.Rmd wrote: each CDS must still
# begin with its designed 513 nt variant, so a spreadsheet edit cannot silently change a sequence.
des <- read_csv(here("output", "predictive_modeling", "nanoluc_codon_designed_sequences.csv"),
                show_col_types = FALSE)
dmap <- c(parent_bidir = "nanoluc_original", `3d-adapted` = "nanoluc_3d_adapted",
          `3d-anti_adapted` = "nanoluc_3d_antiadapted", negative_control_3d = "nanoluc_3d_negctrl")
stopifnot("a CDS no longer starts with its designed sequence from notebook 64" =
            all(map_lgl(pl$construct, ~ startsWith(pl$cds[pl$construct == .x],
                                                   des$bases[des$name == dmap[[.x]]]))))

# --- firefly: locate the TSS, splice the annotated intron, split into regions ---------
ff_raw <- raw %>% filter(construct == FIREFLY)
stopifnot("the firefly row is missing from the spreadsheet" = nrow(ff_raw) == 1)
ff_entry <- ff_raw$tx[1]; ff_cds <- ff_raw$cds[1]

# The TSS, from the CMV landmark. What licenses this motif is that the four NanoLuc transcripts
# are entered FROM their TSS and so begin with it - asserted, not assumed.
stopifnot("the NanoLuc 5'UTRs do not begin at the TSS motif - the landmark no longer anchors the start" =
            all(startsWith(pl$utr5, TSS_MOTIF)))
lm_h <- gregexpr(LANDMARK, ff_entry, fixed = TRUE)[[1]]
stopifnot("the CMV landmark does not occur exactly once in the firefly entry" =
            length(lm_h) == 1 && lm_h[1] > 0)
ff_tss <- lm_h[1] + (regexpr(TSS_MOTIF, LANDMARK, fixed = TRUE)[1] - 1L)
ff_unspliced <- substr(ff_entry, ff_tss, nchar(ff_entry))

ff_cds_h <- gregexpr(ff_cds, ff_unspliced, fixed = TRUE)[[1]]
stopifnot("the firefly CDS does not occur exactly once downstream of the TSS" =
            length(ff_cds_h) == 1 && ff_cds_h[1] > 0)
ff_u5_unspliced <- substr(ff_unspliced, 1, ff_cds_h[1] - 1L)
ff_u3           <- substr(ff_unspliced, ff_cds_h[1] + nchar(ff_cds), nchar(ff_unspliced))
stopifnot(
  "the firefly 5'UTR does not end on the gccacc Kozak - the CDS start is misplaced" =
    endsWith(ff_u5_unspliced, "GCCACC"),
  "the firefly unspliced 5'UTR is not 433 nt" = nchar(ff_u5_unspliced) == 433L)

# The annotated intron, located rather than trusted.
ff_i <- gregexpr(FF_INTRON, ff_u5_unspliced, fixed = TRUE)[[1]]
stopifnot(
  "the annotated firefly intron does not occur exactly once in its unspliced 5'UTR" =
    length(ff_i) == 1 && ff_i[1] > 0,
  "the annotated firefly intron does not begin with a GT donor" =
    startsWith(FF_INTRON, "GT"),
  "the annotated firefly intron must lie wholly inside the 5'UTR" =
    ff_i[1] + nchar(FF_INTRON) - 1L < nchar(ff_u5_unspliced),
  "the annotated firefly intron is identical to the NanoLuc intron - wrong sequence pasted in" =
    FF_INTRON != INTRON)

# The CAGG deletion is ASSERTED, not excepted: firefly's intron tail plus CAGG must reproduce
# the NanoLuc tail plus its first exon base. This is the whole reason the acceptor is CA, and a
# discriminative check - a different 3'-end edit would not satisfy it.
W <- 20L
nl_h <- gregexpr(INTRON, ff_raw_nl <- raw$tx[raw$construct == "parent_bidir"], fixed = TRUE)[[1]][1]
stopifnot(
  "firefly's 3' end is not the NanoLuc tail minus CAGG" =
    paste0(substr(FF_INTRON, nchar(FF_INTRON) - W + 1L, nchar(FF_INTRON)), CAGG) ==
    substr(ff_raw_nl, nl_h + nchar(INTRON) - (W + 3L), nl_h + nchar(INTRON)),
  "the annotated firefly intron ends in AG - the CAGG deletion is absent" =
    substr(FF_INTRON, nchar(FF_INTRON) - 1L, nchar(FF_INTRON)) != "AG")

ff_u5 <- paste0(substr(ff_u5_unspliced, 1, ff_i[1] - 1L),
                substr(ff_u5_unspliced, ff_i[1] + nchar(FF_INTRON), nchar(ff_u5_unspliced)))
ff <- tibble(construct = FIREFLY_ID, tx = paste0(ff_u5, ff_cds, ff_u3), cds = ff_cds,
             utr5 = ff_u5, utr3 = ff_u3,
             protein = as.character(translate(DNAString(ff_cds))),
             gc3_pct = gc3(ff_cds))
stopifnot(
  "splicing did not remove exactly the annotated intron length" =
    nchar(ff_u5_unspliced) - nchar(ff_u5) == nchar(FF_INTRON),
  "the annotated intron survives in the spliced firefly transcript" =
    n_hits(FF_INTRON, ff$tx) == 0,
  "the firefly mature 5'UTR is not 207 nt" = nchar(ff_u5) == 207L)
universal_checks(ff, "firefly")

# Biological spot-check: it must be luc2P - luc2 plus the mouse ODC PEST degron - not the
# PEST-less luc2, which is a prefix of it and passes every length and frame check.
stopifnot(
  "firefly does not encode the luc2 N-terminus" =
    startsWith(ff$protein, "MEDAKNIKKGPAPFYPLEDGTAGEQLHKAMKRYALVPGTIAFTDAHIE"),
  "firefly is missing the C-terminal ODC PEST degron (luc2, not luc2P?)" =
    endsWith(ff$protein, "SHGFPPEVEEQAAGTLPMSCAQESGMDRHPAACASARINV*"))
cat(sprintf("firefly: TSS at %d (%d nt upstream trimmed), 5'UTR %d -> %d nt, CDS %d nt, 3'UTR %d nt\n",
            ff_tss, ff_tss - 1L, nchar(ff_u5_unspliced), nchar(ff_u5),
            nchar(ff_cds), nchar(ff_u3)))

out <- bind_rows(pl, ff) %>%
  transmute(construct, utr5, cds, utr3, transcript = tx,
            utr5_nt = nchar(utr5), cds_nt = nchar(cds), utr3_nt = nchar(utr3),
            total_nt = nchar(tx), codons = nchar(cds) / 3, gc3_pct = round(gc3_pct, 2),
            cds_gc_pct = round(100 * map_dbl(cds, ~ mean(strsplit(.x, "")[[1]] %in% c("G","C"))), 2),
            source = basename(XLSX),
            note = if_else(construct == FIREFLY_ID,
                           sprintf("%d nt 5'UTR intron spliced out (Promega annotation); %d nt upstream of the TSS trimmed; CAGG deletion at the 3' splice site",
                                   nchar(FF_INTRON), ff_tss - 1L),
                           sprintf("%d nt 5'UTR intron spliced out", nchar(INTRON))))
stopifnot("expected five constructs" = nrow(out) == 5,
          # ids become FASTA names and filenames downstream.
          "a construct id contains whitespace" = !any(grepl("\\s", out$construct)))
write_csv(out, OUT)
print(as.data.frame(out %>% dplyr::select(construct, utr5_nt, cds_nt, utr3_nt, total_nt, gc3_pct, cds_gc_pct)))
cat("\nall checks passed; wrote", OUT, "\n")
