# 33a_parse_plasmid_sequences.R - parse accessories/transcribed_plasmids.xlsx into the canonical
# per-construct 5'UTR / CDS / 3'UTR table that 33_ turns into model features.
#
# The spreadsheet is Kate's record and is never modified; this writes a derived CSV.
#
# SCOPE: the four synonymous NanoLuc constructs only. The firefly entry is parsed and checked but
# deliberately NOT emitted - see "firefly" below.
#
# SPLICING. These transcripts carry a 229 nt intron in the 5'UTR (sequence supplied by Kate,
# 2026-10-06; verified here as an exact single-copy match at 236-464 in all four, obeying GT..AG).
# The models describe mature mRNA, so the intron is removed: the 5'UTR goes 554 -> 325 nt. The
# intron lies wholly inside the 5'UTR, which is asserted, so the CDS and 3'UTR are untouched.
#
# FIREFLY (excluded, 2026-10-06). Two unresolved issues, neither affecting the NanoLuc set:
#   1. Its transcript as entered begins 644 nt upstream of the transcription start site, inside the
#      CMV promoter. The start is locatable from the CMV landmark "GTTTAGTGAACCGTCAGAT" - the
#      "GTCAGAT" inside it is exactly where the NanoLuc transcripts begin - giving a 433 nt 5'UTR.
#      Corroborated independently: the 90 nt by which its intron's upstream flank is shorter than
#      the NanoLuc's equals the 90 nt difference in their TSS-to-shared-element linkers.
#   2. Its copy of the intron is missing CAGG at the 3' splice site, deleting the acceptor AG. It
#      therefore either retains the intron (5'UTR 433 nt) or splices to the next AG at 1027
#      (5'UTR 194 nt). Sequence alone cannot decide; it needs RNA-seq, RT-PCR or a sequenced
#      transcript. To include it, resolve that and add it to CONSTRUCTS with its own intron coords.
#   Rscript code/translational_control_overview/33a_parse_plasmid_sequences.R
suppressPackageStartupMessages({ library(tidyverse); library(readxl); library(here); library(Biostrings) })

XLSX <- here("accessories", "transcribed_plasmids.xlsx")
OUT  <- here("accessories", "transcribed_plasmids_parsed.csv")
CONSTRUCTS <- c("parent_bidir", "negative_control_3d", "3d-adapted", "3d-anti_adapted")
INTRON <- paste0("GTGAGTACTCCCTCTCAAAAGCGGGCATGACTTCTGCGCTAAGATTGTCAGTTTCCAAAAACGAGGAGGATTTGATAT",
                 "TCACCTGGCCCGCGGTGATGCCTTTGAGGGTGGCCGCGTCCATCTGGTCAGAAAAGACAATCTTTTTGTTGTCAGCTT",
                 "GAGGTGTGGCAGGCTTGAGATCTGGCCATACACTTGAGTGACAATGACATCCACTTTGCCTTTCTCTCCACAG")
PAS <- c("AATAAA", "ATTAAA")

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
stopifnot(
  "regions do not reconstitute the transcript" = all(paste0(pl$utr5, pl$cds, pl$utr3) == pl$tx),
  "CDS length is not a multiple of 3"          = all(nchar(pl$cds) %% 3 == 0),
  "CDS does not start with ATG"                = all(substr(pl$cds, 1, 3) == "ATG"),
  "CDS does not end in a stop codon"           = all(substr(pl$cds, nchar(pl$cds) - 2, nchar(pl$cds)) %in% c("TAA","TAG","TGA")),
  "CDS contains an internal stop codon"        = all(map_int(pl$cds, ~ sum(head(codon_vec(.x), -1) %in% c("TAA","TAG","TGA"))) == 0),
  "a 5'UTR or 3'UTR is empty"                  = all(nchar(pl$utr5) > 0), all(nchar(pl$utr3) > 0),
  "a 3'UTR has no polyadenylation signal"      = all(map_lgl(pl$utr3, ~ any(map_int(PAS, n_hits, s = .x) > 0))),
  # These four are synonymous variants of one reporter in one backbone: everything but the
  # codons must be shared, or they are not the controlled comparison they are taken to be.
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

out <- pl %>%
  transmute(construct, utr5, cds, utr3, transcript = tx,
            utr5_nt = nchar(utr5), cds_nt = nchar(cds), utr3_nt = nchar(utr3),
            total_nt = nchar(tx), codons = nchar(cds) / 3, gc3_pct = round(gc3_pct, 2),
            cds_gc_pct = round(100 * map_dbl(cds, ~ mean(strsplit(.x, "")[[1]] %in% c("G","C"))), 2),
            source = basename(XLSX), note = sprintf("%d nt 5'UTR intron spliced out", nchar(INTRON)))
write_csv(out, OUT)
print(as.data.frame(out %>% dplyr::select(construct, utr5_nt, cds_nt, utr3_nt, total_nt, gc3_pct, cds_gc_pct)))
cat("\nall checks passed; wrote", OUT, "\n")
