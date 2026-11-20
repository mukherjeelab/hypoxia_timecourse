# Reporter controls: eIF3d-sensitive and eIF3d-insensitive 5'UTRs

A running list of mRNAs whose 5'UTR may make translation depend on eIF3d (positive controls),
and of mRNAs whose translation does not depend on eIF3d (negative controls, in their own section
at the end), for reporter experiments. Newest entries go at the bottom of each section. Started
2026-11-17. Candidates checked and rejected are in `CHECKED_NOT_SELECTED.md`.

## Definitions used throughout

- **eIF3d target:** translation efficiency (TE = ribosome footprints / mRNA) drops on eIF3d
  knockdown, TE log2 fold change (sictrl / si3d) > 0.5 and padj < 0.05. This is the same rule
  behind `output/genesets/hypoxia_3d_promotes_TE_1hr_lfc0.5.csv` (MDA-MB-231) and
  `mcf7six1_hypoxia_3d_promotes_TE_1hr_lfc0.5.csv` (MCF7-SIX1). Positive LFC = eIF3d promotes.
- **Category** (from `categorize_translation_changes()`, significance = padj < 0.05, no LFC
  cutoff): **Exclusive** = footprints and TE change, mRNA does not (a purely translational
  effect, the best kind for a 5'UTR reporter). **Buffered** = mRNA and TE both change, so the
  protein-level effect is smaller than the TE number suggests.
- **Expressed in our 231 / MCF7-SIX1 cells:** the most abundant transcript of the gene in
  that cell line's siCTRL normoxia RNA-seq, as used in our feature matrices.
- **Caveat that applies to every entry:** TE is measured per gene, not per transcript. It
  reflects the mix of isoforms the cells express, which is mostly the expressed transcript
  above. A rare isoform's 5'UTR is not directly tested in our data.
- **MANE Select:** the single reference transcript per gene agreed by NCBI and Ensembl, read
  from `accessories/human/mane_select_v49.csv` (extracted from the GENCODE v49 annotation).
- Length and GC come from GENCODE v49 on hg38, via `transcript_length_gc.R`
  (output: `output/finding_positive_controls/transcript_length_gc.csv`).

---

# Positive controls

## BCL2L12

**Status: strongest candidate from set 1. Both 5'UTRs listed below.**

eIF3d target and **Exclusive** in hypoxia 1hr in **both** cell lines:

| cell line | TE LFC (padj) | RNA LFC (padj) | ribo LFC (padj) | category |
|---|---|---|---|---|
| MDA-MB-231 | 0.89 (0.007) | 0.26 (0.09) | 1.27 (4e-6) | Exclusive |
| MCF7-SIX1 | 0.72 (2e-5) | -0.15 (0.26) | 0.60 (2e-4) | Exclusive |

Not significant at 4hr or 24hr hypoxia in MDA-MB-231.

### Transcripts

| | ENST00000246784 | ENST00000441864 |
|---|---|---|
| role | expressed in our 231 and MCF7-SIX1 cells | tested in the Teleman data (per Kate; set 1 list) |
| MANE Select? | **yes** | no |
| share of gene's mRNA, 231 / MCF7-SIX1 (hypoxia 1hr siCTRL) | 42.5% / 60.9% | 0.8% / 2.9% |
| transcript length, GC | 1,064 nt, 61.2% | 1,371 nt, 62.5% |
| 5'UTR length, GC | 143 nt, 62.2% | 510 nt, 60.2% |
| CDS length, GC | 753 nt, 65.2% | 750 nt, 65.3% |
| 3'UTR length, GC | 168 nt, 42.3% | 111 nt, 54.1% |
| upstream ATGs in 5'UTR | 0 | 5 |

**How the two 5'UTRs relate.** The 143-nt 5'UTR expressed in our cells is identical to the last 143 nt of the
510-nt one. The 510-nt version adds 367 nt at its 5' end (exon 1 starts at chr19:49,665,566
instead of 49,665,933; exon 2 is the same 8 nt). Most expressed BCL2L12 isoforms (111-146 nt)
share this same region next to the start codon.

**All 5 upstream ATGs are in the extra 367 nt.** Three start short uORFs that stop inside the
5'UTR (36, 19, 9 codons). Two (252 and 174 nt upstream) are in frame with the CDS with no stop
before the CDS start codon, so they would produce N-terminally extended BCL2L12. In a reporter, these
would change the readout independently of eIF3d. None of these ATGs are in the 143-nt
5'UTR expressed in our cells.

**Note on the Teleman source.** Our own Teleman HeLa transcript map
(`output/predictive_modeling/precomputed_teleman_tx.rds`, the dominant isoform in Roiuk control
total RNA) picks **ENST00000246784** for BCL2L12, not ENST00000441864. So the 510-nt transcript
is not the dominant HeLa isoform in that RNA-seq; "tested in the Teleman data" refers to
whichever Teleman dataset the set 1 list came from, which is not recorded here yet.

The two transcripts also differ downstream (CDS 3 nt shorter, different 3'UTR), which
matters only if more than the 5'UTR is cloned.

**Source:** `01_eif3d_hypoxia_target_lookup.Rmd`, `00_mcf7six1_hypoxia_si3d_categories.Rmd`,
`02_bcl2l12_5utr_check.Rmd` (this folder).

---

## TPRA1

**Status: MCF7-SIX1 only, Buffered. The tested 5'UTR is the one MCF7-SIX1 expresses, not the one
231 expresses.**

eIF3d target in hypoxia 1hr in MCF7-SIX1 only; just misses in MDA-MB-231:

| cell line | TE LFC (padj) | RNA LFC (padj) | ribo LFC (padj) | category |
|---|---|---|---|---|
| MDA-MB-231 | 0.52 (0.097) | 0.003 (0.99) | 0.61 (0.051) | Not_sig |
| MCF7-SIX1 | 0.81 (2e-4) | -0.46 (0.002) | 0.41 (0.12) | Buffered |

**Buffered** in MCF7-SIX1: the mRNA *rises* on eIF3d knockdown while TE falls, and footprints
do not change significantly, so the protein-level effect is likely smaller than the TE number
suggests. Not significant at 4hr or 24hr hypoxia in MDA-MB-231.

**Low abundance:** the whole gene is ~3 TPM in both lines (hypoxia 1hr siCTRL: 231 2.9, MCF7-SIX1
3.5), vs ~16-19 TPM for BCL2L12, so the transcript shares below rest on few reads.

### Transcripts

| | ENST00000296210 | ENST00000355552 | ENST00000965987 |
|---|---|---|---|
| role | tested in the Teleman data (per Kate; set 1 list) | expressed in our MCF7-SIX1 cells | expressed in our 231 cells |
| MANE Select? | no | **yes** | no |
| share of gene's mRNA, 231 / MCF7-SIX1 (hypoxia 1hr siCTRL) | 10.4% / 29.8% | 2.5% / 25.4% | 38.7% / not quantified |
| transcript length, GC | 1,725 nt, 62.5% | 3,752 nt, 56.3% | 1,652 nt, 59.8% |
| 5'UTR length, GC | 315 nt, 75.9% | 341 nt, 74.2% | 82 nt, 62.2% |
| CDS length, GC | 882 nt, 59.6% | 1,122 nt, 59.7% | 1,122 nt, 59.7% |
| 3'UTR length, GC | 528 nt, 59.3% | 2,289 nt, 52.0% | 448 nt, 59.6% |
| upstream ATGs in 5'UTR | 1 | 1 | 1 |

In MCF7-SIX1 the gene's mRNA is split roughly evenly between the tested transcript (29.8%), the
MANE Select (25.4%) and ENST00000450633 (31.4%, 291-nt 5'UTR, not compared here). The
"expressed in MCF7-SIX1" pick (made from normoxia RNA-seq) is therefore not dominant in hypoxia.

**How the 5'UTRs relate.**
- **Tested vs MANE Select (expressed in MCF7-SIX1):** the tested 315-nt 5'UTR is identical to the
  last 315 nt of the MANE 341-nt one, which adds only 26 nt at its 5' end (exon 1 extends to
  chr3:127,590,733 instead of 127,590,707; minus strand). Both are very GC-rich (74-76%). Both
  have the same single upstream ATG, 72 nt before the CDS start codon, in frame, starting a
  12-codon uORF that stops inside the 5'UTR.
- **Tested vs expressed in 231:** different leaders. The 231 transcript uses a different first
  exon about 1.2 kb upstream (chr3:127,591,888-127,591,952); only the 17-nt exon 2 next to the
  start codon is shared. Its one upstream ATG (61 nt before the CDS start, out of frame) starts
  an ORF that runs past the CDS start codon (134 codons), overlapping the main ORF.
- So the tested 5'UTR represents what MCF7-SIX1 expresses, where TPRA1 is a target, but not what
  231 expresses.

**The tested transcript also encodes a shorter protein:** CDS 882 nt vs 1,122 nt for both other
transcripts (different splicing downstream of the start codon). This matters only if more than
the 5'UTR is cloned.

**Source:** `01_eif3d_hypoxia_target_lookup.Rmd`, `00_mcf7six1_hypoxia_si3d_categories.Rmd`,
`02_bcl2l12_5utr_check.Rmd` knit with `symbol = "TPRA1"` -> `02_TPRA1_5utr_check.html`
(this folder).

---

## Set 1, other genes (for reference)

Set 1 = the 19-transcript "5' UTR" list, `set1_5utr_transcripts.csv`. Other eIF3d targets in
hypoxia 1hr, all MCF7-SIX1 only (TPRA1 has its own entry above). Their isoforms and 5'UTRs have **not** been checked yet.

| gene | MCF7-SIX1 TE LFC (padj) | category | given = transcript expressed in our cells? | given = MANE Select? | note |
|---|---|---|---|---|---|
| CCDC85B | 0.63 (0.01) | Exclusive | yes, both lines | yes | not significant in 231 (TE -0.45) |
| HNRNPA0 | 0.56 (3e-9) | Buffered | yes, both lines | yes | mRNA rises on knockdown |
| CDKN2B | 0.45 (0.02) | Exclusive | no (MCF7-SIX1) | no (MANE = ENST00000276925) | below the 0.5 cutoff; absent from 231 results |

---

# Negative controls

A good negative control: TE flat on eIF3d knockdown (project rule: |TE LFC| < 0.5 and
padj > 0.3) in both cell lines, well expressed, and a simple 5'UTR.

## ACTB

**Status: selected as the negative control (Kate, 2026-11-18). Not neutral at 231 hypoxia 24hr.**

Named in Roiuk et al. 2024 (*Nat Commun* 15:6692, doi 10.1038/s41467-024-51027-z, Fig. 5E) as
an mRNA not affected by eIF3d knockdown. That paper's reporter negative control was the HBB
5'UTR, which our cells do not express (0 reads in both lines), so our data cannot speak to HBB.

TE flat in every eIF3d knockdown contrast except one:

| contrast | TE LFC (padj) |
|---|---|
| MDA-MB-231 hypoxia 1hr | 0.06 (0.82) |
| MDA-MB-231 hypoxia 4hr | -0.43 (0.44) |
| **MDA-MB-231 hypoxia 24hr** | **-1.14 (0.009), eIF3d inhibits** |
| MDA-MB-231 normoxia 1+4hr (time-regressed) | 0.08 (0.73) |
| MDA-MB-231 normoxia 1hr / 4hr / 24hr | 0.04 / 0.09 / 0.18 (all padj >= 0.50) |
| MCF7-SIX1 hypoxia 1hr | -0.05 (0.70) |
| MCF7-SIX1 normoxia 1hr | -0.11 (0.55) |
| HeLa normoxia (steady state) | -0.20 (0.59) |

In MDA-MB-231 hypoxia 1hr ACTB mRNA falls on knockdown (RNA LFC 0.95) together with footprints
(Forwarded), but TE, the quantity a 5'UTR reporter reports, does not change. Abundant: ~600-800
TPM.

### Transcripts

| | ENST00000646664 | ENST00000331789 (GENCODE v26 only) |
|---|---|---|
| role | MANE Select; expressed in our 231 cells | expressed in our MCF7-SIX1 cells (99% of ACTB in the v26 quantification) |
| MANE Select? | **yes** | no (retired before GENCODE v49) |
| transcript length, GC | 1,812 nt, 55.2% | 1,917 nt (v26 model) |
| 5'UTR length, GC | 84 nt, 76.2% | 192 nt, 72.4% (v26 model) |
| CDS length, GC | 1,128 nt, 60.5% | same start codon (chr7:5,529,655-5,529,657) |
| 3'UTR length, GC | 600 nt, 42.2% | not computed |
| upstream ATGs in 5'UTR | 0 | 1, in the extra 108 nt |

**These are the same mRNA leader.** MCF7-SIX1 RNA-seq was quantified against GENCODE v26, where
ACTB's dominant transcript is ENST00000331789. Its 5'UTR is the MANE 84 nt plus 108 nt further
upstream (exon 1 to chr7:5,530,709 vs 5,530,601; minus strand). Those 108 nt contain ACTB's
CCAAT box (`GGCCAATCAG`) and TATA box (`TATAAAA`, about 25 nt upstream of the MANE start), so
they are promoter, mis-annotated as transcribed in v26. The real leader in both lines is the
84-nt MANE 5'UTR. (The MCF7-SIX1 feature matrix lists ENST00000425660 for ACTB, which carries
0.1% of ACTB mRNA; that is an artifact of the v26 / v49 mismatch, not biology.)

**Source:** `01_eif3d_hypoxia_target_lookup.Rmd` knit on `set4_roiuk2024_controls.csv`
(-> `01_eif3d_hypoxia_target_lookup_set4.html`); `02_bcl2l12_5utr_check.Rmd` with
`symbol = "ACTB"` (-> `02_ACTB_5utr_check.html`); all-contrast TE from the
`output/translation_categories_si3d_vs_sictrl_*.csv` files; v26 model from the GENCODE release 26
annotation (ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/).
