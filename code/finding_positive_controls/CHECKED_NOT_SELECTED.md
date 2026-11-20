# Candidates checked and not selected as positive controls

Transcripts checked for eIF3d-dependent translation in our hypoxia data and **not** added to
`POSITIVE_CONTROLS.md` (Kate's decision). Recorded so they are not re-checked, and so the
results can be found again. Definitions (eIF3d target, Exclusive, Buffered, MANE Select) are
at the top of `POSITIVE_CONTROLS.md`. Values are TE log2 fold change (padj), sictrl / si3d,
positive = eIF3d promotes, hypoxia 1hr.

## Set 2: JUN and BTG1 (checked 2026-11-17, not selected 2026-11-18)

RefSeq IDs given by Kate, mapped through NCBI MANE v1.5
(`MANE.GRCh38.v1.5.summary.txt.gz`, ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/),
which pairs one RefSeq with one Ensembl transcript; agrees with
`accessories/human/mane_select_v49.csv`.

| gene | RefSeq | Ensembl (MANE Select) | MDA-MB-231 | MCF7-SIX1 |
|---|---|---|---|---|
| JUN | NM_002228.4 | ENST00000371222.4 | 0.02 (0.94), not significant | **-0.62 (1e-4), eIF3d inhibits**, Buffered; ribo -0.97 (2e-11) |
| BTG1 | NM_001731.3 | ENST00000256015.5 | 0.13 (0.66), not significant | **0.85 (1e-5), target**, Exclusive; ribo 0.76 (4e-4) |

Neither is significant at 4hr or 24hr hypoxia in MDA-MB-231.

- **JUN** goes the opposite way from a positive control in MCF7-SIX1 (footprints rise on eIF3d
  knockdown). The given transcript is the only one expressed (>= 99% of JUN mRNA in both
  lines). 5'UTR 977 nt, 65.0% GC, no upstream ATGs; transcript 3,257 nt, 56.1% GC.
- **BTG1** is an eIF3d target in MCF7-SIX1 only, where the given transcript is 94% of BTG1 mRNA.
  5'UTR 311 nt, 65.6% GC, 3 upstream ATGs (uORFs of 24, 47, 2 codons); transcript 4,629 nt,
  38.6% GC. **Open question in 231:** Salmon assigns 97% of BTG1 reads to ENST00000673901
  (`protein_coding_CDS_not_defined`; own first exon inside intron 1, shares part of exon 2,
  lacks the start codon), which is why BTG1 is missing from the 231 feature matrix. Real or a
  quantification artifact is unresolved; read coverage over BTG1 exon 1 would settle it.

**Where the data is**
- Input list: `set2_jun_btg1.csv` (this folder)
- TE / category calls: `output/finding_positive_controls/set2_jun_btg1_eif3d_hypoxia_calls.csv`;
  report `01_eif3d_hypoxia_target_lookup_set2.html`
- Transcript and 5'UTR checks: `02_JUN_5utr_check.html`, `02_BTG1_5utr_check.html`
  (from `02_bcl2l12_5utr_check.Rmd` with `symbol` / `given_tx` params), isoform shares in
  `output/finding_positive_controls/{JUN,BTG1}_isoform_share_hypoxia_1hr.csv`
- Length / GC: `output/finding_positive_controls/transcript_length_gc.csv`

## Set 3: ALKBH5 (checked and not selected 2026-11-18)

Given transcript ENST00000399138 (Ensembl ID from Kate) is MANE Select and the transcript
expressed in both our 231 and MCF7-SIX1 cells.

| | MDA-MB-231 | MCF7-SIX1 |
|---|---|---|
| TE | 0.10 (0.71) | **0.50 (6e-4)**, target, right at the 0.5 cutoff |
| RNA | 0.23 (0.14) | -0.50 (1e-10), rises on eIF3d knockdown |
| ribo | 0.38 (0.14) | 0.03 (0.91), unchanged |
| category | Not_sig | Buffered |

Not significant at 4hr or 24hr hypoxia in MDA-MB-231.

- **Why not selected:** a marginal target in MCF7-SIX1 only, and the TE drop comes entirely from
  the mRNA rising while footprints stay flat, so protein output likely does not change.
- **Transcript:** 99.9% of ALKBH5 mRNA in MCF7-SIX1, 65% in 231 (the other 35% is
  ENST00000490106, `retained_intron`, non-coding). The only other annotated 5'UTR belongs to
  ENST00000541285, expressed at <= 0.2%.
- **Length / GC:** transcript 3,159 nt, 57.6%; 5'UTR 416 nt, 75.7%; CDS 1,185 nt, 63.3%;
  3'UTR 1,558 nt, 48.4%.
- **5'UTR:** 4 upstream ATGs, all starting uORFs that stop inside the 5'UTR (28, 6, 6, 12
  codons).

**Where the data is**
- Input list: `set3_alkbh5.csv` (this folder)
- TE / category calls: `output/finding_positive_controls/set3_alkbh5_eif3d_hypoxia_calls.csv`;
  report `01_eif3d_hypoxia_target_lookup_set3.html`
- Transcript and 5'UTR check: `02_ALKBH5_5utr_check.html`, isoform shares in
  `output/finding_positive_controls/ALKBH5_isoform_share_hypoxia_1hr.csv`
- Length / GC: `output/finding_positive_controls/transcript_length_gc.csv`

## Set 5: RPTOR (checked and not selected 2026-11-18)

Given transcript ENST00000306801 (Ensembl ID from Kate) is MANE Select and the transcript
expressed in both our 231 and MCF7-SIX1 cells.

- **Why not selected:** not eIF3d-dependent in any knockdown contrast in either cell line.

| contrast | TE LFC (padj) |
|---|---|
| MDA-MB-231 hypoxia 1hr / 4hr / 24hr | 0.17 (0.53) / -0.02 (0.98) / 0.00 (0.99) |
| MCF7-SIX1 hypoxia 1hr | 0.03 (0.91) |
| MCF7-SIX1 normoxia 1hr | 0.23 (0.26) |
| MDA-MB-231 normoxia 1+4hr (time-regressed) | 0.35 (0.23) |
| HeLa normoxia (steady state) | 0.30 (0.36) |

- **Expression:** low, ~1-2 TPM in both lines, so the TE estimate is noisy.
- **Length / GC:** transcript 6,838 nt, 60.0%; 5'UTR 792 nt, 60.2%; CDS 4,008 nt, 59.0%;
  3'UTR 2,038 nt, 62.0%.
- **5'UTR:** 6 upstream ATGs, all starting uORFs that stop inside the 5'UTR (8, 3, 3, 15, 4, 20
  codons).

**Where the data is**
- Input list: `set5_rptor.csv` (this folder)
- TE / category calls: `output/finding_positive_controls/set5_rptor_eif3d_hypoxia_calls.csv`;
  report `01_eif3d_hypoxia_target_lookup_set5.html`
- Transcript and 5'UTR check: `02_RPTOR_5utr_check.html`, isoform shares in
  `output/finding_positive_controls/RPTOR_isoform_share_hypoxia_1hr.csv`
- All-contrast TE: `output/translation_categories_si3d_vs_sictrl_*.csv`
- Length / GC: `output/finding_positive_controls/transcript_length_gc.csv`

## Set 6: ITGB1, CD101, CD103 / ITGAE (checked and not selected 2026-11-18)

No transcripts given; each gene's MANE Select used (`accessories/human/mane_select_v49.csv`).
CD103 is the integrin alphaE chain, gene symbol ITGAE.

| gene | MANE Select | 231 hypoxia 1hr | 231 hypoxia 4hr / 24hr | MCF7-SIX1 hypoxia 1hr |
|---|---|---|---|---|
| ITGB1 | ENST00000302278 | **0.74 (0.04), target, Exclusive**; ribo 0.98 (5e-4), RNA 0.10 (0.56) | 0.06 (0.91) / 0.08 (0.86) | 0.04 (0.84) |
| CD101 | ENST00000682167 | not expressed | not expressed | not expressed |
| ITGAE | ENST00000263087 | 0.38 (0.18), Forwarded | 0.03 (0.97) / -0.10 (0.71) | 0.17 (0.42), Forwarded |

- **ITGB1:** eIF3d target only in 231 hypoxia 1hr (also 231 normoxia 4hr, 0.74, padj 0.007;
  normoxia 1+4hr time-regressed 0.49, padj 0.12; HeLa -0.72, padj 0.14). Well expressed (~100 TPM
  231, ~45 MCF7-SIX1). In MCF7-SIX1 the MANE transcript is 68% of ITGB1 mRNA. **In 231 two
  leaders are nearly equal:** MANE 51%, ENST00000966600 41%, whose 26-nt 5'UTR shares only the
  last 2 nt before the start codon with the MANE one, so the gene-level 231 TE averages two
  different 5'UTRs. MANE 5'UTR: 86 nt, 80.2% GC, no upstream ATGs.
- **CD101 / ITGAE:** immune-cell genes. CD101 absent from both lines' TE results; ITGAE
  expressed but its changes are mRNA-level (Forwarded), TE not significant.

**Length / GC**

| transcript | whole | 5'UTR | CDS | 3'UTR |
|---|---|---|---|---|
| ITGB1 ENST00000302278 | 3,735 nt, 39.5% | 86 nt, 80.2% | 2,397 nt, 41.3% | 1,252 nt, 33.3% |
| CD101 ENST00000682167 | 3,542 nt, 50.2% | 58 nt, 46.6% | 3,066 nt, 50.0% | 418 nt, 52.4% |
| ITGAE ENST00000263087 | 3,803 nt, 54.6% | 44 nt, 70.5% | 3,540 nt, 55.3% | 219 nt, 40.2% |

**Where the data is**
- Input list: `set6_itgb1_cd101_itgae.csv` (this folder)
- TE / category calls: `output/finding_positive_controls/set6_itgb1_cd101_itgae_eif3d_hypoxia_calls.csv`;
  report `01_eif3d_hypoxia_target_lookup_set6.html`
- ITGB1 transcript and 5'UTR check: `02_ITGB1_5utr_check.html`, isoform shares in
  `output/finding_positive_controls/ITGB1_isoform_share_hypoxia_1hr.csv`
- All-contrast TE: `output/translation_categories_si3d_vs_sictrl_*.csv`
- Length / GC: `output/finding_positive_controls/transcript_length_gc.csv`
