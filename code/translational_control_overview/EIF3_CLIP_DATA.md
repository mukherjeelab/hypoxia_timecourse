# eIF3 CLIP datasets: inventory and tabled plan

Recorded 2026-09-30. **Tabled** - nothing built or modelled. This records what was verified so the
plan can resume without re-deriving it.

**Summary.** Three eIF3 CLIP studies are usable from deposited processed files (no raw-read
reprocessing). Only the Cate-lab NPC Quick-irCLIP gives a continuous per-gene signal; Lee 2015 and
De Silva 2021 report only genes that passed a cluster threshold. **Open decision:** what the outcome
is for the thresholded two (bound yes/no vs strength among bound), which Jurkat bands and states,
NPC region resolution, and where the heatmap goes.

Terms: **CLIP** has no input or IgG control, so "enrichment" here means CLIP signal per gene divided
by that gene's mRNA abundance in the same cells - crosslinks per mRNA molecule, relative to other
genes. **Band** - De Silva excised gel bands that contain crosslinked mixtures of subunits.

## Datasets

| Study | Cells | eIF3 subunit(s) | Deposited | RNA-seq denominator |
|---|---|---|---|---|
| Lee, Kranzusch, Cate 2015 *Nature* 522:111, GSE65004 | HEK293T | a, b, d, g (separate IPs, 3 replicates each) | `GSE65004_PARCLIP.txt.gz` = 1,477 clusters reproducible across replicates, reads + T>C, annotation, **hg19**. Identical to `mpra_preliminary/data/lee_parclip.csv` | none deposited; HEK293T RNA-seq in GSE236188 (`GSE236188_RNAseq_rpkm.txt.gz`, WT DMSO), or the subunit-seq inputs |
| De Silva et al. 2021 *eLife* 10:e74272, GSE191306 | Jurkat, +/- PMA + ionomycin 5 h | bands A/C/B, B/C, D/L, F/E, G/E/F | `..._PAR-CLIP-clusters.xlsx`: per-gene summaries, 20 sheets (5 bands x 2 states x 2 reps), **gene symbol only**; PARpipe min 7 reads/cluster; hg38. `..._global-stats.xlsx`. PARalyzer `act.clusters` / `non.clusters` (~118k clusters each, hg38) in `mpra_preliminary/data/jurkat_t_cell/` | `GSE191306_DeSilva-Jurkat-eIF3-RNASeq.xlsx`: transcript TPM, activated and non-activated |
| Mattick / Cate lab, *eLife* 2025 (bioRxiv 2023.11.11.566681), GSE246727 (+ APA-seq GSE246786) | WTC-11 NPCs, differentiated / undifferentiated | whole eIF3 (antibody not yet confirmed) | per-replicate CPM bigWigs GSM7876069-74 (3 diff + 3 undiff, 12-21 MB each); pooled bigWigs; HOMER peaks = `S2/S3_cate_2023.xlsx` (Peak Score all 0) | `GSE246727_JHDC001{G,H}-*_Transcriptome_Hits.xlsx` = **kallisto mRNA-seq** (counts + TPM), not CLIP |

Downloaded to the session scratchpad only (re-fetch from GEO): De Silva RNA-seq + global stats,
Lee GSE65004 table, Cate kallisto tables.

## Facts established

- **GEO sample-label correction (2025-12-02):** titles swapped between GSM5743374/6 and GSM5743375/7
  (replicate 2, activated vs non-activated raw reads). The processed clusters xlsx on GEO is
  md5-identical to the local copy (`637de73e...`). Its labels are **correct**: same-condition
  replicates agree best in every band (Spearman of per-gene ReadCountSum, e.g. band A activated
  0.91 / non-activated 0.85 vs 0.64 / 0.46 across conditions; B 0.87 / 0.82 vs 0.63 / 0.59).
- **Cited Lee cluster counts differ from the file:** a secondary source quotes a 328, b 264, d 356;
  the file has a 375, b 298, d 403, g 401. Reconcile against the paper's own table (probably an
  mRNA-only filter) before use.
- **Detection follows expression.** Lee-bound genes: median HEK293T input RPKM ~14 vs ~5 unbound
  (all four subunits); non-activated Jurkat bound genes: TPM 10 vs 4. Length barely differs
  between bound and unbound. Among bound genes, reads per mRNA rise with length (rho +0.29 to +0.45
  in Jurkat).
- **Activated-Jurkat RNA-seq looks abnormal:** only 2,776 matrix genes at TPM >= 1 vs 6,969
  non-activated; ARF5 17.8 -> 2.1 TPM. Must be understood before it serves as a denominator.
- **Cate headline:** eIF3 crosslinks mostly at 3'UTR termini next to the poly(A) tail, so
  whole-transcript NPC enrichment will largely mean 3'-end binding. bigWigs are labelled "BOTH" -
  check whether strands are merged.
- The matrix's existing `clip_*` features came from a string-parsed score card and are lossy (42 of
  285 eIF3d-bound genes missing) - rebuild from sources.

## Tabled plan

1. Prep notebook(s): copy sources to `accessories/eif3_clip/` with md5; replication checks
   (Lee counts; De Silva global stats, replicate labels, TCR TRAC/TRBC 3'UTR binding in activated;
   Cate CPM sums, 3'UTR-end concentration, peak counts); gene mapping (Lee Ensembl; De Silva
   symbol -> GENCODE v49 with ambiguity counts; Cate coverage on the matrix transcript's hg38
   exons); expression filter + denominators; outcome tables and bound / unbound gene sets.
2. Models: the regression binding mode (`outcome_kind = "binding"`), generalized beyond HEK293T to
   Jurkat and NPC; intrinsic features only.
3. Heatmap: separate from the TE heatmap; placement undecided.
4. Running log with every check and caveat, as in `CAP_BINDING_DATA.md`.

## Open questions (asked, not answered)

1. Lee / De Silva outcome: bound yes/no among expressed genes, reads per mRNA among bound genes,
   or both.
2. Jurkat columns: A/C/B + B/C + D/L in both states (6), D/L only (2), all (10), or activated only.
3. NPC: whole transcript + 3'UTR end, full region split, or whole transcript only.
4. Heatmap: one eIF3-binding heatmap with the subunit-seq columns, a CLIP-only one, or both.
