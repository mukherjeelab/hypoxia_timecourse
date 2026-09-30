# eIF3d cap-binding (subunit-seq) and PAR-CLIP data: what we have and what is wrong with it

Recorded 2026-09-30. **Status:** provisional models from the deposited RPKM are running into a
**separate** heatmap (section 9). Nothing cap-binding enters the combined TE heatmap. Reprocessing
from raw reads is deferred - every result below carries that caveat.

**Summary.** Complete all-gene data for eIF3d cap binding exist on GEO and reproduce both papers'
published tables, so a continuous cap-binding outcome is possible. But raw enrichment (IP over input)
rises steeply with transcript length (Spearman 0.58-0.78), neither paper sequenced a background IP
(IgG or untagged) that could show whether that is biology or pull-down background, and our own
RIP-seq shows IP-vs-RNA picks up the same length signal while IP-vs-IgG removes it. **Open decision:**
model raw enrichment (length-confounded), the stress-vs-baseline ratio (cancels shared capture bias),
or reprocess the raw reads first to get real DESeq2 output for every gene and read-coverage shapes.

Terms used below:

- **Subunit-seq** - the Lee lab eIF3d-TEV cap-binding assay (not PAR-CLIP). Described in full below.
- **Enrichment** - IP over input for one gene; how much of that mRNA survived the pull-down relative
  to the typical mRNA.
- **Input** - 10% of the cytoplasmic lysate, i.e. total cytoplasmic RNA. It is not a control IP.
- **Background IP** - a pull-down that should capture nothing specific (IgG, untagged cells). It
  measures what sticks regardless of the bait. Neither study has one.
- **Stress-vs-baseline ratio** - enrichment in the stressed condition divided by enrichment in its
  own baseline (glucose deprivation / complete media; thapsigargin / DMSO). Any capture bias shared
  by both conditions cancels.
- **RPKM** - reads per kilobase per million: a gene's reads / (gene length in kb x library size in
  millions). Makes genes of different length comparable *within* a library.

---

## 1. The data

All from HEK293T. Gene-level only - no transcript IDs anywhere, so the feature matrix's MDA-MB-231
isoform stands in for each gene.

| Study | Conditions | Published table (hits only) | All-gene data |
|---|---|---|---|
| Lamper et al. 2020, *Science* 370:853 | complete media (CM), glucose deprivation (GLU) | S1 = GLU hits (DESeq2, 668; **not on disk**), S2 = CM hits (DESeq2, 664; `~/Downloads/abb0993-lamper-sm-table-s2 (4).xlsx`), S3 = IP/input ratios + DESeq2 log2FC/padj for the 1,141-gene union (`mpra_preliminary/data/lamper-s3.xlsx`) | GEO **GSE158249**: `GSE158249_Complete_RPKM.xlsx`, `GSE158249_Glucose_RPKM.xlsx` |
| Mukhopadhyay, Amodeo, Lee 2023, *Mol Cell* 83(18) | DMSO, thapsigargin (Tg, chronic ER stress, 600 nM 16 h) | S2 = Tg hits (DESeq2, 563; `mpra_preliminary/data/integrated_stress_response.xlsx`) | GEO **GSE236188**: `GSE236188_SubunitSeqRep{1,2}_rpkm.txt.gz` |
| Lee et al. 2015 PAR-CLIP | eIF3 subunits a, b, d, g | cluster calls only (`mpra_preliminary/data/lee_parclip.csv`, hg19) | **none** - bound yes/no only |

The GEO files were downloaded to a session scratchpad, not into the repo. Re-fetch from
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158249/suppl/` and
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236188/suppl/`.

Each GEO file holds IP and input RPKM, two biological replicates, for ~58,000 genes.
Raw reads: Lamper SRA SRP284032; Mukhopadhyay under GSE236188.

### Do we have DESeq2 for IP over input?

**Only for the published hits.** For every other gene - and for DMSO entirely - we have per-sample
RPKM, from which IP/input can be computed but with no padj, no standard error and no count-based
noise model. RPKM cannot be converted back to counts (library sizes and the authors' gene lengths
are not deposited). Real all-gene DESeq2 needs the raw reads re-aligned and counted.

---

## 2. Verified against the papers

| Check | Result |
|---|---|
| Lamper S3 per-replicate columns = IP RPKM / input RPKM from GEO | every numeric column reproduced to ~1e-14, all 1,141 genes. The S3 columns labelled `RPKM_*_rep*` are really IP/input **ratios** |
| Lamper S3 `IP/IN (AVG)` and `Fold Enrichment (GLU/CM)` | reproduced to ~1e-14 |
| Lamper DESeq2 log2FC vs our log2(IP/input) | GLU slope 1.002, r 1.000; CM slope 0.989-1.034, r 0.983-0.987; constant offset -1.5 to -2.1 (library-size normalisation) |
| Lamper hit counts at padj < 0.005, log2FC > 2 | 668 GLU (paper text: "668 transcripts"), 664 CM (= Table S2 exactly, set-equal) |
| CM and GLU overlap | 191 genes - matches S2's own `eIF3d_glucose_dep_target` column (191) |
| Lamper "~70% of GLU targets >= 2-fold GLU/CM" | 81% from RPKM ratios, 64% from DESeq2 log2FC difference - the paper's figure sits between; cannot tell which they used |
| Mukhopadhyay S2 contrast | IP_Tg / Input_Tg: Pearson 0.997 (other contrasts 0.03-0.73); slope 0.906, r 0.950 against their DESeq2 |
| Mukhopadhyay named targets | ALKBH5 and EIF2AK4 (GCN2) present in S2 |
| Mukhopadhyay hit count | S2 has 563; text says 564 |
| Two copies of Lamper S3 (repo vs `~/Downloads`) | identical values (the Downloads copy has an extra header row) |

---

## 3. How subunit-seq works (Lamper 2020 supplementary methods)

1. **Cell line.** HEK293T with an extra CMV-driven copy of eIF3d knocked into the AAVS1 safe harbor.
   A TEV protease site sits inside the protein just before the cap-binding domain; an HA tag is on
   the C-terminus, the cap-binding-domain side. Endogenous eIF3d is untouched.
2. **UV 254 nm crosslink** in living cells, 2.5 min - covalent protein-RNA bonds at direct contacts.
3. **Lyse** (NP-40), spin out nuclei; **10% of the cytoplasmic lysate = input**.
4. **Anti-HA IP**: pulls down the whole eIF3 complex and whatever it sits on. High-salt washes.
5. **TEV cleavage on the beads, then wash.** The rest of eIF3 is washed away. **The cap-binding
   domain stays on the beads** with any RNA crosslinked to it - at the cap or anywhere else.
6. **On-bead cleanup** to destroy RNA whose cap was not shielded:
   Antarctic phosphatase (uncapped ends -> 5'-OH) -> Cap-Clip (removes any reachable cap) ->
   T4 PNK (5'-monophosphate on everything exposed) -> XRN1, 2 h (5'->3' digestion of everything
   with a 5'-monophosphate). RNA whose cap sat in eIF3d's pocket survives intact.
7. Proteinase K, phenol-chloroform; **libraries from poly(A) RNA** - reads span whole mRNAs,
   so no 5'-end pile-up is expected.
8. Counts per gene (summarizeOverlaps) -> filter (Lamper: < 50 input reads; Mukhopadhyay: < 60
   reads) -> **DESeq2 IP vs input** -> hits (Lamper padj < 0.005 & log2FC > 2; Mukhopadhyay padj <
   0.1 & FC > 2).

**Controls.** Sequenced: IP and input only. **No IgG, no untagged parental cells, no no-TEV IP.**
Specificity was shown only by qPCR on a few genes (Jun, Raptor, ALKBH5, GCN2): m7GpppG but not GpppG
competes binding away, and PSMB6 serves as a non-bound mRNA.

**Where non-cap signal could get through** (hypotheses, untested):

- incomplete Cap-Clip / XRN1 digestion - a partly-chewed long mRNA keeps its poly(A) tail and is
  still sequenced;
- internal crosslinks - an mRNA crosslinked to the domain away from its cap stays on the beads; its
  cap is removed, but XRN1 may stall at the crosslink, leaving a poly(A) 3' piece;
- caps shielded by something other than eIF3d;
- nonspecific sticking to beads or domain.

The first two would make enrichment rise with length.

---

## 4. Our processing vs theirs

| Step | Papers | Ours | Same? |
|---|---|---|---|
| Length normalisation | **none** - DESeq2 on raw counts; size factors correct depth per sample, not length per gene | RPKM divides by length, but it **cancels** in IP / input for the same gene | yes - neither can create or remove a length effect |
| Enrichment value | DESeq2 log2FC | log2(IP RPKM / input RPKM), mean of 2 replicates | yes up to a constant (slope ~1) |
| Pseudocount | none (count model) | +0.1 RPKM | negligible - length rho 0.781 without vs 0.780 with (CM) |
| **Expression filter** | **>= 50 / 60 input reads** | **input RPKM >= 1** in both replicates | **no** |
| Hit thresholds | as above | same | yes |

**The filter is the one place length enters.** At equal expression a long gene has more reads, so a
read-count filter keeps long, low-expression genes that an RPKM filter drops. 73% of the published
GLU hits have input RPKM < 1 (median 0.59), so our filter removes most of them. Published GLU hits
have median transcript length 4,415 nt vs 2,559 nt for other matrix genes.

---

## 5. The length problem

**Why RPKM does not fix it.** RPKM corrects "a long mRNA gives more reads per molecule" - identical in
IP and input, so it cancels in the ratio regardless. What we see is "more *molecules* of a long mRNA
were captured", which no length normalisation can undo.

**Size of the effect** (log2 IP/input, input RPKM >= 1 both replicates, lengths from the MDA-MB-231
isoform):

| Condition | n | rho with log2 total length | rho with input RPKM | replicate agreement |
|---|---|---|---|---|
| CM | 3,925 | +0.78 | -0.43 | 0.87 |
| GLU | 3,246 | +0.71 | -0.54 | 0.98 |
| Tg | 5,898 | +0.58 | -0.51 | 0.56 |
| DMSO | 4,347 | +0.70 | -0.44 | 0.51 |

Conditions agree with each other at 0.80-0.88, including CM vs DMSO (two unstressed HEK293T
experiments three years apart) - most of the enrichment is a per-gene property shared across
conditions. GC correlates at <= 0.16.

**Not specific to one region.** Partial Spearman, each region's length with the other two held:

| Condition | 5'UTR | CDS | 3'UTR | R2, three regions | R2, total length alone |
|---|---|---|---|---|---|
| CM | +0.42 | +0.60 | +0.54 | 0.63 | 0.61 |
| GLU | +0.39 | +0.52 | +0.54 | 0.55 | 0.52 |
| Tg | +0.27 | +0.33 | +0.45 | 0.35 | 0.34 |
| DMSO | +0.22 | +0.55 | +0.51 | 0.47 | 0.48 |

All three regions contribute; total length explains as much as the three together; the **5'UTR -
the cap-proximal region - has the smallest unique share in every condition** (3-8% of variance vs
8-22% for CDS and 15-19% for 3'UTR). Holding input expression and GC barely changes this. Slope of
log2 enrichment on log2 total length: 1.40 CM, 1.20 GLU, 0.66 Tg, 0.89 DMSO.

**Our own RIP-seq says a background IP removes it.** Spearman with log2 total length, ~9,850 genes:

| Bait | IP vs RNA | IP vs IgG |
|---|---|---|
| eIF3d (normoxia 1 h; hypoxia 1 / 4 / 24 h) | +0.31, +0.30, +0.47, +0.15 | -0.15, -0.21, -0.31, -0.34 |
| eIF4E (same timepoints) | +0.28, +0.30, +0.49, +0.24 | -0.15, -0.16, -0.31, -0.25 |

Both baits favour long mRNAs against RNA - eIF4E too, which has no reason to prefer long 3'UTRs - and
the sign flips against IgG. Caveats: a different assay (no crosslink, no TEV, no enzymatic cleanup,
MDA-MB-231), and the RIP-seq feature notebook notes its IP-vs-RNA contrast is confounded with
sequencing batch. Subunit-seq's correlations (0.58-0.78) are *stronger* than RIP's IP-vs-RNA ones,
which suggests the cleanup does not remove it - suggestive across assays, not proof.

---

## 6. Other problems found

**Feature-matrix bug (unfixed).** `clip_lee_glu_dep_bound` holds the **complete-media** set: all 480
flagged genes are CM-bound, only 154 are GLU-bound. The `_lfc` columns are correct (r = 1.0 with
source). So `clip_cap_binding_count_v2` = ISR + CM, and glucose-deprivation binding has never entered
any model. The earlier audit saw the CM and GLU bound columns were identical and dropped one, without
noticing the kept one was the wrong condition. Affects only the full-feature heatmap (0.05% SHAP
share); the intrinsic-only runs exclude every `clip_*` column.

**Existing PAR-CLIP flag is lossy.** `clip_lee_parclip_d_only_bound` misses 42 of the 285 eIF3d-bound
genes in the matrix (19 with 5'UTR clusters). Rebuild from `lee_parclip.csv`, not the matrix.

**`clip_eif3_count`** (3.45% SHAP share, rank 11 in the full model) mixes Lee PAR-CLIP with the Cate
NPC CLIP, so the Lee data is diluted inside it.

---

## 7. Plan as agreed so far

- Cap binding and PAR-CLIP become **outcome columns** (what predicts eIF3d binding), not feature rows.
- **Intrinsic features only**, to avoid circularity with the `clip_*` features.
- Cap binding: three conditions - CM, GLU, Tg. PAR-CLIP: subunits a, b, d and g as separate columns,
  mature-mRNA regions only (5'UTR, start codon, CDS, 3'UTR; no intron / lincRNA / miRNA / rRNA
  clusters). PAR-CLIP is bound yes/no.
- Columns go into both the combined all-cell-line intrinsic heatmap and a new binding-only heatmap.

## 8. Open questions

1. **Which cap-binding quantity to model**:
   (a) raw enrichment, reading its length share as at least partly technical;
   (b) stress-vs-baseline ratio - cancels shared capture bias, is what Lamper used for their stress
   claim (fig. S4F), but loses the unstressed-binding column;
   (c) reprocess raw reads first.
2. **Reprocess from SRA?** Download, align, count, then DESeq2 with the papers' filters. Gives
   all-gene log2FC / padj / lfcSE computed their way, and the IP vs input coverage shape. A genuine
   cap-protected mRNA should have IP coverage shaped like input; internal-crosslink survivors would
   be 3'-shifted. Several hours of compute.
3. **Expression filter** - match the papers' read-count filter (keeps long low-expression genes) or
   keep RPKM >= 1.
4. **PAR-CLIP background** - "unbound" = every matrix gene without a cluster, or only genes expressed
   in HEK293T (subunit-seq input RPKM gives a same-cell-line expression measure).
5. **Fix the `clip_lee_glu_dep_bound` bug** - add a corrected variant alongside, per the variant-pair
   convention, and decide whether the full-feature heatmap columns are re-run.

---

## 9. Run log: provisional models from deposited RPKM (started 2026-09-30)

**Decision (Kate, 2026-09-30):** do not reprocess raw reads now. Model the RPKM-derived outcomes,
keep them in their own heatmap, and record that **all of it must be rerun from the raw reads**
(SRA SRP284032; GSE236188) with DESeq2 on counts. DESeq2 stress-over-unstressed does not exist
anywhere - neither paper ran it - so the combined heatmap gets no cap-binding column until that
reprocessing is done.

### What was built

| Piece | File |
|---|---|
| Source files, copied into the repo (md5-checked against the GEO downloads) | `accessories/subunit_seq/` - both GEO RPKM sets, Lamper S2 + S3, Mukhopadhyay S2 |
| Outcome prep: replication checks, five outcomes, hit / non-hit gene sets | `10s_subunit_seq_outcome_prep.Rmd` -> `output/subunit_seq_hek293t_capbind_<k>_rpkm.csv` |
| Regression binding mode | `10_rf_regression.Rmd`: `outcome_kind = "binding"` (requires `cell_line = "HEK293T"`), `min_spot_sep`, `heldout_sep_check`. Defaults leave every TE run unchanged |
| Jobs (intrinsic features only, 57) | `jobs_subunit_seq_rpkm.txt` |
| Separate heatmap | `15c_render_subunit_seq.R` -> `15_shap_heatmap_subunit_seq_rpkm.html`; the MDA-MB-231 3d hypoxia 1hr intrinsic TE column is included as a reference |

Outcome definitions: enrichment = mean over two replicates of log2((IP + 0.1) / (input + 0.1)),
RPKM units; a gene is kept only if input RPKM >= 1 in every input sample the outcome uses. The
stress ratios are the stressed enrichment minus the baseline enrichment (replicate means).
Pseudoautosomal genes appear as X and `_PAR_Y` copies in the Mukhopadhyay files (reads split
between them; ZBED1 is a published hit only under its `_PAR_Y` id) - the copies are summed. Three or four
duplicated ids in the Lamper files (none a hit) are dropped.

### Checks passed (all asserted)

- Lamper S3 IP/input ratios reproduced from GEO, max relative deviation ~1e-14.
- GLU hits = 668 (paper text); CM hits set-equal to Table S2 (664); overlap 191 = Table S2's own flag.
- Mukhopadhyay S2 = IP_Tg / Input_Tg (r 0.950, against 0.757 for the DMSO ratio); all 563 found in GEO.
- ALKBH5: thapsigargin enrichment +2.71, 97th percentile.
- `10_`: permuted-outcome null ~0 (Spearman +0.0001 to +0.023); 57 features, identical to the
  MDA-MB-231 intrinsic reference (parity asserted).

### The five outcomes

| Outcome | Genes in matrix | Published hits kept | Hits minus non-hits | rho with length | Test-render held-out Spearman |
|---|---|---|---|---|---|
| complete media, raw | 3,942 | 141 of 664 | 3.25 | +0.78 | 0.83 |
| glucose deprivation, raw | 3,261 | 145 of 668 | 2.96 | +0.71 | |
| thapsigargin, raw | 5,942 | 361 of 563 | 1.53 | +0.58 | |
| glucose / complete media | 2,953 | 140 of 668 | 0.53 | **-0.23** | |
| thapsigargin / DMSO | 4,287 | 160 of 563 | **0.11** | **-0.43** | 0.56 |

(Held-out Spearman is filled in from the test renders so far; the full runs are in progress.)

### New caveats found while building

- **The stress ratios do not cancel length - they reverse it.** Stressed enrichment rises less
  steeply with length than baseline enrichment (slopes: Tg 0.66 vs DMSO 0.89), so the ratio
  carries a *negative* length correlation. Neither ratio is length-free.
- **Most published thapsigargin targets are not stress-induced.** On Tg / DMSO they sit only 0.11
  above non-hits; they are nearly as enriched in DMSO. On held-out genes the model predicts them
  *below* non-hits (-0.18 vs -0.14), because they are long and the ratio falls with length. The
  held-out separation check is therefore reported, not asserted, for both ratio outcomes.
- **The RPKM filter keeps only 21-22% of Lamper's published hits (141 of 664 CM, 145 of 668 GLU, counted within the matrix)** (they are long,
  low-expression genes the papers' read-count filter admits). The Tg set keeps 64%.
- Raw complete-media enrichment is predicted from sequence alone at Spearman 0.83 - higher than
  any TE column - which is what a length-dominated outcome would give.
- Gene level, HEK293T, MDA-MB-231 isoforms; no padj; the gene sets written for `10_` carry `_TE_`
  in their filenames only because `10_` builds the path that way.

### Progress

- [x] outcome prep knitted, all checks green
- [x] binding mode test-rendered on complete media (raw) and thapsigargin / DMSO (ratio)
- [x] five regression + SHAP runs
- [x] separate heatmap rendered: `15_shap_heatmap_subunit_seq_rpkm.html` (via `15c_render_subunit_seq.R`)
- [x] results recorded here

### Permutation check changed for binding runs

Glucose deprivation first failed `10_`'s single-shuffle bound (|shuffled Spearman| < 0.1). Ten
shuffled forests per column scattered -0.21 to +0.18 around zero with random sign (glucose mean
|0.10|): small test sets (590-1,190 genes) plus a strongly feature-determined outcome let a
meaningless forest correlate by chance. Not leakage. Binding runs now use `perm_n = 20`: the
shuffles must centre within 0.1 of zero and the real model must beat their maximum. TE runs keep
`perm_n = 1`, unchanged. All five columns pass (`10s_rerun_binding_permcheck.sh`).

### Results (intrinsic features)

| Outcome | Held-out Spearman | Top SHAP families (share %) |
|---|---|---|
| complete media, raw | 0.833 | **length 56.4**, g4 12.0, structure_utr3 6.3 |
| glucose deprivation, raw | 0.787 | **length 51.4**, g4 15.5, gc 8.1 |
| thapsigargin, raw | 0.649 | **length 47.9**, g4 12.7, structure_utr3 9.9 |
| glucose / complete media | 0.633 | **gc 33.0**, length 20.9, structure_cds 10.2 |
| thapsigargin / DMSO | 0.561 | **length 48.9**, codon_optimality 12.2, structure_cds 7.5 |
| *reference: MDA-MB-231 3d hypoxia 1hr TE* | *0.471* | *gc 33.0, length 19.5, codon_optimality 11.9* |

- **Raw enrichment is a length readout** (48-56% of SHAP), as the rho 0.58-0.78 predicted.
- **The glucose / complete-media ratio has almost exactly the TE reference's family profile**
  (gc 33.0 vs 33.0, length 20.9 vs 19.5). Suggestive, not yet a finding: no noise floor, RPKM
  only, and the ratio still carries a negative length correlation (-0.23).
- **Thapsigargin / DMSO stays length-dominated**, with the opposite sign (rho -0.43).
