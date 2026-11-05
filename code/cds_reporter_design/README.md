# CDS reporter design

Where a candidate reporter CDS falls in the MDA-MB-231 coding landscape, as the
baseline for choosing or recoding a reporter.

Nothing here writes to `output/predictive_modeling/`. Outputs go to
`output/cds_reporter_design/`; figures to `plots/` with a `cds_reporter_` prefix.

## Notebooks

| Notebook | Role |
|---|---|
| `reporter_helpers.R` | shared, cell-line agnostic machinery: reporter sequence import and verification, CDS metrics, the density plotter, percentiles |
| `01_cds_length_gc3_distributions.Rmd` | MDA-MB-231: CDS length, GC% and GC3% densities for the expressed transcriptome and the eIF3d hypoxia-promoted set, with NanoLuc reference lines and percentiles |
| `02_cds_length_gc3_distributions_mcf7six1.Rmd` | the MCF7-SIX1 counterpart; also writes the two-cell-line medians table. Run `01_` first, since section 8 reads its exported medians |
| `03_cds_length_gc3_fusion_constructs.Rmd` | MDA-MB-231 only: the eGFP / mCherry / TagRFP -P2A-NanoLuc fusions assembled from parts, with NanoLuc alone as the reference |
| `04_cds_length_gc3_flag_constructs.Rmd` | MDA-MB-231 only: 1x / 2x / 3x FLAG-P2A-NanoLuc, tags cut nested from the supplied `flag3x` |
| `05_cds_length_flag_copy_ladder.Rmd` | MDA-MB-231 only: a FLAG copy-number ladder as a length series; the step size is the `flag_step` parameter and the output filenames carry it |
| `06_cds_length_smbit_p2a_ha.Rmd` | MDA-MB-231 only: SmBiT-P2A-1xHA (120 nt), -1xFLAG (117 nt) and -FLAG-eGFP (840 nt), no NanoLuc; three cumulative views of each panel, and the length panel is also drawn zoomed because ~120 nt is unreadable on the standard 6 kb axis |
| `07_cds_length_egfp_copy_ladder.Rmd` | MDA-MB-231 only: an eGFP copy ladder on `06_`'s SmBiT-P2A-FLAG backbone, Waldo linker included (0/1/2/3 copies = 117/876/1,599/2,322 nt by default); `max_gfp_copies` and `linker_between_copies` are parameters and the output filenames carry both |
| `08_gc_matched_recode_egfp_ladder.Rmd` | the same four rungs as `07_`, with the FLAG tag and every eGFP copy **recoded synonymously** onto the 231 composition: GC3 93% to 59%, GC 61% to 50%, protein unchanged. Copies are recoded independently, which also removes `07_`'s tandem repeat |

**The FLAG ladder confounds length with GC.** It is a clean length series on the
GC3 axis - GC3 moves 1.6 points across 1x to 51x, 5% of the transcriptome IQR,
against a length change of 90% of its IQR. But the FLAG epitope is GC-poor (40.7%
vs NanoLuc's 53.9%), so total GC falls 7.6 points over the same span, 53% of the
GC IQR: nearly as large a move, in distribution units, as the length change the
ladder exists to produce. Use it to vary length against GC3 and wobble effects;
do not use it to vary length against total GC. `05_` asserts the GC3 property and
prints a warning with the GC drift rather than failing on it.

**The bare-tag SmBiT constructs are off the bottom of the length distribution;
adding eGFP fixes length and ruins GC3.** This is the sharpest version of the
fluorophore result below, because here the same backbone is measured with and
without eGFP:

| construct | nt | pct expressed | pct promoted | GC3 | pct expressed | pct promoted |
|---|---|---|---|---|---|---|
| SmBiT-P2A-HA | 120 | 0.3 | 0.5 | 67.5 | 66.9 | 82.1 |
| SmBiT-P2A-FLAG | 117 | 0.3 | 0.4 | 71.8 | 74.9 | 87.5 |
| SmBiT-P2A-FLAG-link-eGFP | 876 | 31.2 | **49.1** | 91.1 | **99.2** | **99.8** |

FLAG-link-eGFP is the only construct in this folder that lands essentially **on**
the promoted set's own length median (876 nt against 891, the 49th percentile of
promoted). It buys that by adding 723 nt of eGFP at GC3 96.7%, which puts the
assembled construct at the **99.1st percentile of the expressed transcriptome on
GC3** - further into the tail than any other reporter here. The Waldo linker, at
GC3 41.7%, pulls it back 2.1 points from where it sat without one, which is real
but nowhere near enough. Against a 231 target set that is short *and* GC3-poor,
the bare tags are wrong on length and the eGFP version is wrong on GC3; nothing
in this family is right on both. A codon-shuffled eGFP would be the way to keep
the length and lose the GC3.

The two bare-tag constructs are not points inside the distribution at all, so the
percentile is the number to quote and their reference lines are unreadable on the
wide panel - which is why `06_` draws the length panel twice, once at 6 kb for
comparability and once zoomed to 2 kb.

**The HA and FLAG versions are one point for this analysis.** They differ by 3 nt,
0.2% of the transcriptome's length IQR, so the two reference lines are
superimposed on the length panels. GC3 separates them most (4.3 points, 15% of
the IQR) because the tags differ more in composition than in length: HA is GC3
44.4% over 9 codons, FLAG 62.5% over 8. `06_` section 6 reports the spread for
the bare pair and for all three separately, which is the only place it is
readable.

**The GFP chromophore sits at residues 66-68, not 65-67.** avGFP numbers it
T65-Y66-G67, but eGFP carries an extra Val at position 2, so everything shifts by
one. `06_` also asserts the `LTYGVQCF` context rather than the bare `TYG`,
because `TYG` occurs a second time at residue 39 and a motif-presence check would
pass on a sequence with no chromophore in place.

**SmBiT is not the NanoLuc C-terminus, and `06_` pins the difference.** SmBiT
(`VTGYRLFEEIL`) aligns against NanoLuc's C-terminal 11-mer (`VTGWRLCERIL`) at 8
of 11 residues - W4Y, C7F, R9E - so it is recognisably derived from it but is not
a substring of it, and its codons are different too. A check that only confirmed
similarity would pass if the split were accidentally reverted, so `06_` asserts
exactly 3 mismatches and asserts the peptide is *not* found in the NanoLuc
protein. The spreadsheet holds two versions, `smBit n` (36 nt, with the initiator
ATG) and `smbit c` (33 nt, no ATG, for C-terminal fusion); they differ by one
synonymous base in the final Leu codon (CTC vs CTG) as well as the ATG.

**P2A puts the tag on the downstream product, and only the eGFP version has mass
to run on a gel.** Skipping happens between the final Gly and Pro of P2A, so
every `06_` construct yields the *identical* 30-aa upstream product (SmBiT + 18
aa of P2A remnant - identical because the tag is entirely downstream of the skip)
and a downstream product of Pro + tag: **10 aa for HA, 9 for FLAG, 250 aa
(~27.5 kDa) for FLAG-eGFP**. A free 9-10 aa peptide has no mass to resolve on a
gel, so the bare tags work as a concept for "the second cistron was translated"
but not as western targets. `06_` reports every product rather than leaving the
reader to work out where the tag ended up, and prints which downstream products
clear 50 aa.

**The Waldo linker sits between the FLAG tag and eGFP**, in `06_` and `07_` both.
It is 36 nt / `GSAGSAAGSGEF`, GC3 41.7%, and it **ends in an in-frame `GAATTC`
EcoRI site** (Glu-Phe) - that is how folding-reporter constructs are built, but
it leaves one internal EcoRI site per linker in the ORF. Both notebooks count
`GAATTC` per construct and assert the count equals the linker count, with the
linker-free rungs asserted to carry none so the check is discriminative.

**The linker precedes every eGFP copy**, so `07_`'s repeat unit is linker + eGFP
= 759 nt. `linker_between_copies` is kept as a parameter (`TRUE` default) because
`FALSE` - one linker at the FLAG-eGFP junction only, eGFP copies fused directly,
723 nt per rung - is a reasonable reading of the same design and its outputs are
still reproducible under their own filename tag. `TRUE` is the default because
directly fused tandem fluorophores are the folding problem the Waldo linker was
designed to relieve, and because the extra linker dilutes GC3 further (3x: 93.0%
against 94.6%). The 0-copy rung never carries a linker under either setting.

**The eGFP ladder is the cleanest length series in this folder, and it runs
along the far GC3 tail for its whole range.** `07_` stacks exact tandem copies of
the linker + eGFP unit on `06_`'s SmBiT-P2A-FLAG backbone:

| construct | copies | nt | pct expressed | pct promoted | GC | GC3 | GC3 pct expr |
|---|---|---|---|---|---|---|---|
| SmBiT-P2A-FLAG | 0 | 117 | 0.3 | 0.3 | 54.7 | 71.1 | 73.7 |
| -link-1xeGFP | 1 | 876 | 31.2 | 49.1 | 61.0 | 91.1 | 99.2 |
| -link-2xeGFP | 2 | 1,635 | 63.7 | 76.8 | 61.4 | 92.5 | 99.6 |
| -link-3xeGFP | 3 | 2,394 | 81.2 | 87.9 | 61.6 | 93.0 | 99.7 |

Across the three eGFP rungs **length moves 114% of the expressed IQR while GC
moves 4.3% and GC3 6.6%** - length dominates composition by 17x, where `05_`'s
FLAG ladder managed only ~2x on GC. That is because the repeat unit is a whole
protein rather than a 24 nt epitope, so above one copy the construct is
asymptotically just eGFP and its composition stops moving. The ladder spans the
31st to 81st percentile of the expressed transcriptome, straddling both medians,
and the 1x rung sits essentially **on** the promoted median (49.3rd percentile of
promoted).

**But the level, not the drift, is the problem**: every rung sits above the
99.7th percentile of expressed GC3 and at or above the 99.9th of promoted. The
ladder holds composition still at a value no mRNA in the comparison set has. Use
it to vary length when GC3 is explicitly not the variable; do not use it to argue
anything about a typical mRNA's length response. A **codon-shuffled eGFP** - the
recoding machinery `64_` used for NanoLuc - would keep every length and move GC3
to any target, and `07_` is the baseline such a design would be compared against.

**Tandem identical repeats are a construction risk, so `07_` measures them.** The
2x and 3x constructs carry 723 and 1,446 nt of exactly repeated sequence, which
recombines in *E. coli* and is hard to sequence through. The 1x construct is
asserted to carry **no** internal repeat, so the measure is discriminative rather
than reporting a number that would look plausible either way. Codon-divergent
copies fix the repeat and the GC3 problem at the same time.

**`05_` uses a different repeat unit from `04_`.** Beyond three copies there is
nothing left to nest, so `05_` tandemly repeats the canonical DYKDDDDK epitope
(24 nt, exactly linear), while `04_`'s 3x uses the three distinct repeats of the
supplied `flag3x` (66 nt). A 3x differs by 6 nt between the two schemes. Both take
their codons from the spreadsheet.

**FLAG is composition-neutral where a fluorophore is not.** The three FLAG
fusions span 720-762 nt at GC3 64-65%, against NanoLuc alone at 636 nt / 65.1% -
so tagging with FLAG leaves the construct essentially where NanoLuc already sat,
while any of the fluorophores moves GC3 by 9-16 points. The three FLAG constructs
also barely differ from *each other*: their whole range is 3-5% of the
transcriptome's own IQR on every axis, so they overlap on the panels and
`04_` section 6 reports the separation numerically instead.

**A sequence can be in frame, the right length, and still be backwards.** The
`flag3x` row was first pasted as the reverse complement: 66 nt, divisible by 3,
and it passed every length and composition check while translating to nonsense.
`04_` therefore asserts the *peptide*, cross-checks it against the 3xFLAG embedded
in `gfp3xflagdestab`, and asserts that the reverse complement does **not** encode
3xFLAG, so the orientation check is discriminative rather than vacuous. The two
copies differ by one synonymous base, so they are compared as peptides, not DNA.

**The fluorophore dominates a P2A fusion's composition.** eGFP and mCherry are
near-maximally GC3-rich as supplied (96.7% and 97.0%), TagRFP less so (83.1%), and
each is larger than NanoLuc itself. A fusion therefore sits far above where
NanoLuc alone sits: GC3 65.1% alone, 81.1 / 81.2 / 74.1% fused. The consequence
for codon design is that **the fluorophore can cancel the recoding**: the
3d-adapted NanoLuc is GC3 12.7% alone, but 57.6% behind eGFP, which is the
transcriptome median. Behind TagRFP it is 50.4%, essentially the eIF3d-promoted
median. Always evaluate a design on the assembled construct.

**P2A is not a stop.** The ribosome skips and the proteins separate
co-translationally, so the mRNA carries one ~1.4 kb ORF, and it is measured as one
CDS.

**The GC3 result is opposite in the two cell lines**, and the figures should never
be read as if it were not. In MDA-MB-231 the eIF3d hypoxia-promoted set is
GC3-*poor* relative to its transcriptome (50.0% vs 57.3%); in MCF7-SIX1 it is
GC3-*rich* (67.8% vs 59.0%). Both are clean monotone gradients against their own
eIF3d-inhibited sets (231: 50.3 promoted vs 71.4 inhibited; MCF7-SIX1: 67.8 vs
42.2), so neither is noise. This is not a sign error in the MCF7-SIX1 pipeline:
see the three rulers in `02_` section 2, and note that the negated MCF7-SIX1
contrast correlates *positively* with the matched 231 contrast, so the two
orientations agree while the GC3 association does not. A reporter tuned for one
cell line is therefore tuned against the other.

## Inputs

| Input | Used for |
|---|---|
| `accessories/cds seq.xlsx` | the NanoLuc construct sequence, imported not transcribed |
| `output/translation_categories_si3d_vs_sictrl_hypoxia_1hr.csv` | the eIF3d hypoxia-promoted set (`te_lfc > 0.5`, `te_padj < 0.05`, n = 1,604 genes) |
| `output/predictive_modeling/feature_matrix_external_stability.rds` | the expressed transcriptome pool and its one-isoform-per-gene transcript choice |
| `2025_MDAcells_riboseq/salmon/NM2025_009{4,5,6}` | hypoxia siCTRL fallback isoform for promoted genes absent from the matrix |
| `accessories/human/txdb.gencode49.sqlite` | CDS ranges |
| `output/predictive_modeling/nanoluc_codon_designed_sequences.csv` | the four designed NanoLuc variants from notebook 64 |
| `output/translation_categories_si3d_vs_sictrl_mcf7six1_hypoxia_1hr.csv` | MCF7-SIX1 promotes-positive TE, written by `10m_mcf7six1_outcome_prep.Rmd` |
| `output/genesets/mcf7six1_hypoxia_3d_promotes_TE_1hr_lfc0.5.csv` | the committed MCF7-SIX1 set (853 genes) that `02_` must reproduce |
| `output/translational_control_overview/feature_matrix_tco_mcf7six1.rds` | MCF7-SIX1 expressed transcriptome, own dominant-isoform selection |

## Conventions specific to this folder

**The reporter sequence is imported and cross-checked, never retyped.** `01_` reads
the construct from `accessories/cds seq.xlsx` and asserts it is identical to the
string `code/predictive_modeling/62_cds_length_gc_hypoxia_promotes.Rmd` hard-codes,
by parsing notebook 62's own source. If either changes, the check fires rather than
the two drifting apart silently. Any new reporter added here should get the same
treatment.

**Every legend entry carries the number it stands for.** `plot_density_luc()`
appends each density curve's group median and each reference line's own value to
its legend label, so a panel can be read without the percentile table beside it.
Units follow the **metric, not the axis**: the length axis is drawn in kb because
kb ticks are readable, but a construct is specified in nt, so the labels say nt.
This is in the shared helper, so all six notebooks carry it and all their figures
were re-knit together when it was added.

**All CDS metrics are computed on stop-codon-free ORFs**, genes and reporter alike.
GENCODE txdb CDS ranges include the terminal stop, so it is stripped from the gene
sequences first. The reporter's own TAA is stripped the same way.

**The comparator is the expressed transcriptome, not the screened negative control
set.** The question is where a reporter sits relative to a typical 231 mRNA. The
negative control set is genes selected for having no eIF3d effect, which makes it a
contrast rather than a background; notebook 62 carries that version. The grey curve
includes the promoted genes, since it describes the whole population, but the
Wilcoxon compares promoted against expressed-minus-promoted so the two samples are
disjoint.

**The C-terminal extension is 41 aa, not 40.** All three lab reporters carry the
mouse ODC PEST degron, which is 40 aa, but one asparagine sits between the ORF and
the degron, encoded by the `AAT` inside the `GAATTC` EcoRI site at the cloning
junction (`...CTGGCG|AATTCTCACGGC...`). Splitting a construct on the 40-mer degron
puts the boundary 3 nt late and leaves the junction codon on the reporter side.
`01_` asserts both the 40-mer degron and the 41-mer extension.

**MCF7-SIX1 directionality is opposite and its usual check silently skips.** The
MCF7-SIX1 DESeq2 interaction coefficient is log2(KD / sictrl), so positive means
eIF3d *inhibits* — the reverse of the 231 files. `02_` reads the already-negated
promotes-positive tables from `10m_mcf7six1_outcome_prep.Rmd` so the sign is
applied exactly once in the repo. The standing project check, "the knockdown
target's own `rna_lfc` must be strongly positive", cannot run: the MCF7-SIX1
source has no `rna_lfc` column, so that check *skips* rather than fails. `02_`
replaces it with three rulers that can all fail — the table equals the negation of
its source column, EIF3D's own `te_lfc` is strongly negative (−1.79) in the
promotes-positive orientation, and the derived set reproduces the committed
853-gene set exactly.

**MCF7-SIX1 uses its own dominant isoforms.** 34.8% of shared genes pick a
different dominant transcript in MCF7-SIX1 than in 231, so the 231 matrix would
hand many genes the wrong isoform's CDS. There is no salmon fallback for promoted
genes missing from the matrix, unlike `01_`: the MCF7-SIX1 quant files are in the
other repo and cover siCTRL *normoxia* only, so they cannot supply a hypoxia
dominant isoform. Those genes are reported and excluded rather than filled in
under the wrong condition.

**Set sizes are not comparable between cell lines.** 231 is spike-in normalized
and MCF7-SIX1 is not.

**Each reporter carries the tag in its own codons, so never transplant it between
constructs.** PEST GC3 is 52.5% in NanoLuc, 75.0% in Firefly, 82.5% in Renilla.
For NanoLuc the tag is substantially less GC3-biased than the ORF it is fused to
(51.2% vs 68.4%), which is why the full construct and the ORF alone are reported
separately.

## Related work elsewhere in the repo

- `code/predictive_modeling/62_cds_length_gc_hypoxia_promotes.Rmd` — the source of
  this analysis; carries Renilla and Firefly reference lines and the negative
  control comparator
- `code/predictive_modeling/64_nanoluc_codon_adaptation.Rmd` — bg-RSCU codon
  adaptation and the adapted / anti-adapted NanoLuc variant designs
- `code/predictive_modeling/68_cds_length_gc_red_beetle_luciferase.Rmd` — the same
  density framework for the red beetle luciferase
- `code/translational_control_overview/33a_parse_plasmid_sequences.R` — parses the
  lab plasmid constructs and cross-checks their CDS against notebook 64's output

## The codon basis every metric here uses

**Composition is measured on internal codons, matching the modelling.**
`trim_internal()` in `reporter_helpers.R` is ported verbatim from
`code/translational_control_overview/01_build_feature_matrix.Rmd`, which defines
it for `gc3_internal`, so the reporter design and the modelling describe one
quantity. Both trims are CONDITIONAL, and both conditions earn their place:

- the initiator is dropped **only when it really is ATG**. 25 of the 10,054
  expressed transcripts start on a near-cognate codon instead - CTG 14, GTG 6,
  ATT 2, ACG 1, ATA 1, TTG 1 - which is real initiation biology, not annotation
  noise, and they keep their first codon;
- the terminal stop is dropped **only when one is present**.

The principle is *trim only what you can positively identify*. Blindly removing
position 1 would silently delete a real codon from anything non-canonical.

**Transcripts with no terminal stop are dropped.** These are GENCODE
`cds_end_NF`: the CDS starts at ATG but the annotation runs out. 6,745 of
214,541 in-frame transcripts genome-wide, **53** inside the expressed
background, which is exactly the count CLAUDE.md documents as "in-frame and
passing every existing filter". They are compositionally skewed (GC3 median 48.9
against 57.2), so leaving them in drags the target down.

**This is a deliberate divergence from the modelling, which masks rather than
drops.** `01_build_feature_matrix.Rmd` sets the codon-derived features to NA for
these and keeps the rows, so they still contribute to length. Here the legend
prints one n per group across three panels, so masking would leave that n right
for length and silently overstated for GC and GC3. Dropping 53 of 10,107 costs
nothing measurable - the length median is unchanged at 1.28 kb.

**GC3 and total GC sit on different bases, and that is the modelling's pairing,
not an oversight.** GC3 is internal-codon (`gc3_internal`); total GC is the full
stop-stripped CDS, because the modelling has **no** internal variant of
`cds_gc` - so this folder's GC column stays comparable with every `cds_gc`
elsewhere in the project. CDS length also keeps the stop-stripped basis. The
helpers are named `orf_gc_pct` / `orf_gc3_pct` for this reason; an earlier
`gc_pct_internal` claimed a trim it no longer performed.

**Part-composition tables keep the plain basis** - they describe fragments as
ordered, not ORFs being compared against the gene distributions. `trim_internal`
is conditional, so on a part like P2A or the Waldo linker it is the identity.

**The initiator exclusion is ours, not the collaborator's.** The reconciliation
inventory in `00_feature_reconciliation.Rmd` records the gc3 difference as
"FIX: exclude the stop codon; keep percent units" - the initiator is not
mentioned, and CAI/FOP are logged `ours_only`, "keep, whole-CDS only", yet
`cai_internal`/`fop_internal` were built removing the initiator AUG anyway. The
one place "internal" is the collaborator's own word is `tai_genomic_internal`,
which was ADOPTed. Their source is not mirrored in this repo - only their data
tables are - so **whether they use the same conditional-on-ATG guard is unknown
and should not be asserted.**

**The literature convention is stricter than either.** Classical codon-bias work
drops every codon with no synonymous third-position choice, i.e. **all Met and
Trp**, not just the initiator - that is the GC3 vs **GC3s** distinction (CodonW;
Wright's Nc, *Gene* 1990), and the same logic behind CAI (Sharp & Li, *NAR*
1987) and RSCU (Sharp, Tuohy & Mosurski, *NAR* 1986) excluding single-codon
families. eGFP carries 6 Met and 1 Trp per copy, all forced GC-ending and all
still in the denominator here. Not adopted: consistency with the modelling was
the stated goal. Recorded because the rationale for going internal does not
actually stop at the initiator.

## Recoding the ladder onto the 231 composition (`08_`)

**The design is a set of counts: how many codons end in each of A, C, G and T.**
Synonymous choice moves only the third base of a codon, plus the first base for
Leu, Arg and Ser. So third-base composition is the only thing that can be set,
and total GC is not a second target - it is a consequence.

**GC3 alone is not enough, which is the non-obvious part.** A construct can hold
57.3% GC3 while being all C and no G. Whether a GC-ending residue offers G, C or
both is fixed by the amino acid: Glu, Lys, Gln, Met and Trp offer only G; Cys,
Asp, Phe, His, Ile, Asn and Tyr offer only C; only the 4- and 6-fold families
offer both. eGFP is Asp/Asn/Tyr/Phe/Ile-rich, so its C-only residues outnumber
its G-only ones **77 to 51**. A GC3-only recoder inherits that: an earlier
version landed at **G:C 0.65** against the transcriptome's 1.01, with C3 at the
70th percentile and G3 at the 25th. So GC3 sets the G+C total and the
transcriptome's own G:(G+C) and A:(A+T) ratios split each class.

The recoder is one pass with no iteration:

```
1. How many codons should end in each of A, C, G, T?
2. Give every position a third base its amino acid allows, in two stages:
   first decide which positions end G/C at all (this fixes GC3),
   then split those between G and C, and the rest between A and T.
3. At each position, pick a real codon with that third base,
   weighted by how often the 231 transcriptome uses it.
```

Delivered, against a target of GC3 57.3% (background n = 10,054):

| rung | nt | GC3 before | GC3 after | GC3 pctile | GC after | GC pctile | length pctile |
|---|---|---|---|---|---|---|---|
| SmBiT-P2A-FLAG | 117 | 71.0 | 57.9 | 51.0 | 50.4 | 46.5 | 0.3 |
| + link-1xeGFP | 876 | 91.1 | 57.4 | 50.1 | 48.7 | 40.8 | 31.2 |
| + link-2xeGFP | 1,635 | 92.5 | 57.4 | 50.0 | 48.6 | 40.4 | 63.7 |
| + link-3xeGFP | 2,394 | 93.0 | 57.3 | 50.0 | 48.7 | 40.6 | 81.2 |

**All four rungs sit at the GC3 median while length spans the 0.3rd to 81st
percentile.** That separation is the whole point: the ladder varies length
against a fixed, ordinary composition. All four third bases sit between the 20th
and 80th percentile on every eGFP rung.

**The bare 117 nt rung is effectively matched too, which the full-CDS basis hid.**
It reads 57.9% against the 57.3% target - 0.6 points, the 51st percentile. On
the old basis the same construct read 59.0% against 57.5%, and the gap was
attributed to SmBiT and P2A being 31 of its 39 codons. Most of it was really the
initiator ATG counted in a 39-codon denominator, where one codon is 2.6%. The
rung is still the one with least freedom, but it is not off target.

**Total GC lands at 48.7%, about 2.8 points below the 51.5% median**, and no
recoding can close it: eGFP's amino acids do not carry enough G/C in codon
positions 1 and 2. Against a GC IQR of ~14 points that is ~0.2 IQR. Reported,
not optimized.

**Stage the assignment, or GC3 quietly loses to the base balance.** An earlier
version assigned all four bases at once by swapping to reduce the *total*
deviation across A/C/G/T, which treats "3 codons too many on C" as
interchangeable with "3 codons short on A". On the bare rung it settled at GC3
66.7% when 59.0% was reachable at the same total deviation, because FLAG's Asp
and Tyr can only end C or T while its Lys can only end A or G. Deciding the
GC/AT split first and the within-class split second makes the priority
structural.

**Each eGFP copy is drawn independently**, which removes `07_`'s 723 nt and
1,446 nt exact tandem repeats; the copies share at most ~30 nt. The ladder stays
**nested**: rung *k* is built with rung *k-1* held fixed, so each rung is still
its predecessor plus one linker and one copy, and adding a copy never forces a
re-synthesis.

**Motif constraints are handled by redrawing, not repairing.** The Waldo
linker's `GAATTC` must stay the only EcoRI site and homopolymers are capped at
6. A fresh draw is cheap and cannot distort the composition the way an in-place
repair can.

**The skew check carries a discriminative companion.** Asserting all four bases
sit between the 20th and 80th percentile is worth something only if a skewed
sequence fails it, so the *original* eGFP construct is scored by the same
measure and asserted to fail: A at the 1st percentile, C at the 100th, T at the
2nd.

**`as.character(GENETIC_CODE)` drops the codon names.** Every amino-acid lookup
then returns `NA` silently rather than erroring. Restored explicitly and
asserted. `purrr::reduce` likewise needs its namespace - IRanges exports an S4
`reduce` that masks it and fails on a plain list.

## Which average is "the transcriptome"? (`09_`)

**"Match the 231 transcriptome" is not one number.** GC3 spans ~3.5 points
across five defensible ways of averaging, on siCTRL **normoxia** 1+4hr:

| weighting | GC3 | G:(G+C) |
|---|---|---|
| per transcript, median *(the target `08_` uses)* | 57.2 | 0.508 |
| per transcript, mean | 57.8 | 0.505 |
| pooled over codons | 55.9 | 0.503 |
| mRNA weighted (RNA-seq) | 59.4 | 0.494 |
| ribosome weighted (riboseq) | 57.3 | 0.514 |

**The ribosome-weighted answer lands on the unweighted one, not the mRNA one.**
The per-transcript target is therefore also, within a rounding error, the
composition the translation machinery actually encounters, and **mRNA weighting
is the outlier** rather than the more realistic summary. That defuses the
tRNA-competition argument for reweighting, so `08_` is unchanged. The
within-class split is indifferent throughout: G:(G+C) moves only 0.494 to 0.514,
about 1.5 codons in a 241-codon eGFP.

**Transcript-level riboseq quantification is unusable here; any future notebook
weighting by riboseq must sum to gene level first.** At transcript level **48 of
86 ribosomal protein transcripts carry exactly zero** footprint TPM; summed to
gene level the same genes sit 92-fold above the median and the top translated
genes are RPL35, TMSB10, PTMA, RPS3, TXN. Footprints cover the CDS only, so
isoforms differing in their UTRs are not separable and salmon's EM scatters
their reads. **RNA-seq passes the same check at transcript level**, and that
asymmetry is what pins the cause on footprints rather than on the isoform picks.
Section 3b keeps the transcript-level version as a worked counter-example and
asserts that it fails.

**A median per base is not a composition.** `median(G%) + median(C%)` is not
`median(G% + C%)`, so the median row's four bases sum to 98.65 rather than 100.
GC3 is normalised by each row's own total, or that row reads ~0.8 points low
against the others as a pure artefact of how it was summarised.
