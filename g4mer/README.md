# G4mer rG4 scoring

RNA G-quadruplex prediction across 5'UTR, CDS, and 3'UTR using
[G4mer](https://huggingface.co/Biociphers/g4mer) ([Nat Commun 2025](https://www.nature.com/articles/s41467-025-65020-7)),
a 6-mer tokenized RNA language model fine-tuned from mRNAbert.

Follows the same pattern as `cluster_scripts/viennarna/`: an R script generates
FASTA, a non-R tool scores it, results return as flat files, and a notebook ingests
them. Scripts are tracked in git; score tables are not (`output/` is gitignored).

## Model constraints that shape the design

- **Maximum 70 nt per input.** Not 512. Longer sequences must be tiled.
- The authors recommend scanning with a **70 nt sliding window** and taking the
  **maximum** score across windows. That is the default aggregation here.
- Output is a 2-class softmax; `P(rG4)` is index 1.

Because the aggregation is a maximum over windows, `g4mer_max` rises with region
length by construction, and rG4 formation requires G-runs so it also tracks GC and
G content. **Any downstream association must control for both.** Notebook `01u_`
reports every result raw and after residualizing on log2 length and GC%.

## Transcript universe

All three regions subset to `output/predictive_modeling/precomputed_most_abundant_tx.rds`
(11,879 transcripts), the same set used for the ViennaRNA structure features. This
means G4mer scores join to the structure features on `transcript_id_clean` with no
reconciliation. Do not substitute a different transcript selection.

Region counts after the `>= 10 nt` filter:

| region | n | median nt | max nt | no annotation |
|---|---|---|---|---|
| 5'UTR | 9,545 | 129 | 3,347 | 2,213 (18.6%) |
| CDS | 10,025 | 1,227 | 20,724 | 1,854 (15.6%) |
| 3'UTR | 9,697 | 763 | 17,423 | 2,165 (18.2%) |

## One-time setup

The model is a **gated** HuggingFace repo. Accept the terms at
https://huggingface.co/Biociphers/g4mer, then authenticate from a normal terminal
(not through an agent session — the token should not enter a transcript):

```bash
conda create -y -n g4mer python=3.10
conda activate g4mer
pip install torch transformers
hf auth login          # 'huggingface-cli' is deprecated in huggingface_hub 1.x
```

## Running

Generate FASTA for all three regions (needs the v49 TxDb from
`code/predictive_modeling/00_build_txdb_v49.R`):

```bash
Rscript g4mer/01_generate_fasta.R
```

Score the 5'UTR locally — this is the primary analysis and small enough for a laptop:

```bash
conda run -n g4mer python g4mer/02_run_g4mer.py --region utr5 --selftest
```

`--selftest` scores the model card's known rG4 example and hard-fails below 0.5.
It is a wiring check for the model id, tokenization, and softmax index — not a
validation of the predictions. Use `--limit 50` for a fast pilot.

Score CDS and 3'UTR on the cluster (~31 Mb of sequence, ~3M windows at stride 10):

```bash
sbatch g4mer/03_run_g4mer_cluster.sbatch
```

Then copy `g4mer_summary_*.tsv` and `g4mer_windows_*.tsv.gz` back into `output/g4mer/`.

## Deviation from the reference implementation

The authors' own `predict.py` (Bitbucket `biociphers/g4mer`, `src/g4mer/predict.py`
and `tutorial/src/predict.py`) tiles with **stride 1** — every 1 nt offset — and
collapses with `groupby('ind')['g4mer'].idxmax()`. The 70 nt window and the max
aggregation here match that. **The stride does not:** this pipeline uses 10 (5'UTR)
and 20 (CDS, 3'UTR) to keep the compute tractable. `--stride` is not part of the
authors' design.

Measured on 400 randomly sampled 5'UTRs (`set.seed(9)`), scored at stride 1 and
compared against the project values:

| | median `g4mer_max` | Spearman vs stride 1 | Pearson | rG4 calls at 0.5 that flip |
|---|---|---|---|---|
| stride 1 (reference) | 0.290 | - | - | - |
| stride 10 | 0.211 | 0.9924 | 0.9860 | 2.2% |
| stride 20 | 0.187 | 0.9786 | 0.9652 | 4.0% |

**Implication:** rank-based analyses (Spearman, quartiles, Wilcoxon) are unaffected.
Absolute scores are systematically ~0.08-0.10 lower than the reference and should not
be compared against published values or the web tool output.

Cost of matching the reference: stride 1 took 11.2 min for 400 5'UTRs on CPU, so the
full 5'UTR set is roughly 4.5 h and CDS + 3'UTR roughly 30M windows. A GPU node makes
it feasible if exact reproduction is ever needed.

## Stride used per region

The 5'UTR was scored locally at **stride 10**; CDS and 3'UTR were scored on the
cluster at **stride 20**. Denser sampling raises the expected maximum, so
`g4mer_max` is not directly comparable across regions. Notebook `01u_` recomputes the
5'UTR maximum on the stride-20 subset of its saved per-window scores
(`g4mer_max_harmonized`) for any cross-region comparison. Within-region results are
unaffected: Spearman between the stride-10 and stride-20 5'UTR maxima is 0.985,
though the median shifts 0.189 to 0.152.

## Outputs

`output/g4mer/g4mer_summary_<region>.tsv` — one row per transcript:
`transcript_id_clean, region, region_length, n_windows, g4mer_max, g4mer_mean,
g4mer_frac_above, g4mer_argmax_start`

`output/g4mer/g4mer_windows_<region>.tsv.gz` — one row per window:
`transcript_id_clean, region, window_start, window_end, g4_prob`

Per-window scores are kept in full so alternative aggregations (mean, fraction above
threshold, position-resolved) can be explored without re-running inference.

## Notes

- No MPS on the current Mac: torch 2.11 reports `mps.is_available() == False` on
  macOS 13 (Darwin 22.6). The 5'UTR run is CPU-only, roughly 20 min at ~115 win/s.
  This is the main argument for putting CDS/3'UTR on a GPU node.
- `extractTranscriptSeqs()` in `01_generate_fasta.R` is strand-corrected and
  exon-spliced. Do not replace it with `bedtools getfasta` on genomic intervals,
  which would need manual strand handling and exon concatenation.
- Installed `transformers` is 5.x; the model card example targets 4.x. If model
  loading breaks after an upgrade, pin `pip install 'transformers<5'`.

## Downstream

`code/predictive_modeling/01u_feature_extraction_g4mer.Rmd` validates the scores,
writes `output/predictive_modeling/feature_matrix_g4mer_utr5.rds`, and tests the
association between 5'UTR rG4 propensity and eIF3d-dependent TE in hypoxia at 1hr.
