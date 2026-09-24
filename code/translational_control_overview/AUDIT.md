# Feature-matrix audit: findings, fixes, and what is now guaranteed

Scope: `output/translational_control_overview/feature_matrix_tco.rds` (10,172 transcripts x
173 columns), built by `01_build_feature_matrix.Rmd`. Nothing in
`code/predictive_modeling/` was changed; that pipeline remains frozen.

## 1. How the audit was done

An independent agent recomputed features from primary sources — GENCODE v49 TxDb,
`BSgenome.Hsapiens.UCSC.hg38`, and the published reference tables — with its own
implementation, and re-joined every reference table by hand. It did not re-run our notebooks
and did not rely on our assertions, because a check written by the author of a feature cannot
independently confirm that feature.

That distinction turned out to be the whole point. **Every defect below passed every check we
had at the time.** None was found by an assertion.

## 2. Findings

All seven are fixed. Severity as assessed at the time of discovery.

| # | Defect | Evidence | Fix |
|---|---|---|---|
| **F1** | 1,689 transcripts had every `tco_struct_*` column **zero-filled instead of NA** (`replace_na(., 0)` in `01d_`/`01k_`). Zero reads as "maximally structured"; because the NA was gone, the sparsity filter and median imputation both missed them | the zero set is *identical* to the set where G4mer is NA — same missing transcripts, handled correctly there | restored to NA on an unambiguous sentinel: region mean 0 **and** sd 0 |
| **F2** | `tco_struct_num_structured_regions*` is region length in disguise | Spearman **0.99998** with `cds_length`; the "count" *exceeds* the length by 1 for **7,653** transcripts, impossible for a count of positions | deprecated |
| **F3** | `hia2026_dhx29_occupancy` joined at **gene** level from a transcript-keyed source, so ~1,013 transcripts carry a different isoform's occupancy | source has 9,740 transcript-keyed rows; transcript join gives 8,529; stored column had 9,542 non-NA | `hia2026_dhx29_occupancy_tx` added as the corrected variant |
| **F4** | `clip_lee_complete_media_bound` is **bit-identical** to `clip_lee_glu_dep_bound` (already identical upstream in `score_card.csv`), so `clip_cap_binding_count` counted one dataset twice | all 356 transcripts with count == 2 are that pair | duplicate deprecated; `clip_cap_binding_count_v2` counts the 2 distinct conditions |
| **F5** | `tco_struct_accessibility_*_min` identically 0 in all 10,172 rows — `min()` taken after NA became 0 | sd exactly 0 | deprecated; `02_` now drops zero-variance columns generally |
| **F6** | six `tco_struct_*` CDS/3'UTR columns filed under `structure_utr5`: `assign_family()` anchored the region tag with `$`, but those names end in `_mean`/`_min`/`_sd` | 12 columns in `structure_utr5`, only 6 of them 5'UTR | regex matches the tag anywhere; asserted 6/6/6 |
| **F7** | `te_lfc` / `te_padj` / `te_lfc_bin` are the **siCTRL hypoxia-vs-normoxia** contrast, not a knockdown contrast | bit-for-bit `translation_categories_1hr.csv` (r = 1.0000, max dev 0); EIF3D's own `te_lfc` is **negative**; only 11 of 1,446 `sig_hypoxia_3d_promotes` transcripts have `te_lfc > 0.5` | no model impact (`block = label`, asserted out of features), but the identity is now pinned by assertion and stated in the definition |

### Effect on results already reported

- **F6 corrupted the SHAP family attribution.** "`structure_utr5` = 9.97%" was a mixture of
  all three regions. Corrected: **CDS 5.78%, 5'UTR 2.62%, 3'UTR 1.49%** — CDS structure is the
  dominant structural contribution, not 5'UTR.
- The model is **66 features** after the deprecations (was 72).
- Single-run test AUC rose 0.8083 → 0.8276, but **one draw means nothing here**: the negative
  draw alone moves test AUC by ~0.07. This needs the paired sweep before it is quotable.

## 3. What now guarantees this does not recur

Each defect was a *class*, so each became a standing check rather than a one-off repair.

### Structural invariants (`01_`, all asserted, all passing)

| | catches | status |
|---|---|---|
| **A** no two modelled features bit-identical | F4 | 0 |
| **B** no modelled feature constant | F5 | 0 |
| **C** no near-perfect length proxy outside the `length` family | F2 | 0 |
| **D** no unexplained sentinel block (one exact value over >5% of rows), with a named allowlist | F1 | 0 unexpected |
| **E** corrected occupancy matches a fresh transcript-level join against its source | F3 | max dev 0 |

Invariant E is stated carefully, because the obvious version of it does not work: this matrix
holds **exactly one transcript per gene**, so there are no sibling isoforms to compare and a
gene-level join leaves no duplicate-value fingerprint. That is precisely why F3 was silent.
The only real check is against the source file.

### Mutation testing (`01c_mutation_test_checks.R`) — 12 of 12 caught

A check that has never been shown to **fail** is uninformative: the `01h_` Kozak PWM passed its
own tier-monotonicity check for years while being wrong at position −3. So each check is now
proven against a deliberate corruption of its input. A mutation run never writes the matrix.

`kozak_minus3_bug` · `noderer_window_shift` · `cridge_permute` · `term_stop_shuffle` ·
`peptide_frame_shift` · `g4_region_swap` · `struct_zero_fill` · `gene_level_join` ·
`duplicate_feature` · `constant_feature` · `length_proxy_feature` · `golden_drift`

The suite earned its place on its first run: `g4_region_swap` (exchange the 5'UTR and 3'UTR
G4mer summaries) originally passed **every** check, because the motif spot-check asserted only
an ordering *within* each region and G-richness is correlated across regions of the same
transcript. The CDS/5'UTR-vs-3'UTR result rested on a join nothing verified. It now fails
against a discriminative check.

### Golden transcripts (`golden_transcripts.csv`, 20 x 99)

A frozen snapshot, deliberately including awkward cases: the named spot-check genes, the F1
group with NA structure, truncated CDS, and non-AUG starts.

**A snapshot alone is a regression detector, not a correctness detector** — had it been made
before the audit it would have frozen the zero-filled values and asserted them forever. So each
feature carries how far it has actually been checked:

| status | features | meaning |
|---|---|---|
| `recomputed` | **24** | independently re-derived from GENCODE + hg38 + reference tables. **0 numeric disagreements** across 480 values |
| `traced` | 10 | matched by key to its source file; provenance verified, measurement not |
| `unverified` | 65 | neither possible from here |

The NA test is asymmetric because the directions mean opposite things. **Builder holds a value
the primary data cannot support: 0** (the dangerous direction). Builder masks where
recomputation succeeds: 36, every one explained by a truncated CDS or an absent 3'UTR.

### Reproduction gates

The repairs change values, so **bit-exact baseline reproduction is deliberately no longer the
default**. `01_` takes `apply_audit_repairs` (default TRUE); build with it **FALSE** to
reproduce the historical matrix and re-run `02_`'s gate. `02_` skips its zero-variance drop on
a gate run for the same reason. Both gates are green on that path.

## 4. Known limitations

- **A recomputation written by the same author shares that author's assumptions.** It proves
  internal consistency, not truth. Independent re-derivation by someone else is the only thing
  that has actually found defects here — two audits, seven confirmed, versus zero found by
  assertions.
- **65 of 99 modelled features remain unverified**, chiefly the inherited `clip_*` (23),
  `tco_struct_*`, and codon-optimality columns. Their checks live in their original
  `code/predictive_modeling/` notebooks; `FEATURES.md`'s `assigned in` column names where.
- **RNAplfold and G4mer values themselves were never verified** — both were computed off-repo
  and only summary files are present. Joins and derived columns are verified; the underlying
  numbers are not.
- **RESOLVED — the RNAplfold parameters differ by region, which neither document said.**
  Verified by re-running RNAplfold 2.7.2 on ENST00000000233 and reproducing the stored
  `_lunp` files bit-for-bit:

  | region | parameters |
  |---|---|
  | 5'UTR | `-W 80 -L 40 -u 30` |
  | CDS | `-W 150 -L 100 -u 30` |
  | 3'UTR | `-W 150 -L 100 -u 30` |

  The annotation claimed `-W 150 -L 100` for all three; `01d_`'s inline docs claimed
  `-W 80 -L 40 -u 30` for all three. Each was right about its own region. **A larger window
  sees more folding context, so accessibility is not comparable across regions** — the same
  caveat the `g4` family carries for its differing stride. Any statement ranking CDS
  structure against 5'UTR structure (e.g. the SHAP shares 5.78% vs 2.62%) is confounded by
  window size and must not be read as biology. The annotation now records the parameters and
  the non-comparability per column.

- **RESOLVED — the reported `_lunp` off-by-one is NOT confirmed.** A cut-off audit claimed
  one without evidence. The primary accessibility column is `df[[2]]`, the `l=1` column, and
  for a 1-mer "ending at position i" *is* position i, so that parse is correct. What is real
  is a different defect in the `l=30` column: it is NA for rows 1-29 (a 30-mer cannot end
  before position 30) and the parser zero-fills those, which reads as "maximally structured".

- **Root cause of F2 established.** The "structured region" count thresholds the 30-nt window
  at `< 0.2`, but the probability of a full 30-mer being unpaired is around **1e-9**, so
  **100% of positions pass**. The feature never counted structure; it counted positions, which
  is why it equalled region length. Already deprecated.
- **Not re-derived at all:** `kozak_score` / `kozak_optimal` tiers, `cnot3_weighted_codon_score`,
  positional CSC, and `clip_lee_glu_vs_cm_enrichment` — which is *not*
  `glu_dep_lfc − complete_media_lfc` despite the name (max dev 21.8) and whose definition is
  unknown.

## 5. Reproducing all of it

```
Rscript code/translational_control_overview/01b_feature_dictionary.R    # FEATURES.md
Rscript code/translational_control_overview/01c_mutation_test_checks.R  # 12 mutations, ~28 min
Rscript code/translational_control_overview/01d_golden_transcripts.R    # snapshot + recompute
```

`01_build_feature_matrix.Rmd` runs every invariant and the golden check on each knit.
