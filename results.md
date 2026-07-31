# Results: Benchmarking Published HRD Gene Expression Signatures

## Question

Do published RNA-based homologous recombination deficiency (HRD) signatures
generalize to independent data? Specifically: do high signature scores predict
PARP inhibitor (olaparib) sensitivity in DepMap cancer cell lines, once cancer
lineage is controlled for?

## Data

- **Expression, cell line metadata, mutations:** DepMap Public 26Q1
- **Drug response (olaparib, IC50/AUC):** PRISM Secondary Repurposing 20Q2
- **Final analysis set:** 180 cell lines with both expression and olaparib
  response data, across 3 lineages: Ovary/Fallopian Tube (n=59), Pancreas (n=67),
  Breast (n=54)

## Signatures tested

| Signature | Source | Genes | Type |
|---|---|---|---|
| PARPi-7 | Peng et al., *Breast Cancer Res Treat* (2012) | 7 | HRD-specific, hypothesis-driven panel |
| BRCAness | Konstantinopoulos et al., *JCO* (2010) | 60 | HRD-specific, discriminant classifier |
| CIN70 | Carter et al., *Nat Genet* (2006) | 72 | Negative control (general proliferation/chromosomal instability, NOT HRD-specific) |

Full provenance, exact gene lists, and known deviations documented in
`docs/provenance_table.md`.

## Methods

1. Computed each signature's score per cell line from log2(TPM+1) expression data.
2. Tested Spearman correlation between signature score and olaparib AUC,
   overall and within each lineage.
3. Fit a lineage-adjusted regression (`auc ~ signature_score + C(lineage)`)
   for each signature.
4. Ran a permutation test: compared each real signature's correlation against
   500 random gene sets of matched size, computing an empirical p-value.
5. Tested inter-signature correlation, to assess whether the three signatures
   agree with each other.
6. Tested BRCA1/2 mutation status (likely-LoF variants only) as a baseline
   predictor, using the same lineage-adjusted regression structure.
7. Repeated the core regression using IC50 as an alternate drug-response
   metric, as a robustness check.

## Results

### 1. Basic correlation (no lineage adjustment)

No signature showed a significant correlation with olaparib AUC, overall
(n=180) or within any single lineage.

| Signature | Overall rho | Overall p |
|---|---|---|
| PARPi-7 | -0.019 | 0.81 |
| BRCAness | 0.003 | 0.96 |
| CIN70 | -0.068 | 0.37 |

### 2. Lineage-adjusted regression

Adjusting for lineage did not change the conclusion. All p-values remained
well above 0.05; full-model R² was low (0.04-0.05) for all three signatures.

| Signature | Coefficient | p-value | R² (full model) |
|---|---|---|---|
| PARPi-7 | 0.0025 | 0.546 | 0.046 |
| BRCAness | 0.0211 | 0.181 | 0.053 |
| CIN70 | -0.0004 | 0.274 | 0.050 |

### 3. Permutation testing

This was the most decisive check. For each signature, we compared its real
correlation against 500 random gene sets of the same size, drawn from the
same expression data.

| Signature | Real rho | Empirical p-value | Interpretation |
|---|---|---|---|
| PARPi-7 | -0.019 | 0.846 | 85% of random gene sets performed as well or better |
| BRCAness | 0.003 | 0.958 | 96% of random gene sets performed as well or better |
| CIN70 | -0.068 | 0.366 | 37% of random gene sets performed as well or better |

None of the three signatures outperformed random gene sets of matched size.

### 4. Do the signatures agree with each other?

| Pair | rho | p-value |
|---|---|---|
| PARPi-7 vs BRCAness | -0.208 | 0.005 |
| PARPi-7 vs CIN70 | -0.282 | 0.0001 |
| BRCAness vs CIN70 | **0.419** | **<0.0001** |

BRCAness correlates strongly with CIN70 -- a signature that is explicitly
*not* designed to detect HRD, only general proliferation/chromosomal
instability. This is consistent with BRCAness substantially capturing
non-HRD-specific genomic instability rather than the specific DNA-repair
defect it claims to measure. (PARPi-7's negative correlations with the other
two are likely a sign-convention artifact -- PARPi-7 is oriented toward
"sensitivity," while BRCAness/CIN70 trend toward "instability" -- rather than
biological disagreement, and should not be over-interpreted.)

### 5. BRCA1/2 mutation status baseline

Of 180 cell lines, 24 carried a likely-damaging BRCA1/2 mutation. This
gold-standard genetic marker also failed to predict olaparib response in
this dataset:

| Test | Coefficient | p-value |
|---|---|---|
| Simple correlation | rho = -0.067 | 0.369 |
| Lineage-adjusted regression | -0.008 | 0.833 |

### 6. Robustness check: IC50 vs AUC

Repeating the lineage-adjusted regression using IC50 (n=70, smaller subset
with IC50 data available) instead of AUC:

| Predictor | Coefficient | p-value |
|---|---|---|
| PARPi-7 | -0.137 | 0.189 |
| BRCAness | +0.788 | **0.048** |
| CIN70 | +0.009 | 0.428 |
| BRCA-mutant | -0.477 | 0.613 |

BRCAness crosses the conventional p<0.05 threshold with IC50, but this
should not be read as a positive result: (a) the coefficient sign is
*reversed* from what the signature predicts (higher BRCAness -> higher IC50
-> *less* sensitive, not more), (b) the sample size is much smaller (n=70
vs 180), and (c) 8 total tests were run across both metrics, so one
borderline result is consistent with chance (~1 in 20 expected) rather than
signal.

## Conclusion

Across three independent statistical approaches (basic correlation,
lineage-adjusted regression, and permutation testing against random gene
sets), **none of the three published HRD expression signatures showed a
robust, direction-consistent relationship with PARP inhibitor sensitivity**
in this independent DepMap/PRISM dataset.

Critically, the gold-standard genetic marker (BRCA1/2 mutation status) also
failed to predict response in this same dataset. This is an important
caveat: it means we cannot fully separate "these RNA signatures do not
generalize" from "this specific drug-response assay may lack sufficient
sensitivity to detect even a known, validated effect." A dataset in which a
true positive control also returns null limits how strongly we can
conclude the signatures themselves are invalid, versus the assay being
underpowered or noisy for this purpose.

The strongest and most interpretable finding is the inter-signature
correlation result: BRCAness's substantial overlap with a non-HRD
proliferation signature (CIN70) offers a plausible mechanistic explanation
for its lack of specific predictive power, independent of any assay
limitations.

## Limitations

- Cell-line drug response (in vitro) does not directly equate to patient
  clinical benefit.
- PRISM's pooled screening methodology differs from single-agent dose-response
  assays used in the original signature-development papers.
- Per-lineage sample sizes (54-67 cell lines) are modest; some effects may be
  underpowered.
- The BRCA-mutant baseline null result limits confidence that this dataset
  can detect known HRD-related drug sensitivity effects at all.
- CIN70's gene list required best-guess mapping of several outdated (2006-era)
  gene symbols to current nomenclature; a handful of genes could not be
  mapped and were excluded (see `docs/provenance_table.md`).
- PARPi-7 was scored using equal weighting as a placeholder; the original
  paper's exact numeric weights were not independently confirmed.

## Reproducing this analysis

See `README.md` for the full pipeline. Code: `merge_data.py` ->
`scores_signatures.py` -> `correlation_analysis.py` -> `lineage_regression.py`
-> `permutation_test.py` -> `signature_agreement.py` -> `brca_baseline.py` ->
`robustness_check.py`.