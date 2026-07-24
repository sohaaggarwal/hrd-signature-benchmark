# Signature Provenance Table

This table documents each candidate signature's origin. All three are now fully sourced from primary literature.

| Signature | Source publication | Training data / platform | Endpoint predicted | Gene list status | Notes |
|---|---|---|---|---|---|
| **PARPi-7** | Peng et al., *Breast Cancer Res Treat* (2012). PMID: 22875744 | 22 breast cancer cell lines, olaparib response, cross-platform mRNA | Olaparib sensitivity (continuous) | CONFIRMED. 7 genes: BRCA1, MRE11A, NBN, TDG, XPA (resistance-associated); CHEK2, MAPKAPK2 (sensitivity-associated) | Deliberately small, hypothesis-driven panel from 118 candidate DNA-repair genes. Exact numeric weights (vs. equal weighting) still need confirming from paper before final scoring. MRE11A appears as MRE11 in current DepMap release (updated nomenclature). |
| **BRCAness (Konstantinopoulos)** | Konstantinopoulos et al., *J Clin Oncol* 28:3555-3561 (2010). PMID: 20547991 | 61 EOC patients (34 BRCA1/2-mutant, 27 sporadic), microarray | BRCA-like classification; platinum/PARPi response | CONFIRMED. 60 genes + weights, full list in `src/signatures.py` (from Appendix Table A1) | Diagonal linear discriminant classifier, 94% cross-validation accuracy in original paper. One gene symbol (#43, "P11") is ambiguous/uncertain from the source table and needs manual double-check. |
| **CIN70** | Carter et al., *Nat Genet* 38:1043-1048 (2006) | Multiple cancer types, 3 datasets (breast, ovarian, small-cell lung) used for gene ranking | NOT an HRD signature — proliferation/chromosomal instability control | CONFIRMED. 72 gene symbols (70 ranks, 2 ranks each cover 2 genes), full list in `src/signatures.py` (from Supplementary Table 1) | Deliberately included as a negative control: if this "wrong" signature predicts PARPi response as well as the HRD-specific ones, that's evidence the real signatures are tracking general genomic instability, not HRD biology specifically. Several gene symbols are outdated 2006-era names (e.g. STK6→AURKA, ch-TOG→CKAP5, TOPK→PBK) and some (FLJ10036, MTB, KIAA0286) may not map cleanly to current DepMap gene symbols — expect some genes to come up missing when scoring; log any that do rather than guessing a substitute. |

## Rules for this table
- No signature is scored until its row is fully filled in with a verifiable source.
- Any deviation (missing gene, outdated nomenclature, platform mismatch, ambiguous weight) gets logged here, not silently patched.
