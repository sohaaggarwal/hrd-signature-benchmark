# Benchmarking Published HRD Gene Expression Signatures Against PARP Inhibitor Sensitivity

**Question:** Do published RNA-based homologous recombination deficiency (HRD) signatures
generalize to independent data — i.e., do high signature scores actually predict PARP
inhibitor sensitivity in DepMap cancer cell lines, once cancer lineage is controlled for?

This is a reproducibility / benchmarking study, not a new-signature discovery project.

---

## Exact step-by-step plan

### Phase 0 — Setup (do this first, before touching data)
1. [ ] Fill in docs/provenance_table.md completely for all 3 signatures — get the
   BRCAness (Konstantinopoulos 2010, Appendix Table A1) and CIN70 (Carter 2006) gene
   lists from the primary sources. Do not proceed to Phase 1 until this is done.
2. [ ] Set up Python environment (see requirements.txt — pandas, numpy,
   scipy, statsmodels, scikit-learn).
3. [ ] Create a GitHub repo, push this skeleton, add a LICENSE and this README.

### Phase 1 — Data acquisition
4. [ ] Download DepMap data (current release) from depmap.org/portal/download:
   - Gene-level RNA-seq expression file, log2(TPM+1)
   - Cell line metadata file (lineage, subtype, BRCA1/2 mutation status)
   - Somatic mutations file, for BRCA1/2 mutation status if not in metadata
   - Drug sensitivity data: PRISM or GDSC screen results for olaparib (primary) and
     one secondary PARPi (rucaparib or niraparib — check current coverage first)
5. [ ] Verify exact current DepMap release version and file names before downloading —
   these change between releases.
6. [ ] Document DepMap release version, download date, and file checksums in
   docs/data_provenance.md — standard reproducibility practice.

### Phase 2 — Signature reconstruction
7. [ ] Implement each signature's exact scoring formula in src/signatures.py.
8. [ ] For each signature, check what fraction of its genes are present in the DepMap
   expression matrix. Log any missing genes in the provenance table.
9. [ ] Compute a score per cell line per signature.

### Phase 3 — Primary association analysis
10. [ ] Spearman correlation: signature score vs. drug response, across all lineages.
11. [ ] Repeat within each of your 3 chosen lineages separately (ovarian, breast, +1).

### Phase 4 — Lineage-adjusted model (the methodologically honest core result)
12. [ ] Fit: drug_response ~ signature_score + lineage. Compare to the unadjusted
    correlation from Phase 3.

### Phase 5 — Permutation null model
13. [ ] Generate ~1,000 random gene sets matched on mean expression and variance.
14. [ ] Score cell lines with each random set the same way as the real signature.
15. [ ] Compare the real signature's association strength to this null distribution.

### Phase 6 — BRCA-mutant baseline
16. [ ] Test BRCA1/2 mutation status alone as a predictor, using the same
    lineage-adjusted model structure.
17. [ ] Compare: does any expression signature outperform this simple genetic marker?

### Phase 7 — CIN70 specificity check
18. [ ] Run CIN70 through the identical pipeline (Phases 2-6).

### Phase 8 — Robustness + writeup
19. [ ] Sensitivity check: does the result hold with AUC vs IC50? With/without outliers?
20. [ ] Write up methods, results, and limitations.
21. [ ] Push final analysis to GitHub. Consider a bioRxiv preprint once complete.

---

## Repo structure
