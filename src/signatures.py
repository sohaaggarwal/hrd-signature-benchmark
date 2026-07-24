"""
Signature scoring functions.

IMPORTANT: Do not fill in gene lists or weights here from memory or assumption.
Every value must trace back to the source cited in docs/provenance_table.md.
"""

import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# PARPi-7 (Peng et al., Breast Cancer Res Treat 2012, PMID 22875744)
# CONFIRMED gene list -- direction of association with olaparib SENSITIVITY:
#   resistance-associated (negative weight): BRCA1, MRE11A, NBN, TDG, XPA
#   sensitivity-associated (positive weight): CHEK2, MAPKAPK2
#
# NOTE: gene list and reported direction are confirmed, but the exact numeric
# weights from the original paper's weighted-voting formula have NOT yet been
# pulled. Before using this for real scoring, check the original paper
# (PMID 22875744) and confirm:
#   1. exact weight for each gene (or whether it's an unweighted vote)
#   2. exact z-scoring / normalization method used before combining
# Do not assume equal weights -- that is a simplification, not a documented fact.
# ---------------------------------------------------------------------------

PARPI7_GENES = {
    "resistance_associated": ["BRCA1", "MRE11A", "NBN", "TDG", "XPA"],
    "sensitivity_associated": ["CHEK2", "MAPKAPK2"],
}


def score_parpi7(expr_df: pd.DataFrame, weights: dict = None) -> pd.Series:
    """
    Compute PARPi-7 score per sample.

    Parameters
    ----------
    expr_df : DataFrame, samples x genes, values = log2(TPM+1) or similar
    weights : dict of gene -> weight. If None, uses PLACEHOLDER equal weighting
              (+1 for sensitivity-associated, -1 for resistance-associated) --
              THIS IS A SIMPLIFICATION until real weights are confirmed from
              the source paper. Do not treat this default as validated.

    Returns
    -------
    Series of scores, higher = predicted more sensitive to olaparib
    """
    genes = PARPI7_GENES["sensitivity_associated"] + PARPI7_GENES["resistance_associated"]
    missing = [g for g in genes if g not in expr_df.columns]
    if missing:
        raise ValueError(f"Missing PARPi-7 genes in expression matrix: {missing}")

    if weights is None:
        weights = {g: 1 for g in PARPI7_GENES["sensitivity_associated"]}
        weights.update({g: -1 for g in PARPI7_GENES["resistance_associated"]})

    z = (expr_df[genes] - expr_df[genes].mean()) / expr_df[genes].std()
    score = sum(z[g] * weights[g] for g in genes)
    return score


# ---------------------------------------------------------------------------
# BRCAness (Konstantinopoulos et al., JCO 2010, PMID 20547991)
# STATUS: NOT YET IMPLEMENTABLE.
# Need: the 60 gene symbols + diagonal linear discriminant weights from
# Appendix Table A1 of the original paper. Fill in below before use.
# ---------------------------------------------------------------------------

BRCANESS_GENES = None  # TODO: fill in from Appendix Table A1

def score_brcaness(expr_df: pd.DataFrame) -> pd.Series:
    raise NotImplementedError(
        "BRCAness gene list not yet sourced from Konstantinopoulos et al. 2010 "
        "Appendix Table A1. See docs/provenance_table.md."
    )


# ---------------------------------------------------------------------------
# CIN70 (Carter et al., Nat Genet 2006) -- negative control, NOT an HRD signature
# STATUS: NOT YET IMPLEMENTABLE. Need gene list from original paper.
# ---------------------------------------------------------------------------

CIN70_GENES = None  # TODO: fill in from Carter et al. 2006

def score_cin70(expr_df: pd.DataFrame) -> pd.Series:
    raise NotImplementedError(
        "CIN70 gene list not yet sourced from Carter et al. 2006. "
        "See docs/provenance_table.md."
    )
