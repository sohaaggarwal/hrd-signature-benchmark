PY
"""
Signature scoring functions.
 
All three signatures below are sourced from primary literature -- see
docs/provenance_table.md for full citations and known deviations/caveats.
"""
 
import pandas as pd
import numpy as np
 
 
# ---------------------------------------------------------------------------
# PARPi-7 (Peng et al., Breast Cancer Res Treat 2012, PMID 22875744)
# Direction of association with olaparib SENSITIVITY confirmed.
# NOTE: exact numeric weights from the paper's weighted-voting formula are
# NOT yet confirmed -- default below uses equal weighting as a placeholder.
# ---------------------------------------------------------------------------
 
PARPI7_GENES = {
    "resistance_associated": ["BRCA1", "MRE11", "NBN", "TDG", "XPA"],  # MRE11A -> MRE11 in current nomenclature
    "sensitivity_associated": ["CHEK2", "MAPKAPK2"],
}
 
 
def score_parpi7(expr_df: pd.DataFrame, weights: dict = None) -> pd.Series:
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
# 60-gene diagonal linear discriminant classifier, weights from Appendix Table A1.
# NOTE: gene #43 "P11" is an ambiguous/uncertain symbol from the source table.
# ---------------------------------------------------------------------------
 
BRCANESS_WEIGHTS = {
    "DAD1": 0.0997, "RAD21": 0.1743, "LDHA": 0.0165, "SPARC": 0.1571,
    "SKP1": 0.137, "PPP1CC": 0.0249, "RAN": 0.1043, "USP9X": -0.0066,
    "ECHS1": -0.0495, "GNAI3": 0.1404, "HDAC1": -0.0272, "LGMN": -0.0683,
    "CYR61": -0.0454, "SH3BGRL": 0.0659, "MGST3": -0.0162, "RCN2": 0.2087,
    "SMC1A": 0.0313, "BST2": 0.2118, "ENG": 0.058, "SRGN": -0.1577,
    "GBP1": 0.0778, "UNC119B": 0.139, "RGS1": -0.2573, "CDC2": 0.1101,
    "SAT1": -0.0386, "GGH": -0.0169, "GUCY1B3": 0.0248, "WFDC2": 0.0719,
    "NMI": 0.0873, "PRAME": 0.0267, "CCL4": -0.1545, "MGST2": -0.0045,
    "MMP7": 0, "SNCA": -0.1605, "VTN": -0.0399, "ALPP": -0.0274,
    "MTAP": -0.0511, "IDUA": 0.097, "SERPINF2": -0.0053, "WAS": -0.0191,
    "CD1D": -0.1431, "GFI1": -0.0318, "P11": 0.0091,  # ambiguous symbol, verify
    "SEMA3F": 0.0805, "TNF": -0.024, "ROS1": 0.1667, "MADCAM1": -0.0322,
    "PDIA4": 0.1564, "HMGN2": -0.0247, "HLA-B": 0.0956, "TM9SF1": 0.1125,
    "CCDC93": -0.0865, "APEX1": 0.0212, "VEGFA": 0.0612, "POSTN": 0.2142,
    "PSTPIP1": -0.041, "PMS1": -0.0579, "HLA-A": 0.0202, "PCTP": -0.0536,
    "SEH1L": 0.1127,
}
 
 
def score_brcaness(expr_df: pd.DataFrame) -> pd.Series:
    genes = list(BRCANESS_WEIGHTS.keys())
    present = [g for g in genes if g in expr_df.columns]
    missing = [g for g in genes if g not in expr_df.columns]
    if missing:
        print(f"Warning: {len(missing)} BRCAness genes missing from expression matrix: {missing}")
 
    z = (expr_df[present] - expr_df[present].mean()) / expr_df[present].std()
    score = sum(z[g] * BRCANESS_WEIGHTS[g] for g in present)
    return score
 
 
# ---------------------------------------------------------------------------
# CIN70 (Carter et al., Nat Genet 2006) -- negative control, NOT an HRD signature
# 70-rank / 72-gene list from Supplementary Table 1.
# NOTE: several symbols are outdated 2006-era gene names; some may not map
# to current DepMap gene symbols. Log missing genes rather than substituting.
# ---------------------------------------------------------------------------
 
CIN70_GENES = [
    "TPX2", "PRC1", "FOXM1", "CDC2", "C20orf24", "TGIF2", "MCM2", "H2AFZ",
    "TOP2A", "PCNA", "UBE2C", "MELK", "TRIP13", "CNAP1", "MCM7", "RNASEH2A",
    "RAD51AP1", "KIF20A", "CDC45L", "MAD2L1", "ESPL1", "CCNB2", "FEN1", "TTK",
    "CCT5", "RFC4", "ATAD2", "CKAP5", "NUP205", "CDC20", "CKS2", "RRM2",
    "ELAVL1", "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1",
    "OIP5", "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "NCAPH", "MCM10",
    "LSM4", "NCAPH2", "ASF1B", "ZWINT", "PBK", "FLJ10036", "CDCA3", "ECT2",
    "CDC6", "UNG", "MTCH2", "RAD21", "ACTL6A", "GPI", "MGC13096", "SFRS2",
    "HDGF", "NXT1", "NEK2", "DHCR7", "AURKA", "NDUFAB1", "KIAA0286", "KIF4A",
]
# Notes on renamed genes used above: ch-TOG->CKAP5, TOPK->PBK, STK6->AURKA,
# BRRN1->NCAPH, MTB->NCAPH2 (best-guess mappings, not independently verified).
 
def score_cin70(expr_df: pd.DataFrame) -> pd.Series:
    present = [g for g in CIN70_GENES if g in expr_df.columns]
    missing = [g for g in CIN70_GENES if g not in expr_df.columns]
    if missing:
        print(f"Warning: {len(missing)} CIN70 genes missing from expression matrix: {missing}")
 
    z = (expr_df[present] - expr_df[present].mean()) / expr_df[present].std()
    score = z.sum(axis=1)  # unweighted sum, per original paper's simple classifier approach
    return score
