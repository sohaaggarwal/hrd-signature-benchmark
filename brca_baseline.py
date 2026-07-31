import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

print("Loading mutations file (this is the big 580MB one, may take a bit)...")
mut = pd.read_csv("data/OmicsSomaticMutations.csv",
                   usecols=["ModelID", "HugoSymbol", "LikelyLoF"])
print(f"Loaded {mut.shape[0]} mutation records\n")

# Filter to likely damaging BRCA1/BRCA2 mutations
brca_mut = mut[
    (mut["HugoSymbol"].isin(["BRCA1", "BRCA2"])) &
    (mut["LikelyLoF"] == True)
]
brca_mutant_ids = set(brca_mut["ModelID"].unique())
print(f"Found {len(brca_mutant_ids)} cell lines with a likely damaging BRCA1/2 mutation\n")

print("Loading scored signature data...")
df = pd.read_csv("data/scored_signatures.csv")
df["brca_mutant"] = df["ModelID"].isin(brca_mutant_ids).astype(int)
print(f"Of {df.shape[0]} rows, {df['brca_mutant'].sum()} are BRCA-mutant\n")

print("=" * 60)
print("PHASE 6: BRCA-mutant status as a baseline predictor")
print("=" * 60)

# Simple comparison: does BRCA-mutant status alone correlate with AUC?
valid = df[["brca_mutant", "auc"]].dropna()
rho, pval = stats.spearmanr(valid["brca_mutant"], valid["auc"])
print(f"\nSimple correlation (BRCA-mutant vs AUC):")
print(f"  rho = {rho:.4f}, p = {pval:.4g}, n = {len(valid)}")

# Lineage-adjusted regression, same structure as Phase 4
df["lineage_clean"] = df["OncotreeLineage"].str.replace("[^A-Za-z0-9]", "_", regex=True)
valid = df[["brca_mutant", "auc", "lineage_clean"]].dropna()
model = smf.ols("auc ~ brca_mutant + C(lineage_clean)", data=valid).fit()

print(f"\nLineage-adjusted regression (BRCA-mutant vs AUC):")
print(f"  Coefficient on brca_mutant: {model.params['brca_mutant']:.5f}")
print(f"  P-value on brca_mutant:     {model.pvalues['brca_mutant']:.4g}")
print(f"  R-squared (full model):     {model.rsquared:.4f}")
print(f"  N = {len(valid)}")

print("\n" + "=" * 60)
print("Comparison to signature scores (from Phase 4):")
print("  parpi7_score p-value:   0.546")
print("  brcaness_score p-value: 0.181")
print("  cin70_score p-value:    0.2741")
print("=" * 60)