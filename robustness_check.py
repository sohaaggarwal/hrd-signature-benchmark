import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

print("Loading data...")
df = pd.read_csv("data/scored_signatures.csv")
df["lineage_clean"] = df["OncotreeLineage"].str.replace("[^A-Za-z0-9]", "_", regex=True)

# Re-add BRCA mutant status
mut = pd.read_csv("data/OmicsSomaticMutations.csv", usecols=["ModelID", "HugoSymbol", "LikelyLoF"])
brca_mut = mut[(mut["HugoSymbol"].isin(["BRCA1", "BRCA2"])) & (mut["LikelyLoF"] == True)]
brca_mutant_ids = set(brca_mut["ModelID"].unique())
df["brca_mutant"] = df["ModelID"].isin(brca_mutant_ids).astype(int)

predictors = ["parpi7_score", "brcaness_score", "cin70_score", "brca_mutant"]

print("=" * 60)
print("ROBUSTNESS CHECK: does the conclusion hold using IC50 instead of AUC?")
print("=" * 60)

for metric in ["auc", "ic50"]:
    print(f"\n{'='*20} Using {metric.upper()} {'='*20}")
    for pred in predictors:
        valid = df[[pred, metric, "lineage_clean"]].dropna()
        if len(valid) < 10:
            print(f"  {pred}: too few valid rows, skipping")
            continue
        formula = f"{metric} ~ {pred} + C(lineage_clean)"
        model = smf.ols(formula, data=valid).fit()
        print(f"  {pred}: coef = {model.params[pred]:.5f}, p = {model.pvalues[pred]:.4g}, n = {len(valid)}")