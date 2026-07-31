import pandas as pd
import statsmodels.formula.api as smf

print("Loading scored data...")
df = pd.read_csv("data/scored_signatures.csv")
print(f"Loaded {df.shape[0]} rows\n")

signature_cols = ["parpi7_score", "brcaness_score", "cin70_score"]

print("=" * 60)
print("PHASE 4: Lineage-adjusted regression")
print("Model: auc ~ signature_score + lineage")
print("This tests whether the signature predicts drug response")
print("AFTER statistically accounting for cancer lineage.")
print("=" * 60)

# Clean up lineage names for use as a formula variable (no slashes/spaces issues)
df["lineage_clean"] = df["OncotreeLineage"].str.replace("[^A-Za-z0-9]", "_", regex=True)

for sig in signature_cols:
    valid = df[[sig, "auc", "lineage_clean"]].dropna()
    formula = f"auc ~ {sig} + C(lineage_clean)"
    model = smf.ols(formula, data=valid).fit()

    print(f"\n--- {sig} ---")
    print(f"  Coefficient on {sig}: {model.params[sig]:.5f}")
    print(f"  P-value on {sig}:     {model.pvalues[sig]:.4g}")
    print(f"  R-squared (full model): {model.rsquared:.4f}")
    print(f"  N = {len(valid)}")