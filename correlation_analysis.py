import pandas as pd
from scipy import stats

print("Loading scored data...")
df = pd.read_csv("data/scored_signatures.csv")
print(f"Loaded {df.shape[0]} rows\n")

signature_cols = ["parpi7_score", "brcaness_score", "cin70_score"]

print("=" * 60)
print("PHASE 3: Basic correlation (all cell lines, no lineage control)")
print("=" * 60)

for sig in signature_cols:
    valid = df[[sig, "auc"]].dropna()
    rho, pval = stats.spearmanr(valid[sig], valid["auc"])
    print(f"\n{sig}: rho = {rho:.4f}, p = {pval:.4g}, n = {len(valid)}")

print("\n" + "=" * 60)
print("Same correlation, broken down by cancer lineage")
print("=" * 60)

for lineage in df["OncotreeLineage"].unique():
    sub = df[df["OncotreeLineage"] == lineage]
    print(f"\n--- {lineage} (n={len(sub)}) ---")
    for sig in signature_cols:
        valid = sub[[sig, "auc"]].dropna()
        if len(valid) < 5:
            print(f"  {sig}: too few samples ({len(valid)}), skipping")
            continue
        rho, pval = stats.spearmanr(valid[sig], valid["auc"])
        print(f"  {sig}: rho = {rho:.4f}, p = {pval:.4g}, n = {len(valid)}")