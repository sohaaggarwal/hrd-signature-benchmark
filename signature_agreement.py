import pandas as pd
from scipy import stats

print("Loading scored data...")
df = pd.read_csv("data/scored_signatures.csv")
print(f"Loaded {df.shape[0]} rows\n")

signature_cols = ["parpi7_score", "brcaness_score", "cin70_score"]

print("=" * 60)
print("DO THE SIGNATURES AGREE WITH EACH OTHER?")
print("If they're all measuring the same underlying HRD biology,")
print("they should correlate with each other -- even if they don't")
print("predict drug response.")
print("=" * 60)

pairs = [
    ("parpi7_score", "brcaness_score"),
    ("parpi7_score", "cin70_score"),
    ("brcaness_score", "cin70_score"),
]

for sig_a, sig_b in pairs:
    valid = df[[sig_a, sig_b]].dropna()
    rho, pval = stats.spearmanr(valid[sig_a], valid[sig_b])
    print(f"\n{sig_a} vs {sig_b}:")
    print(f"  Spearman rho = {rho:.4f}, p = {pval:.4g}, n = {len(valid)}")

print("\n" + "=" * 60)
print("Full correlation matrix (for reference)")
print("=" * 60)
corr_matrix = df[signature_cols].corr(method="spearman")
print(corr_matrix)