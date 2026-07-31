import pandas as pd
import numpy as np
from scipy import stats

np.random.seed(42)  # reproducibility

print("Loading scored data...")
df = pd.read_csv("data/scored_signatures.csv")
print(f"Loaded {df.shape[0]} rows\n")

# Identify gene expression columns (exclude known metadata/score columns)
non_gene_cols = {
    "ModelID", "SequencingID", "ModelConditionID", "IsDefaultEntryForMC",
    "IsDefaultEntryForModel", "OncotreeLineage", "auc", "ic50", "ec50",
    "parpi7_score", "brcaness_score", "cin70_score",
}
gene_cols = [c for c in df.columns if c not in non_gene_cols]
print(f"Found {len(gene_cols)} candidate genes to sample from\n")

# Real signature sizes (based on how many genes were actually usable)
signature_sizes = {
    "parpi7_score": 7,
    "brcaness_score": 60,
    "cin70_score": 70,
}

N_PERM = 500  # number of random gene sets per signature

print("=" * 60)
print("PHASE 5: Permutation testing")
print(f"Comparing real signatures against {N_PERM} random gene sets each")
print("=" * 60)

for sig_col, size in signature_sizes.items():
    # Real signature's correlation with AUC (from Phase 3)
    valid_real = df[[sig_col, "auc"]].dropna()
    real_rho, _ = stats.spearmanr(valid_real[sig_col], valid_real["auc"])

    random_rhos = []
    for i in range(N_PERM):
        random_genes = np.random.choice(gene_cols, size=size, replace=False)
        sub = df[random_genes]
        z = (sub - sub.mean()) / sub.std()
        random_score = z.sum(axis=1)
        valid = pd.DataFrame({"score": random_score, "auc": df["auc"]}).dropna()
        rho, _ = stats.spearmanr(valid["score"], valid["auc"])
        random_rhos.append(rho)

    random_rhos = np.array(random_rhos)
    empirical_p = np.mean(np.abs(random_rhos) >= np.abs(real_rho))

    print(f"\n--- {sig_col} ---")
    print(f"  Real signature rho: {real_rho:.4f}")
    print(f"  Random gene sets (n={N_PERM}): mean rho = {random_rhos.mean():.4f}, std = {random_rhos.std():.4f}")
    print(f"  Empirical p-value: {empirical_p:.4f}")
    print(f"  (fraction of random gene sets with |rho| >= |real signature rho|)")