import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

df = pd.read_csv("data/final_olaparib.csv")
print(f"Loaded {len(df)} cell lines")

sig_cols = ["parpi7_score", "brcaness_score", "cin70_score"]
sig_labels = {"parpi7_score": "PARPi-7", "brcaness_score": "BRCAness", "cin70_score": "CIN70"}

# ---------- FIGURE 1: scatter plots, signature score vs AUC ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

for ax, sig in zip(axes, sig_cols):
    for lineage in df["OncotreeLineage"].unique():
        sub = df[df["OncotreeLineage"] == lineage]
        ax.scatter(sub[sig], sub["auc"], label=lineage, alpha=0.7, s=40)
    ax.set_xlabel(f"{sig_labels[sig]} score")
    ax.set_ylabel("Olaparib AUC")
    ax.set_title(f"{sig_labels[sig]} vs olaparib response")

axes[0].legend(fontsize=8, loc="best")
plt.suptitle("HRD signature scores show no relationship with olaparib sensitivity", y=1.02)
plt.tight_layout()
plt.savefig("results_figure1_scatter.png", dpi=200, bbox_inches="tight")
print("Saved results_figure1_scatter.png")
plt.close()

# ---------- FIGURE 2: correlation heatmap between signatures ----------
corr = df[sig_cols].corr(method="spearman")
corr.index = [sig_labels[c] for c in corr.index]
corr.columns = [sig_labels[c] for c in corr.columns]

plt.figure(figsize=(6, 5))
sns.heatmap(corr, annot=True, fmt=".3f", cmap="RdBu_r", center=0,
            vmin=-1, vmax=1, square=True, cbar_kws={"label": "Spearman rho"})
plt.title("Do the three signatures agree with each other?")
plt.tight_layout()
plt.savefig("results_figure2_heatmap.png", dpi=200, bbox_inches="tight")
print("Saved results_figure2_heatmap.png")
plt.close()

print("\nDone — both figures saved to your project folder.")