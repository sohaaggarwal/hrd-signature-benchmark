
import pandas as pd

print("Loading expression data...")
expr = pd.read_csv("data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")

print("Loading cell line metadata...")
model = pd.read_csv("data/Model.csv")

print("Loading olaparib drug response data...")
drug = pd.read_csv("data/prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv")

expr = expr.rename(columns=lambda c: c.split(" (")[0])

olaparib = drug[drug["name"].str.lower() == "olaparib"].copy()
olaparib = olaparib[["depmap_id", "auc", "ic50", "ec50"]]
olaparib = olaparib.rename(columns={"depmap_id": "ModelID"})

lineage_info = model[["ModelID", "OncotreeLineage"]]

merged = expr.merge(lineage_info, on="ModelID", how="inner")
merged = merged.merge(olaparib, on="ModelID", how="inner")

print(f"\nFinal merged table: {merged.shape[0]} cell lines, {merged.shape[1]} columns")

target_lineages = ["Ovary/Fallopian Tube", "Breast", "Pancreas"]
merged_filtered = merged[merged["OncotreeLineage"].isin(target_lineages)]

print(f"After filtering to target lineages: {merged_filtered.shape[0]} cell lines")
print(merged_filtered["OncotreeLineage"].value_counts())

merged_filtered.to_csv("data/merged_analysis_ready.csv", index=False)
print("\nSaved to data/merged_analysis_ready.csv")