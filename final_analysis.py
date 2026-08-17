import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from scipy import stats

np.random.seed(42)

print("=" * 70)
print("STEP 1: Rebuild merged data with deduplication + both drugs")
print("=" * 70)

expr = pd.read_csv("data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
expr = expr.rename(columns=lambda c: c.split(" (")[0])
model = pd.read_csv("data/Model.csv")
drug = pd.read_csv("data/prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv")

lineage_info = model[["ModelID", "OncotreeLineage"]]
target_lineages = ["Ovary/Fallopian Tube", "Breast", "Pancreas"]

def build_drug_table(drug_name):
    d = drug[drug["name"].str.lower() == drug_name].copy()
    d = d[["depmap_id", "auc", "ic50"]].rename(columns={"depmap_id": "ModelID"})
    # DEDUPLICATE: average across replicates per cell line
    before = len(d)
    d = d.groupby("ModelID", as_index=False).mean(numeric_only=True)
    print(f"  {drug_name}: {before} rows -> {len(d)} unique cell lines after dedup")
    return d

for drug_name in ["olaparib", "rucaparib"]:
    print(f"\n{'='*70}")
    print(f"ANALYSIS FOR: {drug_name.upper()}")
    print(f"{'='*70}")

    d = build_drug_table(drug_name)
    df = expr.merge(lineage_info, on="ModelID", how="inner").merge(d, on="ModelID", how="inner")
    df = df[df["OncotreeLineage"].isin(target_lineages)].copy()
    print(f"  Final analysis set: {len(df)} cell lines")
    print(df["OncotreeLineage"].value_counts().to_string())

    # score signatures
    PARPI7 = {"sens": ["CHEK2", "MAPKAPK2"], "res": ["BRCA1", "MRE11", "NBN", "TDG", "XPA"]}
    BRCANESS_W = {
        "DAD1":0.0997,"RAD21":0.1743,"LDHA":0.0165,"SPARC":0.1571,"SKP1":0.137,"PPP1CC":0.0249,
        "RAN":0.1043,"USP9X":-0.0066,"ECHS1":-0.0495,"GNAI3":0.1404,"HDAC1":-0.0272,"LGMN":-0.0683,
        "CYR61":-0.0454,"SH3BGRL":0.0659,"MGST3":-0.0162,"RCN2":0.2087,"SMC1A":0.0313,"BST2":0.2118,
        "ENG":0.058,"SRGN":-0.1577,"GBP1":0.0778,"UNC119B":0.139,"RGS1":-0.2573,"CDC2":0.1101,
        "SAT1":-0.0386,"GGH":-0.0169,"GUCY1B3":0.0248,"WFDC2":0.0719,"NMI":0.0873,"PRAME":0.0267,
        "CCL4":-0.1545,"MGST2":-0.0045,"MMP7":0,"SNCA":-0.1605,"VTN":-0.0399,"ALPP":-0.0274,
        "MTAP":-0.0511,"IDUA":0.097,"SERPINF2":-0.0053,"WAS":-0.0191,"CD1D":-0.1431,"GFI1":-0.0318,
        "P11":0.0091,"SEMA3F":0.0805,"TNF":-0.024,"ROS1":0.1667,"MADCAM1":-0.0322,"PDIA4":0.1564,
        "HMGN2":-0.0247,"HLA-B":0.0956,"TM9SF1":0.1125,"CCDC93":-0.0865,"APEX1":0.0212,
        "VEGFA":0.0612,"POSTN":0.2142,"PSTPIP1":-0.041,"PMS1":-0.0579,"HLA-A":0.0202,
        "PCTP":-0.0536,"SEH1L":0.1127,
    }
    CIN70 = ["TPX2","PRC1","FOXM1","CDC2","C20orf24","TGIF2","MCM2","H2AFZ","TOP2A","PCNA","UBE2C",
        "MELK","TRIP13","CNAP1","MCM7","RNASEH2A","RAD51AP1","KIF20A","CDC45L","MAD2L1","ESPL1",
        "CCNB2","FEN1","TTK","CCT5","RFC4","ATAD2","CKAP5","NUP205","CDC20","CKS2","RRM2","ELAVL1",
        "CCNB1","RRM1","AURKB","MSH6","EZH2","CTPS","DKC1","OIP5","CDCA8","PTTG1","CEP55","H2AFX",
        "CMAS","NCAPH","MCM10","LSM4","NCAPH2","ASF1B","ZWINT","PBK","FLJ10036","CDCA3","ECT2",
        "CDC6","UNG","MTCH2","RAD21","ACTL6A","GPI","MGC13096","SFRS2","HDGF","NXT1","NEK2",
        "DHCR7","AURKA","NDUFAB1","KIAA0286","KIF4A"]

    def z(cols):
        sub = df[[c for c in cols if c in df.columns]]
        return (sub - sub.mean()) / sub.std()

    w = {g: 1 for g in PARPI7["sens"]}
    w.update({g: -1 for g in PARPI7["res"]})
    zp = z(PARPI7["sens"] + PARPI7["res"])
    df["parpi7_score"] = sum(zp[g] * w[g] for g in zp.columns)

    zb = z(list(BRCANESS_W.keys()))
    df["brcaness_score"] = sum(zb[g] * BRCANESS_W[g] for g in zb.columns)

    df["cin70_score"] = z(CIN70).sum(axis=1)

    df["lineage_clean"] = df["OncotreeLineage"].str.replace("[^A-Za-z0-9]", "_", regex=True)

    print(f"\n  --- Lineage-adjusted regression (outcome = AUC) ---")
    for sig in ["parpi7_score", "brcaness_score", "cin70_score"]:
        valid = df[[sig, "auc", "lineage_clean"]].dropna()
        m = smf.ols(f"auc ~ {sig} + C(lineage_clean)", data=valid).fit()
        print(f"    {sig}: coef={m.params[sig]:.5f}, p={m.pvalues[sig]:.4g}, n={len(valid)}")

    df.to_csv(f"data/final_{drug_name}.csv", index=False)
    print(f"\n  Saved to data/final_{drug_name}.csv")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)