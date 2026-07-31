import pandas as pd

print("Loading mutations file...")
mut = pd.read_csv("data/OmicsSomaticMutations.csv", nrows=5)
print("\nColumns found:")
for col in mut.columns:
    print(f"  {col}")
print("\nFirst few rows:")
print(mut)