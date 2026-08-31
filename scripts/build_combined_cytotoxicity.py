"""
Concatenate the per-source cytotoxicity classification files on their shared normalized
columns. Run scripts/build_tox21_cytotoxicity.py and scripts/build_toxcast_cytotoxicity.py
first.

Usage:
    python scripts/build_combined_cytotoxicity.py
"""
import pandas as pd

COMMON_COLS = [
    "source",
    "biological_context",
    "cid",
    "canonical_smiles",
    "isomeric_smiles",
    "molecular_formula",
    "molecular_weight",
    "cytotoxic_label",
    "potency_value",
    "potency_metric",
    "potency_unit",
    "n_assays_active",
    "n_assays_tested",
    "assay_description",
]

SOURCE_FILES = {
    "Tox21": "data/tox21/tox21_cytotoxicity_classification.csv",
    "ToxCast": "data/toxcast/toxcast_cytotoxicity_classification.csv",
    "GDSC": "data/gdsc/gdsc_cytotoxicity_classification.csv",
}

frames = {}
for name, path in SOURCE_FILES.items():
    df = pd.read_csv(path)
    frames[name] = df
    print(f"{name}: {len(df)} rows")

combined = pd.concat([f[COMMON_COLS] for f in frames.values()], ignore_index=True)
combined.to_csv("data/combined/combined_cytotoxicity_classification.csv", index=False)
print(f"Combined: {len(combined)} rows")

names = list(frames.keys())
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        overlap = set(frames[a]["cid"]) & set(frames[b]["cid"])
        print(f"Compounds (CIDs) present in both {a} and {b}: {len(overlap)}")
