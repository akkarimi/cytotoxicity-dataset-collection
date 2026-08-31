"""
Build a cytotoxicity-only classification dataset from Tox21's real-time
HEK293/HepG2 viability and cytotoxicity assay panel.

Source: PubChem AID 1224867-1224890 (24 assays = 2 assay technologies x
2 cell lines x 6 exposure timepoints), from the Tox21 qHTS real-time
cytotoxicity screen (Huang et al., PLOS ONE 2017,
https://pmc.ncbi.nlm.nih.gov/articles/PMC5439695/).

Produces two output files (see data/tox21/README.md for the schema):
  - tox21_cytotoxicity_long.csv          one row per compound x cell line
                                          x technology x timepoint
  - tox21_cytotoxicity_classification.csv one row per compound, aggregated
                                          to a single binary cytotox label

Usage:
    python scripts/build_tox21_cytotoxicity.py [--out-dir data/tox21] [--raw-dir /tmp/tox21_raw]

Requires: requests, pandas
"""
import argparse
import glob
import os
import re
import time
from io import StringIO

import numpy as np
import pandas as pd
import requests

CYTOTOX_AIDS = list(range(1224867, 1224891))  # 24 AIDs, confirmed boundary at 1224890
CONCISE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/concise/CSV"
PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/"
    "ConnectivitySMILES,SMILES,MolecularFormula,MolecularWeight/CSV"
)


def download_assay_data(raw_dir: str) -> pd.DataFrame:
    os.makedirs(raw_dir, exist_ok=True)
    for aid in CYTOTOX_AIDS:
        path = os.path.join(raw_dir, f"aid_{aid}_concise.csv")
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            r = requests.get(CONCISE_URL.format(aid=aid), timeout=60)
            r.raise_for_status()
            with open(path, "w") as f:
                f.write(r.text)
            time.sleep(0.3)

    frames = [pd.read_csv(f, dtype=str) for f in sorted(glob.glob(os.path.join(raw_dir, "aid_*_concise.csv")))]
    full = pd.concat(frames, ignore_index=True)
    full.columns = [c.strip() for c in full.columns]
    return full


def parse_assay_name(name: str):
    tech = "RealTime-Glo MT Cell Viability" if "RealTime-Glo" in name else "CellTox Green Cytotoxicity"
    cell = "HEK293" if "HEK293" in name else "HepG2"
    m = re.search(r"(\d+)\s*hour", name)
    tp = int(m.group(1)) if m else None
    return tech, cell, tp


def build_long_raw(full: pd.DataFrame) -> pd.DataFrame:
    parsed = full["Assay Name"].apply(parse_assay_name)
    full = full.copy()
    full["technology"] = parsed.apply(lambda x: x[0])
    full["cell_line"] = parsed.apply(lambda x: x[1])
    full["timepoint_h"] = parsed.apply(lambda x: x[2])
    full["cid"] = pd.to_numeric(full["CID"], errors="coerce")
    full["sid"] = pd.to_numeric(full["SID"], errors="coerce")
    full["ac50_um"] = pd.to_numeric(full["Activity Value [uM]"], errors="coerce")
    # "Inconclusive" outcomes are left as NaN rather than forced into a binary class
    full["hit_call"] = full["Activity Outcome"].map({"Active": 1, "Inactive": 0})

    out = full[["cid", "sid", "cell_line", "technology", "timepoint_h", "hit_call", "ac50_um", "AID"]].rename(
        columns={"AID": "aid"}
    )
    out = out.dropna(subset=["cid"])
    out["cid"] = out["cid"].astype(int)
    out["sid"] = out["sid"].astype(int)
    return out


def fetch_smiles(cids: list[int]) -> pd.DataFrame:
    chunks = [cids[i : i + 200] for i in range(0, len(cids), 200)]
    results = []
    for chunk in chunks:
        for attempt in range(3):
            r = requests.post(PROPERTY_URL, data={"cid": ",".join(str(c) for c in chunk)}, timeout=60)
            if r.status_code == 200:
                break
            time.sleep(2)
        else:
            raise RuntimeError(f"Failed to fetch SMILES for chunk starting {chunk[0]}")
        results.append(pd.read_csv(StringIO(r.text), dtype=str))
        time.sleep(0.3)
    smiles = pd.concat(results, ignore_index=True)
    return smiles.rename(
        columns={
            "CID": "cid",
            "ConnectivitySMILES": "canonical_smiles",
            "SMILES": "isomeric_smiles",
            "MolecularFormula": "molecular_formula",
            "MolecularWeight": "molecular_weight",
        }
    )


def build_classification_summary(long_raw: pd.DataFrame) -> pd.DataFrame:
    """Collapse to one row per compound using the 40h (longest-exposure) readouts
    across all 4 assay/cell-line combinations: active if active in ANY of them,
    inactive if inactive in all of them, else inconclusive (NaN label)."""
    primary = long_raw[long_raw["timepoint_h"] == 40].copy()

    def agg(g):
        any_active = (g["hit_call"] == 1).any()
        all_inactive = (g["hit_call"] == 0).all() and g["hit_call"].notna().any()
        label = 1 if any_active else (0 if all_inactive else np.nan)
        ac50 = g.loc[g["hit_call"] == 1, "ac50_um"].min()
        return pd.Series(
            {
                "cytotoxic_label": label,
                "min_ac50_um_40h": ac50,
                "n_40h_assays_active": (g["hit_call"] == 1).sum(),
                "n_40h_assays_tested": g["hit_call"].notna().sum(),
            }
        )

    return primary.groupby("cid").apply(agg, include_groups=False).reset_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="data/tox21")
    ap.add_argument("--raw-dir", default="/tmp/tox21_raw")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    full = download_assay_data(args.raw_dir)
    long_raw = build_long_raw(full)
    smiles = fetch_smiles(sorted(long_raw["cid"].unique().tolist()))

    long_df = long_raw.merge(smiles, on="cid", how="left")
    long_df = long_df[
        [
            "cid",
            "canonical_smiles",
            "isomeric_smiles",
            "molecular_formula",
            "molecular_weight",
            "cell_line",
            "technology",
            "timepoint_h",
            "hit_call",
            "ac50_um",
            "sid",
            "aid",
        ]
    ].copy()
    long_df.insert(0, "source", "Tox21")
    long_df = long_df.sort_values(["cid", "cell_line", "technology", "timepoint_h"]).reset_index(drop=True)
    long_df.to_csv(os.path.join(args.out_dir, "tox21_cytotoxicity_long.csv"), index=False)

    summary = build_classification_summary(long_raw).merge(smiles, on="cid", how="left")
    # Normalized columns shared across all sources in this collection (see data/README.md /
    # scripts/build_combined_cytotoxicity.py) alongside the Tox21-specific originals.
    summary["potency_value"] = summary["min_ac50_um_40h"]
    summary["n_assays_active"] = summary["n_40h_assays_active"]
    summary["n_assays_tested"] = summary["n_40h_assays_tested"]
    summary = summary[
        [
            "cid",
            "canonical_smiles",
            "isomeric_smiles",
            "molecular_formula",
            "molecular_weight",
            "cytotoxic_label",
            "potency_value",
            "n_assays_active",
            "n_assays_tested",
            "min_ac50_um_40h",
            "n_40h_assays_active",
            "n_40h_assays_tested",
        ]
    ].copy()
    summary.insert(0, "source", "Tox21")
    summary.insert(1, "biological_context", "off-target-safety")
    summary["potency_metric"] = "AC50"
    summary["potency_unit"] = "uM"
    summary["assay_description"] = (
        "RealTime-Glo MT Cell Viability + CellTox Green Cytotoxicity, HEK293 & HepG2, 40h exposure (qHTS)"
    )
    summary = summary.sort_values("cid").reset_index(drop=True)
    summary.to_csv(os.path.join(args.out_dir, "tox21_cytotoxicity_classification.csv"), index=False)

    print(f"Long format: {long_df.shape}")
    print(f"Classification format: {summary.shape}")
    print(summary["cytotoxic_label"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
