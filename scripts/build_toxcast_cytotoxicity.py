"""
Build a cytotoxicity-only classification dataset from ToxCast's cytotoxicity-associated
"burst" output, in the same normalized schema as scripts/build_tox21_cytotoxicity.py so
the two can be concatenated (see scripts/build_combined_cytotoxicity.py).

Source: EPA's pre-computed chemical-level cytotoxicity burst table (median/lower-bound
AC50 across the ~86 ToxCast assays EPA flags as measuring cytotoxicity/cell-stress rather
than a specific target), from tcpl::tcplCytoPt() run against invitrodb.
  https://epa.figshare.com/articles/dataset/Effects_of_Cell_Stress_and_Cytotoxicity_on_In_Vitro_Assay_Activity_Data_Associated_with_Publication/6062641
  (file: toxcast_cytotox_table.xlsx)

Produces data/toxcast/toxcast_cytotoxicity_classification.csv — one row per chemical.
Unlike Tox21 (individual qHTS AIDs per assay/cell line/timepoint), this EPA output is
already aggregated across all burst assays, so there is no equivalent "long" file here.

Usage:
    python scripts/build_toxcast_cytotoxicity.py [--out-dir data/toxcast] [--raw-dir /tmp/toxcast_raw]

Requires: requests, pandas, openpyxl
"""
import argparse
import os
import time
from io import StringIO

import pandas as pd
import requests

CYTOTOX_XLSX_URL = "https://ndownloader.figshare.com/files/39823492"
CASRN_XREF_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{casrn}/cids/JSON"
PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/"
    "ConnectivitySMILES,SMILES,MolecularFormula,MolecularWeight/CSV"
)


def download_cytotox_table(raw_dir: str) -> pd.DataFrame:
    os.makedirs(raw_dir, exist_ok=True)
    path = os.path.join(raw_dir, "toxcast_cytotox_table.xlsx")
    if not os.path.exists(path):
        r = requests.get(CYTOTOX_XLSX_URL, timeout=120)
        r.raise_for_status()
        with open(path, "wb") as f:
            f.write(r.content)
    return pd.read_excel(path)


def resolve_casrn_to_cid(casrns: list[str], max_workers: int = 6) -> dict:
    """One PubChem xref lookup per CASRN (no batch endpoint for RegistryID xref),
    parallelized with a small thread pool to keep this from taking ~30+ minutes serial."""
    from concurrent.futures import ThreadPoolExecutor, as_completed

    session = requests.Session()

    def lookup(casrn):
        for attempt in range(3):
            try:
                r = session.get(CASRN_XREF_URL.format(casrn=casrn), timeout=20)
                if r.status_code == 200:
                    cids = r.json().get("IdentifierList", {}).get("CID", [])
                    return casrn, (cids[0] if cids else None)
                elif r.status_code == 404:
                    return casrn, None
            except Exception:
                pass
            time.sleep(1)
        return casrn, None

    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        for fut in as_completed([ex.submit(lookup, c) for c in casrns]):
            casrn, cid = fut.result()
            results[casrn] = cid
    return results


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="data/toxcast")
    ap.add_argument("--raw-dir", default="/tmp/toxcast_raw")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    df = download_cytotox_table(args.raw_dir)

    casrns = df["casrn"].dropna().unique().tolist()
    casrn_to_cid = resolve_casrn_to_cid(casrns)
    df["cid"] = df["casrn"].map(casrn_to_cid)
    df_ok = df.dropna(subset=["cid"]).copy()
    df_ok["cid"] = df_ok["cid"].astype(int)
    print(f"Resolved {len(df_ok)}/{len(df)} chemicals to a PubChem CID")

    smiles = fetch_smiles(sorted(df_ok["cid"].unique().tolist()))
    merged = df_ok.merge(smiles, on="cid", how="left")

    # nhit > 0: chemical triggered the burst in >=1 of the ~86 burst-flagged assays.
    # No "inconclusive" category exists in this pre-aggregated EPA output (unlike Tox21's
    # per-assay qHTS curve calls).
    merged["cytotoxic_label"] = (merged["nhit"] > 0).astype(int)
    # Normalized columns shared across all sources in this collection (see data/README.md /
    # scripts/build_combined_cytotoxicity.py) alongside the ToxCast-specific originals.
    merged["potency_value"] = merged["cytotox_median_um"]
    merged["n_assays_active"] = merged["nhit"]
    merged["n_assays_tested"] = merged["ntested"]

    out = merged[
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
            "cytotox_median_um",
            "cytotox_lower_bound_um",
            "nhit",
            "ntested",
            "dtxsid",
            "casrn",
            "name",
        ]
    ].copy()
    out.insert(0, "source", "ToxCast")
    out.insert(1, "biological_context", "off-target-safety")
    out["potency_metric"] = "burst_median_AC50"
    out["potency_unit"] = "uM"
    out["assay_description"] = (
        "ToxCast cytotoxicity-associated burst (median/lower-bound AC50 across ~86 "
        "burst-flagged assays, invitrodb v3.4)"
    )
    out = out.sort_values("cid").reset_index(drop=True)
    out.to_csv(os.path.join(args.out_dir, "toxcast_cytotoxicity_classification.csv"), index=False)

    print(f"Final shape: {out.shape}")
    print(out["cytotoxic_label"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
