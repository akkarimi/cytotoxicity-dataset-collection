# Tox21 Cytotoxicity Dataset

Cytotoxicity-only slice of Tox21, extracted from PubChem AID 1224867–1224890 — the 24-assay
real-time HEK293/HepG2 viability and cytotoxicity panel run specifically as a cytotoxicity
counter-screen for the rest of the Tox21 10K program (separate from the ~70 nuclear-receptor
and stress-pathway assays). Source study: [Huang et al., "Real-time cell toxicity profiling of
Tox21 10K compounds reveals cytotoxicity dependent toxicity pathway linkage" (PLOS ONE,
2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5439695/).

Built with [`scripts/build_tox21_cytotoxicity.py`](../../scripts/build_tox21_cytotoxicity.py).
Built 2026-08-31 from live PubChem data — rerun the script to refresh.

## Assay design

- **2 technologies:** RealTime-Glo MT Cell Viability Assay (metabolic activity, Promega) and
  CellTox Green Cytotoxicity Assay (membrane integrity / cell death, Promega)
- **2 cell lines:** HEK293 (embryonic kidney), HepG2 (hepatocellular carcinoma)
- **6 timepoints:** 0, 8, 16, 24, 32, 40 hours post-exposure
- 2 × 2 × 6 = 24 PubChem AIDs, one per combination

## Files

### `tox21_cytotoxicity_long.csv` — raw, one row per compound × cell line × technology × timepoint

228,576 rows, 7,671 unique compounds (CIDs).

| Column | Description |
|---|---|
| `source` | Always `Tox21` |
| `cid` | PubChem Compound ID |
| `canonical_smiles` | Connectivity-only SMILES (PubChem `ConnectivitySMILES`) |
| `isomeric_smiles` | Full-stereochemistry SMILES (PubChem `SMILES`) |
| `molecular_formula`, `molecular_weight` | From PubChem |
| `cell_line` | `HEK293` or `HepG2` |
| `technology` | `RealTime-Glo MT Cell Viability` or `CellTox Green Cytotoxicity` |
| `timepoint_h` | Exposure time in hours (0/8/16/24/32/40) |
| `hit_call` | `1` = active (cytotoxic), `0` = inactive, blank = inconclusive (excluded from binary call, not imputed) |
| `ac50_um` | Potency (µM) where a curve was fit; blank for inactive/inconclusive compounds |
| `sid` | PubChem Substance ID for that AID/compound |
| `aid` | Source PubChem AID |

Use this file if you want the full time-course / cross-technology signal, or want to define
your own aggregation rule.

### `tox21_cytotoxicity_classification.csv` — one row per compound, ML-ready

7,671 rows (one per unique compound).

| Column | Description |
|---|---|
| `source` | Always `Tox21` |
| `cid`, `canonical_smiles`, `isomeric_smiles`, `molecular_formula`, `molecular_weight` | As above |
| `cytotoxic_label` | Binary classification target — see derivation below. `NaN` where inconclusive across all 4 assays |
| `min_ac50_um_40h` | Lowest (most potent) AC50 µM among the 40h assays the compound was active in; blank if never active |
| `n_40h_assays_active` / `n_40h_assays_tested` | Out of the 4 assay/cell-line combinations at 40h, how many were active / had a conclusive call |
| `assay_description` | Free-text summary of which assays fed the label |

**Label derivation:** uses the 40-hour timepoint only (longest exposure, most sensitive/consistent
readout per the source study), across all 4 technology × cell-line combinations at that timepoint.
`cytotoxic_label = 1` if active in **any** of the 4; `= 0` if inactive in **all** 4 with no
inconclusive results; `NaN` otherwise. This is a design choice, not the only valid one — the long
format lets you re-aggregate differently (e.g. require agreement across cell lines, use a different
timepoint, or keep technologies separate as independent labels for a multi-task setup).

**Class balance:** 4,788 inactive (62%) / 1,512 active (20%) / 1,371 inconclusive (18%).

## Known limitations

- SMILES are as deposited in PubChem for each CID — not independently re-standardized (salts,
  stereochemistry, mixtures not stripped/canonicalized beyond what PubChem provides).
- "Inconclusive" qHTS calls are dropped from the binary label rather than imputed either way;
  decide explicitly whether to exclude, or treat as a third class, before training on this.
- 40h aggregation is a single reasonable default (see above) — re-run `build_classification_summary`
  in the build script with a different rule if your use case needs one.
