# ToxCast Cytotoxicity Dataset

Cytotoxicity-only slice of ToxCast, extracted from EPA's pre-computed **cytotoxicity-associated
"burst"** table — chemicals showing activation of large numbers of assays over a narrow
concentration range in which general cell stress/cytotoxicity is also seen. Source:
[ToxCast Cytotoxicity-Associated Burst, EPA Figshare](https://epa.figshare.com/articles/dataset/Effects_of_Cell_Stress_and_Cytotoxicity_on_In_Vitro_Assay_Activity_Data_Associated_with_Publication/6062641)
(`toxcast_cytotox_table.xlsx`, from `tcpl::tcplCytoPt()` run against invitrodb v3.4).

Built with [`scripts/build_toxcast_cytotoxicity.py`](../../scripts/build_toxcast_cytotoxicity.py).
Built 2026-08-31 from the live EPA Figshare file — rerun the script to refresh.

## What this is, structurally

Unlike Tox21's dataset here, this is **already a chemical-level aggregate**, not a set of
individual assay AIDs: EPA computed one summary cytotoxicity potency value per chemical
(`cytotox_median_um`, `cytotox_lower_bound_um`) directly from the ~86 ToxCast assays flagged as
measuring cytotoxicity/cell-stress rather than a specific molecular target (`burst_assay == 1` in
ToxCast's `Assay_Summary` files). **There is no equivalent "long" file for ToxCast in this
collection** — the underlying per-assay burst data lives in EPA's full invitrodb release, which
is far larger and not chemical-identifier-clean (DTXSID-keyed, no SMILES) without the same
CASRN→PubChem resolution step done here.

## Files

### `toxcast_cytotoxicity_classification.csv` — one row per chemical, ML-ready

7,805 rows (of 9,116 chemicals in EPA's source table; 1,311 could not be resolved to a PubChem
CID via CASRN and were dropped — see Limitations).

| Column | Description |
|---|---|
| `source` | Always `ToxCast` |
| `cid` | PubChem Compound ID (resolved from CASRN via PubChem xref lookup — not in the EPA source file) |
| `canonical_smiles`, `isomeric_smiles`, `molecular_formula`, `molecular_weight` | From PubChem, same convention as the Tox21 files |
| `cytotoxic_label` | Binary target: `1` if `nhit > 0` (chemical triggered the burst in ≥1 of the ~86 assays), else `0`. No inconclusive category — see below |
| `potency_um` | Alias of `cytotox_median_um`, added for schema alignment with Tox21 (see Combining) |
| `n_assays_active` | Alias of `nhit` |
| `n_assays_tested` | Alias of `ntested` |
| `cytotox_median_um` | EPA's median AC50 (µM) across the chemical's burst-flagged hits |
| `cytotox_lower_bound_um` | EPA's conservative lower-bound cytotoxicity estimate (µM) |
| `nhit` / `ntested` | Number of the ~86 burst assays the chemical was active in / actually tested in |
| `dtxsid`, `casrn`, `name` | Original EPA/DSSTox chemical identifiers, kept for traceability |
| `assay_description` | Free-text summary of the source computation |

**Label derivation:** `cytotoxic_label = 1` iff `nhit > 0`. This is coarser than Tox21's label —
EPA's source table doesn't retain an "inconclusive" state at the chemical level, so every chemical
gets a 0/1 call. If you need a stricter definition (e.g. minimum number of hits, or a potency
cutoff on `cytotox_median_um`), rebuild from `nhit`/`ntested`/`cytotox_median_um` directly.

**Class balance:** 4,184 active / 3,621 inactive (of the resolved 7,805) — close to even, unlike
Tox21's ~20% active rate. This reflects a structurally different labeling rule (any hit among 86
assays vs. a dedicated 4-assay viability panel), not necessarily a difference in true cytotoxicity
prevalence between the two chemical libraries — keep that in mind before pooling label statistics.

## Combining with the Tox21 dataset

Both classification files share the same **normalized column set** — `source`, `cid`,
`canonical_smiles`, `isomeric_smiles`, `molecular_formula`, `molecular_weight`,
`cytotoxic_label`, `potency_um`, `n_assays_active`, `n_assays_tested`, `assay_description` — so
they concatenate directly:

```python
import pandas as pd
tox21 = pd.read_csv("data/tox21/tox21_cytotoxicity_classification.csv")
toxcast = pd.read_csv("data/toxcast/toxcast_cytotoxicity_classification.csv")
combined = pd.concat([tox21[COMMON_COLS], toxcast[COMMON_COLS]], ignore_index=True)
```

or just use the pre-built [`data/combined/combined_cytotoxicity_classification.csv`](../combined/)
(built by [`scripts/build_combined_cytotoxicity.py`](../../scripts/build_combined_cytotoxicity.py)).
**6,286 compounds (CIDs) appear in both sources** — the combined file keeps them as separate rows
(one per source) rather than merging their labels, since the two panels measure cytotoxicity
differently (dedicated real-time viability assay vs. cross-target burst pattern) and shouldn't be
silently averaged. Dedupe/reconcile explicitly if your use case needs one row per compound.

## Known limitations

- **1,311 chemicals dropped** because their CASRN didn't resolve to a PubChem CID (obsolete/
  non-standard CAS numbers, mixtures, or substances not in PubChem). The full 9,116-chemical EPA
  table (DTXSID-keyed, no SMILES) is cached at build time if you need those too.
- `potency_um`/`cytotox_median_um` is set to a default (1000 µM in the source data) for chemicals
  with no burst hit — this is EPA's convention for "no cytotoxicity detected up to the tested
  range," not a measured value. Filter on `cytotoxic_label`/`nhit` before using `potency_um` as a
  continuous target.
- Structures are PubChem's, resolved by CASRN — not independently re-standardized.
