# GDSC Cytotoxicity Dataset

Cytotoxicity/growth-inhibition data from the **Genomics of Drug Sensitivity in Cancer** project —
fitted dose-response (IC50, AUC) for drug–cancer-cell-line pairs across the combined GDSC1 + GDSC2
screens. Source: [cancerrxgene.org bulk download](https://www.cancerrxgene.org/downloads/bulk_download)
(`GDSC1_fitted_dose_response_24Jul22.csv`, `GDSC2_fitted_dose_response_24Jul22.csv`,
`screened_compounds_rel_8.4.csv`).

Built with [`scripts/build_gdsc_cytotoxicity.py`](../../scripts/build_gdsc_cytotoxicity.py).
Built 2026-08-31 from the live Sanger FTP files — rerun the script to refresh.

**This is a different biological context from Tox21/ToxCast** — see the top-level
[`data/README.md`](../README.md) and [`biological_context` column below]. GDSC measures
drug-induced killing of **cancer cell lines** (an efficacy signal — "the drug is working"), not
off-target toxicity in normal cells (a safety signal). Every row here is tagged
`biological_context = "on-target-cancer-efficacy"`.

## Files

### `gdsc_cytotoxicity_long.csv` — one row per compound × cell line × GDSC dataset

One row per drug-cell line pair actually screened. Columns include `gdsc_dataset` (GDSC1/GDSC2 —
different, non-overlapping drug panels and screening technologies), `cell_line`,
`sanger_model_id`/`cosmic_id` (cell line identifiers), `tissue` (TCGA cancer type descriptor),
`putative_target`/`pathway_name` (drug's nominal molecular target, from GDSC annotation — most
compounds here are targeted oncology drugs, not general chemicals), `ic50_um` (= `exp(LN_IC50)`,
GDSC reports natural-log IC50), `auc` (0–1, lower = more sensitive), `z_score` (GDSC's own
per-drug standardized sensitivity score, relative to that drug's own distribution across lines —
not comparable across different drugs).

### `gdsc_cytotoxicity_classification.csv` — one row per compound, ML-ready

**Label derivation — potency + coverage threshold** (a deliberate design choice, not a
GDSC-native "hit call" — GDSC's raw output is continuous, unlike ToxCast's `nhit`/`burst_assay`):
`cytotoxic_label = 1` if the compound's `AUC < 0.8` (a "response") in at least **10%** of the
cancer cell lines it was tested against; `= 0` otherwise. AUC was used as the threshold metric
(not IC50) because it's normalized to each drug's own screening concentration range, making it
comparable across compounds with wildly different potency scales — raw IC50 isn't directly
comparable across drugs for that reason.

These two numbers (`0.8` AUC cutoff, `10%` coverage) are the single biggest editorial decision in
this file. They're deliberately visible as constants (`AUC_ACTIVE_THRESHOLD`,
`MIN_COVERAGE_FRACTION`) at the top of the build script — change them and rerun if your use case
needs a stricter/looser definition of "cytotoxic." A naive "active in ≥1 line" rule was rejected
because nearly every compound kills *some* cell line at a high enough concentration, which would
make the label nearly meaningless (see the discussion this threshold came out of).

| Column | Description |
|---|---|
| `source` | Always `GDSC` |
| `biological_context` | Always `on-target-cancer-efficacy` |
| `cid`, `canonical_smiles`, `isomeric_smiles`, `molecular_formula`, `molecular_weight` | From PubChem, resolved by drug name (with GDSC's listed synonyms as fallback) |
| `cytotoxic_label` | See derivation above |
| `potency_value` | = `median_auc` (median AUC across all tested lines) |
| `potency_metric` | Always `AUC` |
| `potency_unit` | `unitless_0to1_lower_is_more_sensitive` |
| `n_assays_active` / `n_assays_tested` | Cell lines responding (AUC < 0.8) / total cell lines tested — analogous to Tox21/ToxCast's assay counts, but counting **cell lines**, not assay endpoints |
| `frac_lines_responding` | `n_assays_active / n_assays_tested` — the value the 10% coverage threshold is applied to |
| `median_auc`, `median_ic50_um` | Kept as explicit source-native columns alongside the normalized `potency_value` |
| `assay_description` | States the exact threshold/coverage rule used, so it travels with the data |

## Known limitations

- **Name-based compound resolution**: GDSC's compound table has no PubChem/CAS identifier — drugs
  were matched to PubChem by name (with listed synonyms as fallback). Some entries are internal
  Sanger/company codes (e.g. `GSK...`, `AZD...`, numeric IDs) that don't resolve to PubChem at all
  and are dropped; check the build log for the resolved/total count on each run.
- **GDSC1 vs GDSC2 are separate screens** with different drug panels, cell line coverage, and assay
  technology generations — the long file keeps `gdsc_dataset` so you can split them back apart if
  pooling isn't appropriate for your analysis.
- **`z_score` is not cross-drug comparable** — it's standardized within each drug's own line-to-line
  distribution, kept for reference only.
- This is a targeted oncology drug library (~540 compounds, mostly kinase inhibitors and other
  targeted agents plus some chemotherapy), structurally very different from Tox21/ToxCast's broad
  environmental-chemical libraries — expect little compound overlap with those sources, unlike the
  large Tox21↔ToxCast overlap.
