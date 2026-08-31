# Combined Cytotoxicity Dataset

`combined_cytotoxicity_classification.csv` — the union of
[`data/tox21/tox21_cytotoxicity_classification.csv`](../tox21/) and
[`data/toxcast/toxcast_cytotoxicity_classification.csv`](../toxcast/) on their shared normalized
columns. Built by [`scripts/build_combined_cytotoxicity.py`](../../scripts/build_combined_cytotoxicity.py).

15,476 rows = 7,671 (Tox21) + 7,805 (ToxCast). One row per **compound × source** — a compound
tested in both panels appears twice (6,286 CIDs overlap), each row carrying that source's own
label and potency estimate. This is deliberate: see [data/toxcast/README.md § Combining](../toxcast/README.md#combining-with-the-tox21-dataset)
for why the two labels aren't merged into one row automatically.

| Column | Meaning |
|---|---|
| `source` | `Tox21` or `ToxCast` |
| `cid` | PubChem Compound ID |
| `canonical_smiles`, `isomeric_smiles`, `molecular_formula`, `molecular_weight` | From PubChem |
| `cytotoxic_label` | Binary target — derivation differs by source, see each source's README |
| `potency_um` | Primary potency estimate — Tox21: lowest 40h AC50 among active assays; ToxCast: median AC50 across burst-flagged hits |
| `n_assays_active` / `n_assays_tested` | Assay-count context behind the label — not directly comparable across sources (Tox21: out of 4 assays at 40h; ToxCast: out of ~86 burst assays) |
| `assay_description` | Free text identifying which panel produced the row |

**Before treating this as one training set:** the two labeling rules aren't equivalent (Tox21:
active in ≥1 of 4 dedicated viability assays at 40h; ToxCast: active in ≥1 of ~86 assays flagged
as cytotoxicity-associated across many different technologies/cell types), and their class
balance differs sharply (Tox21 ~20% active, ToxCast ~54% active). Treat `source` as a required
feature/stratum, not a nuisance column to drop — at minimum, stratify any train/test split by it,
and consider it a multi-task or domain-adaptation problem rather than one homogeneous label.
