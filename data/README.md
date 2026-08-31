# Built Datasets

Pre-extracted, cytotoxicity-only classification datasets, one directory per source, all sharing a
normalized column set so they can be combined. See each source's own README for how its label and
potency columns were derived — those choices differ meaningfully between sources and are not
interchangeable defaults.

| Source | Dir | Biological context | Compounds | Active | Notes |
|---|---|---|---|---|---|
| Tox21 | [`tox21/`](tox21/) | off-target-safety | 7,671 | 1,512 (20%) | Dedicated real-time HEK293/HepG2 viability panel |
| ToxCast | [`toxcast/`](toxcast/) | off-target-safety | 7,805 | 4,184 (54%) | Cross-target cytotoxicity "burst" (any of ~86 assays) |
| GDSC | [`gdsc/`](gdsc/) | on-target-cancer-efficacy | (see gdsc/README.md) | — | Cancer cell line drug response (AUC-based threshold) |
| Combined | [`combined/`](combined/) | mixed | union | — | Concatenation of the above on shared columns |

## Shared normalized schema

Every `*_classification.csv` file carries at least these columns, which is what makes them
concatenable (see [`scripts/build_combined_cytotoxicity.py`](../scripts/build_combined_cytotoxicity.py)):

| Column | Meaning |
|---|---|
| `source` | Which dataset the row came from |
| `biological_context` | `off-target-safety` (toxicology screening: harm to normal cells) or `on-target-cancer-efficacy` (cancer pharmacology: drug kills cancer cells) — **treat this as a required stratum, not a droppable label**, see below |
| `cid` | PubChem Compound ID |
| `canonical_smiles`, `isomeric_smiles`, `molecular_formula`, `molecular_weight` | From PubChem |
| `cytotoxic_label` | Binary target; derivation differs by source — read that source's README before trusting it as ground truth |
| `potency_value`, `potency_metric`, `potency_unit` | The potency estimate behind the label, with its metric and unit made explicit (AC50/µM, burst-median-AC50/µM, AUC/unitless, ...) — **never compare `potency_value` across sources without checking `potency_metric` first**, they are not the same quantity |
| `n_assays_active` / `n_assays_tested` | Count of positive/tested "contexts" behind the label — an assay endpoint for Tox21/ToxCast, a cancer cell line for GDSC; not directly comparable in magnitude across sources |
| `assay_description` | Free text identifying the panel and (for threshold-based sources) the exact cutoff used |

Each source's own file also keeps its native, source-specific columns (e.g. Tox21's per-40h-assay
counts, ToxCast's `dtxsid`/`casrn`, GDSC's `median_auc`/`gdsc_dataset`) — the normalized columns
are additions for combining, not replacements.

## Why `biological_context` matters more than it looks

Tox21/ToxCast measure whether a chemical **harms normal cells** (a bad outcome — a safety
concern). GDSC (and any future NCI-60/CTRP/CCLE/gCSI/PRISM addition) measures whether a drug
**kills cancer cells** (a good outcome — a measure of efficacy). Both are legitimately
"cytotoxicity" in the raw biological sense (cell death/viability loss), but pooling
`cytotoxic_label=1` across both without keeping `biological_context` would silently merge two
opposite-valence scientific questions. Always stratify splits and models by this column, or filter
to one context, rather than training on the combined file as if it were one homogeneous target.

## Compound overlap across sources

Sources drawing from broad environmental/consumer-chemical libraries (Tox21, ToxCast) overlap with
each other substantially (6,286 of ~7,700–7,800 compounds shared). Sources drawing from targeted
oncology drug libraries (GDSC, and future CTRP/CCLE/gCSI/PRISM/NCI-60) are expected to overlap much
less with the toxicology sources, and heavily with *each other* (many are re-screens of largely the
same oncology compound sets) — check `scripts/build_combined_cytotoxicity.py`'s printed overlap
counts before assuming independence between any two sources here.

## Reproducing

Each source has a `scripts/build_<source>_cytotoxicity.py` that re-downloads from the live
upstream API/file server and rebuilds its `data/<source>/` outputs. Run
`scripts/build_combined_cytotoxicity.py` after any of them to refresh the union.
