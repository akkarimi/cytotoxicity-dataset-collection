# Cytotoxicity Dataset Collection

A curated list of publicly available datasets containing cytotoxicity / cell-viability data, with notes on how to pull the cytotoxicity-specific slice out of each source and where to download it.

Two broad categories, per the general pattern below:

- **Toxicology screening databases** (Tox21, ToxCast, PubChem, CompTox) — cytotoxicity is one endpoint mixed in among hundreds of pathway assays, so you need to **filter** by assay name / target family.
- **Cancer pharmacology databases** (NCI-60, GDSC, CTRP, CCLE, gCSI, PRISM, PharmacoDB, CellMinerCDB, CDRP) — the core screen *is* a cytotoxicity/viability assay, so you don't filter for cytotoxicity — you just avoid pulling in the separate omics (expression, mutation, copy-number) files bundled in the same portal.

---

## Toxicology screening databases (filter required)

### Tox21

> **Already built:** [`data/tox21/`](data/tox21/) has this extracted, ready to use — a long-format
> file (every compound × cell line × technology × timepoint) and a classification-ready summary
> (one row per compound with a binary label), pulled from the AIDs below via
> [`scripts/build_tox21_cytotoxicity.py`](scripts/build_tox21_cytotoxicity.py). See
> [data/tox21/README.md](data/tox21/README.md) for the schema and label derivation.

Cell viability / cytotoxicity assays (e.g. HEK293, HepG2 glo/fluor assays) are separate AIDs from the nuclear-receptor and stress-pathway assays, so filter by assay name. Two real-time assay technologies were run specifically for this purpose, at multiple exposure timepoints (0/8/16/24/32/40 h), each cell line getting its own AID:

- **RealTime-Glo MT Cell Viability Assay** (metabolic activity): e.g. [AID 1224880](https://pubchem.ncbi.nlm.nih.gov/bioassay/1224880) (HEK293, 0h), [AID 1224886](https://pubchem.ncbi.nlm.nih.gov/bioassay/1224886) (HEK293, 24h), [AID 1224889](https://pubchem.ncbi.nlm.nih.gov/bioassay/1224889) (HepG2, 0h), [AID 1224867](https://pubchem.ncbi.nlm.nih.gov/bioassay/1224867) (HepG2, 24h) — full set spans AID 1224867–1224889.
- **CellTox Green Cytotoxicity Assay** (membrane integrity / cell death): e.g. [AID 1224884](https://pubchem.ncbi.nlm.nih.gov/bioassay/1224884) (HEK293, 8h), [AID 1224875](https://pubchem.ncbi.nlm.nih.gov/bioassay/1224875) (HEK293, 24h), [AID 1224878](https://pubchem.ncbi.nlm.nih.gov/bioassay/1224878) (HepG2, 24h), [AID 1224888](https://pubchem.ncbi.nlm.nih.gov/bioassay/1224888) (HEK293, 16h).
- These correspond to the screen published in [Huang et al., "Real-time cell toxicity profiling of Tox21 10K compounds"](https://pmc.ncbi.nlm.nih.gov/articles/PMC5439695/) (PLOS ONE, 2017) — search PubChem AssayName for `"cell viability"` or `"cytotoxicity"` combined with `Tox21` to pull the complete AID list programmatically.
- Apoptosis proxy (related but distinct endpoint): [AID 1347034 — Caspase-3/7 induction in HepG2 cells](https://pubchem.ncbi.nlm.nih.gov/bioassay/1347034).
- **Faster route — CompTox Dashboard:** [Bioactivity search](https://comptox.epa.gov/dashboard/) filtered by assay source "Tox21" and intended target family "Cell Cycle/Viability," or search assay names for `VIABILITY` / `CTOX`.
- **Flat files:** [tox21.gov/data-and-tools](https://tox21.gov/data-and-tools/) — download the Tox21 summary files and filter rows where `assay_name` contains `viability` or `cytotoxicity`. (Combined Tox21/ToxCast summary files are also mirrored on [EPA Figshare](https://epa.figshare.com/articles/dataset/ToxCast_and_Tox21_Summary_Files/6062479).)

### ToxCast

> **Already built:** [`data/toxcast/`](data/toxcast/) has the chemical-level burst output extracted
> and resolved to PubChem/SMILES, in the same normalized schema as [`data/tox21/`](data/tox21/) so
> the two can be combined — see [`data/combined/`](data/combined/) for the pre-built union. Built via
> [`scripts/build_toxcast_cytotoxicity.py`](scripts/build_toxcast_cytotoxicity.py).

The "cytotoxicity burst" is a pre-defined, ready-to-use subset, available at two levels of granularity:

- **Assay-level:** download `Assay_Summary_*.csv` from EPA's [Exploring ToxCast Data](https://www.epa.gov/comptox-tools/exploring-toxcast-data) page and filter rows where the **`burst_assay`** column `== 1` — isolates the ~86 assays EPA classifies as measuring cytotoxicity/cell-stress rather than a specific target.
- **Chemical-level (simplest):** the pre-built **"Download Cytotoxicity Burst Output"** file on the same [Exploring ToxCast Data](https://www.epa.gov/comptox-tools/exploring-toxcast-data) page — one row per chemical, with **`cytotox_median_um`** (median AC50 across burst endpoints) and **`cytotox_lower_bound_um`** (conservative lower-bound estimate) columns already computed. No filtering needed.
- Bulk database (invitrodb, currently v4.2): released via [EPA's high-throughput screening data (gaftp)](https://gaftp.epa.gov/COMPTOX/High_Throughput_Screening_Data/) and mirrored on [EPA Figshare](https://epa.figshare.com/); query the `cytotox` table directly with the [`tcpl`](https://cran.r-project.org/package=tcpl) R package if working from invitrodb/MySQL.

### PubChem BioAssay (general)

- **Search:** use the field-restricted query `"cytotoxicity"[AssayName]` (swap in `"viability"[AssayName]` / `"CC50"[AssayName]` as needed) in [PubChem BioAssay search](https://pubchem.ncbi.nlm.nih.gov/#tab=bioassay) — the `[AssayName]` tag restricts the match to the assay title field rather than the full description, cutting out most false positives.
- **Bulk download:** [PubChem FTP — Bioassay](https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/) has three relevant subfolders per AID batch: `Description/` (assay text — grep this for `cytotoxicity`/`viability`/`CC50`), `Data/` (the actual concentration-response results, CSV/ASN.1/XML), and `Extras/` (deposited metadata). Grep `Description/` first to build your AID list, then pull matching files from `Data/`.

### EPA CompTox Chemicals Dashboard

- Single chemical: on a chemical's page (e.g. [Bisphenol A, DTXSID7020182](https://comptox.epa.gov/dashboard/chemical/details/DTXSID7020182)), open the **Bioactivity** tab → ToxCast/Tox21 summary plot. The **cytotoxicity limit** is drawn as a dashed red vertical line on the assay-hit scatter plot, with the underlying median/lower-bound μM value in the tooltip/legend next to it — this is the same `cytotox_median_um` value as the ToxCast chemical-level file above, just surfaced per-chemical in the UI.
- Bulk/batch chemicals: use [Batch Search](https://comptox.epa.gov/dashboard/batch-search), choose the **"ToxCast Assays: AED"** (or Bioactivity Summary) enhanced data sheet on export, then keep only rows tagged with intended target family "cytotoxicity" or the cytotox burst columns — same fields as the ToxCast bulk file, joined to the Dashboard's chemical identifiers (DTXSID/CASRN).

### NCATS Cyto-Safe dataset

Already cytotoxicity-only by design (3T3 / HEK293 CellTiter-Glo assay) — no extraction needed.

- [AID 1345082 — Cytotoxic profiling in 3T3](https://pubchem.ncbi.nlm.nih.gov/bioassay/1345082)
- [AID 1345083 — Cytotoxic profiling in HEK293](https://pubchem.ncbi.nlm.nih.gov/bioassay/1345083)
- Background: [Cyto-Safe: A Machine Learning Tool for Early Identification of Cytotoxic Compounds in Drug Discovery (J. Chem. Inf. Model., 2025)](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01811)

---

## Cancer pharmacology databases (no filtering — pull the drug-response table)

### NCI-60

Reports GI50 / TGI / LC50 (growth inhibition / total growth inhibition / lethal concentration) per compound–cell line pair — these *are* the cytotoxicity/growth-inhibition endpoints.

- [DTP NCI-60 data download portal](https://wiki.nci.nih.gov/spaces/NCIDTPdata/pages/456425808/NCI-60+Data+Download+-+Previous+Releases) (dtp.cancer.gov data, mirrored on the NCI wiki) — "NCI-60 Growth Inhibition Data" bulk download.
- Alternative interface: [CellMiner](https://discover.nci.nih.gov/cellminer/)

### GDSC (Genomics of Drug Sensitivity in Cancer)

The drug response summary (IC50/AUC per drug–cell line pair) is the cytotoxicity/growth-inhibition readout itself.

- [cancerrxgene.org/downloads](https://www.cancerrxgene.org/downloads/bulk-download) — grab the "Fitted dose response" or "Drug sensitivity IC50" spreadsheet.
- Curve-fitting pipeline (if working from raw plate data): [`gdscIC50` R package](https://github.com/CancerRxGene/gdscIC50)

### CTRP (Cancer Therapeutics Response Portal)

- [CTD² Data Portal, Broad CTRPv2 submission](https://ctd2-dashboard.nci.nih.gov/data/submissions/20151216-broad_ctrpv2/Broad_CTRPv2.html) — raw file: `ftp://caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset.zip`
- R access: [`PharmacoGx`](https://bioconductor.org/packages/release/bioc/html/PharmacoGx.html) → `downloadPSet("CTRPv2")`. Response metric of interest: AUC/IC50 from the viability assay.

### CCLE / gCSI / PRISM

Viability/IC50 screens by nature — pull the drug-response table, not the mutation/expression/copy-number files bundled in the same portal.

- [DepMap portal downloads](https://depmap.org/portal/download/) — "Drug sensitivity" data files (CCLE, and the PRISM Repurposing / Primary Screen datasets).
- [DepMap Repurposing Hub](https://depmap.org/repurposing/) — PRISM drug repurposing screens specifically.
- R access for CCLE/gCSI: [`PharmacoGx::downloadPSet()`](https://bioconductor.org/packages/release/bioc/html/PharmacoGx.html) (e.g. `downloadPSet("CCLE")`, `downloadPSet("gCSI")`).

### PharmacoDB

Aggregates dose-response curves (IC50, EC50, AAC/AUC, Einf, DSS) across CCLE, GDSC, CTRPv2, gCSI, and more, without the omics layers.

- Web interface / batch query: [pharmacodb.pmgenomics.ca](https://pharmacodb.pmgenomics.ca/)
- REST API docs: [pharmacodb.pmgenomics.ca/docs](https://pharmacodb.pmgenomics.ca/docs)
- R access: [`PharmacoGx`](https://bioconductor.org/packages/release/bioc/html/PharmacoGx.html)

### CellMinerCDB

Under "Drug activity" data type (as opposed to expression/proteomics/mutation), export the activity z-scores — cytotoxicity potency values across the integrated cell-line panels.

- [CellMinerCDB (NCATS version)](https://discover.nci.nih.gov/cellminercdb_ncats/)
- [CellMiner — Drug activity z-scores (NCI-60)](https://discover.nci.nih.gov/cellminer/html/drug_zscore.html) / [HTS384 drug activity z-scores](https://discover.nci.nih.gov/cellminer/html/drug_zscore_hts384.html)
- Bulk download: [CellMiner Download Data Sets page](https://discover.nci.nih.gov/cellminer/loadDownload.do) → "Compounds: DTP NCI60"

### NCI Computational Resources — CDRP dataset

Bundles drug response + omics from NCI-60, NCI-ALMANAC, NCI-SCLC, CCLE, GDSC, gCSI, and CTRP together; download the drug-response-only files, not the gene expression files.

- [Cancer Drug Response Prediction Dataset (CDRP)](https://computational.cancer.gov/dataset/cancer-drug-response-prediction-dataset) — check the file manifest on [computational.cancer.gov](https://computational.cancer.gov/view-dataset-finder) for filenames containing `response` or `dose_response` rather than `expression`.

---

## Specialty literature datasets

Datasets from individual papers (ionic liquids, iridium complexes, PHB/PHBV, silver nanoparticles, transcriptome-based cytotoxicity models, etc.) are already cytotoxicity-only by design — the IC50/EC50/CC50 or viability values are the primary output of the study. No filtering step is needed; download the supplementary file or linked repository from each paper's Data Availability statement.

---

## Notes

- Links above point to the primary/official access point for each resource as of 2026-08-19; portal UIs and file layouts change periodically, so if a link moves, search the resource name plus "download" from its home page.
- For R-based workflows, [`PharmacoGx`](https://bioconductor.org/packages/release/bioc/html/PharmacoGx.html) is the single most convenient entry point across GDSC, CTRP, CCLE, gCSI, and PharmacoDB.
