# IL10: Host‑response signatures and outcomes (discovery → validation)

> Multi‑omics analysis pipeline for the role of IL10 in the immunobiology of paediatric infection: discovery RNA‑seq → validation cohort → pathway features → survival models → proteome NMDS (DMSR) → IL10 plasma correlations.

---

## Summary
Run the scripts in order `01 → 07`. They produce (i) DEG tables for discovery & validation, (ii) enriched immune pathways & heatmaps + a small predictive RF model, (iii) cross‑cohort correlation & Venn overlaps, (iv) Kaplan–Meier/Cox survival for IL10 (and ETS1), (v) IL10 plasma–GEX correlations with pathway summaries, and (vi) a DMSR metric from proteomics NMDS with figures and exports.

---

## Repository structure

```
.
├── data/
│   ├── globin.depleted.htseq.count.data.csv      # discovery counts
│   ├── phenodata.csv                              # discovery metadata
│   ├── Read1_host transcriptome.xlsx              # validation expression matrix
│   ├── czh_chain_data.csv                         # validation metadata
│   ├── chain.somascan.csv                         # SomaScan proteomics (validation)
│   ├── chain.tbl.csv                              # SomaScan aptamer ↔ gene mapping
│   ├── data.m.csv                                 # long-format validation GEX (gene-level)
│   ├── dfl.disco.csv                              # written by 01 (discovery DEGs)
│   ├── dfl.chain.csv                              # written by 02 (validation DEGs)
│   ├── datax.disco.csv                            # written by 01 (logCPM + meta; wide)
│   └── output/                                    # written by 07 (DMSR figures/objects)
├── 01_discovery_dge_pipeline.R
├── 02_validation_dge_pipeline.R
├── 03_pathway_enrichment_and_feature_selection.R
├── 04_cross_cohort_comparison.R
├── 05_survival_analysis.R
├── 06_correlation_plasma_IL10_gene_expression.R
└── 07_DMSR_derivation_and_correlation_with_GEX_and_survival_time.R
```

---

## Prerequisites

- **R** ≥ 4.2 with the following packages installed:
  - Core: `tidyverse`, `data.table`, `readr`, `readxl`, `stringr`, `tibble`, `reshape2`, `ggplot2`, `ggrepel`, `viridis`, `cowplot`, `gridExtra`
  - Differential expression & annotation: `edgeR`, `limma`, `biomaRt`
  - Enrichment/annotation: `clusterProfiler`, `org.Hs.eg.db`, `ReactomePA`, `enrichplot`, `AnnotationDbi`, `reactome.db`, `GO.db`
  - Modeling & stats: `randomForest`, `pROC`, `survival`, `survminer`
  - Visuals/extras: `hrbrthemes` (optional), `ggdendro`, `ggpubr`, `VennDiagram`, `ggraph`, `igraph`, `ggbeeswarm`, `vegan`, `wesanderson`, `ggh4x`, `scales`


---

## Quick start

1. Place raw input files under `data/` as shown above.
2. Open R (or RStudio) at the repo root.
3. Run the pipeline, step‑by‑step (each script can also be sourced):
   ```r
   source("01_discovery_dge_pipeline.R")      # writes data/dfl.disco.csv and data/datax.disco.csv
   source("02_validation_dge_pipeline.R")     # writes data/dfl.chain.csv
   source("03_pathway_enrichment_and_feature_selection.R")
   source("04_cross_cohort_comparison.R")
   source("05_survival_analysis.R")
   source("06_correlation_plasma_IL10_gene_expression.R")
   source("07_DMSR_derivation_and_correlation_with_GEX_and_survival_time.R")
   ```
4. Find figures and exports under `data/` and `data/output/` (see below).

---

## Pipeline at a glance

| Step | Script | What it produces (high‑level) |
|---|---|---|
| 01 | `01_discovery_dge_pipeline.R` | Discovery cohort DEGs (`dfl.disco.csv`) and a log‑CPM + metadata matrix (`datax.disco.csv`) for downstream heatmaps/RF. |
| 02 | `02_validation_dge_pipeline.R` | Validation cohort DEGs (`dfl.chain.csv`). |
| 03 | `03_pathway_enrichment_and_feature_selection.R` | GO/Reactome enrichment on discovery DEGs, a heatmap of the top immune pathway, and a small Random‑Forest classifier (AUROC + variable‑importance). |
| 04 | `04_cross_cohort_comparison.R` | Cross‑cohort fold‑change correlation (discovery vs validation) and Venn diagrams (up/down overlaps). |
| 05 | `05_survival_analysis.R` | Kaplan–Meier curves (IL10 high/low) for full cohort and non‑malnourished subgroup; Cox PH models for IL10/ETS1 (age‑adjusted). |
| 06 | `06_correlation_plasma_IL10_gene_expression.R` | Correlation between plasma **IL10** (SomaScan) and gene expression; pathway‑level summaries and plots. |
| 07 | `07_DMSR_derivation_and_correlation_with_GEX_and_survival_time.R` | NMDS of selected immune‑related proteins; **DMSR** (distance to survivor centroid) figures; optional DMSR~GEX correlations; saved R objects. |

---

## Detailed script guide

### 01 — Discovery DGE pipeline (`01_discovery_dge_pipeline.R`)
**Inputs**: `data/globin.depleted.htseq.count.data.csv`, `data/phenodata.csv`  
**Core steps**: sample QC/outlier removal; create `edgeR::DGEList`; TMM normalization; `limma-voom`; contrast *survivor vs non‑survivor*; BH FDR; export DEGs.  
**Outputs**:
- `data/dfl.disco.csv` — discovery DEG table (log2FC, p, Ensembl ID, gene symbol).  
- `data/datax.disco.csv` — wide log‑CPM matrix with `serial`/`outcome` columns for heatmaps/RF.

### 02 — Validation DGE pipeline (`02_validation_dge_pipeline.R`)
**Inputs**: validation expression matrix + metadata (see `data/Read1_host transcriptome.xlsx`, `data/czh_chain_data.csv`).  
**Core steps**: align genes, TMM normalization, `limma-voom`, contrast *alive vs dead*; export DEGs.  
**Outputs**: `data/dfl.chain.csv` — validation DEG table.

### 03 — Pathway enrichment, heatmap & feature selection (`03_pathway_enrichment_and_feature_selection.R`)
**Inputs**: `data/dfl.disco.csv` (DEGs) and curated innate immune gene list; `data/datax.disco.csv` + `phenodata.csv` for plotting.  
**Core steps**: GO BP enrichment (clusterProfiler); select top immune pathway (e.g., **negative regulation of immune response**); build patient × gene heatmap; train a small **Random Forest** on pathway genes, report AUROC, threshold, and importance.  
**Outputs**: enrichment plot(s); heatmap; ROC and variable‑importance plots.

### 04 — Cross‑cohort comparison (`04_cross_cohort_comparison.R`)
**Inputs**: `data/dfl.disco.csv`, `data/dfl.chain.csv`.  
**Core steps**: harmonize by Ensembl ID; merge; adjust p‑values; correlate discovery vs validation fold‑changes (Spearman); compute up/down overlaps; draw Venns.  
**Outputs**: correlation scatter; Venn diagrams (up/down).

### 05 — Survival analysis (`05_survival_analysis.R`)
**Inputs**: `data/Read1_host transcriptome.xlsx`, `data/data.m.csv` (long‑format validation GEX with dates/age/MUAC).  
**Core steps**: define IL10 high/low by median log‑expression; Kaplan–Meier curves for full cohort and **MUAC ≥ 11.5 cm** subgroup; Cox PH for **IL10** and **ETS1** (age‑adjusted).  
**Outputs**: K–M plots; Cox summaries (HRs, CIs).

### 06 — Plasma IL10 ↔ gene expression (`06_correlation_plasma_IL10_gene_expression.R`)
**Inputs**: `data/chain.somascan.csv`, `data/chain.tbl.csv`, `data/czh_chain_data.csv`, `data/data.m.csv`.  
**Core steps**: map aptamers to genes, merge plasma **IL10 RFU** with GEX; compute gene‑wise correlations; summarize by Reactome pathway (immune‑focused); plot pathway‑level slopes/CI and R².  
**Outputs**: correlation tables/plots by pathway.

### 07 — DMSR from proteome NMDS (`07_DMSR_derivation_and_correlation_with_GEX_and_survival_time.R`)
**Inputs**: `data/chain.somascan.csv`, `data/chain.tbl.csv`, `data/czh_chain_data.csv` (+ optional `data/data.m.csv` for DMSR~GEX).  
**Core steps**: curate immune‑related aptamers; per‑record means; **NMDS** (k=2) on scaled proteome; group ellipses & centroids; define **DMSR** = Euclidean distance to the **alive** centroid; compare DMSR across outcomes; (optional) Spearman correlations with GEX and GO annotations.  
**Outputs (in `data/output/`)**:
- `nmds_ordination.png` — NMDS with ellipses/centroids.  
- `dmsr_by_outcome.png` — DMSR vs outcome with CIs.  
- `dmsr_gex_spearman.csv` & `dmsr_gex_spearman.png` (if GEX provided).  
- `dmsr_nmds_objects.rds` — saved `nmds`, `points`, `centroids`, `points_dist` objects.

---

## Reproducing key figures

- **Discovery vs validation correlation** — run step 04 to generate the scatter with Spearman *r* and dashed LM fit.  
- **Immune pathway heatmap** — run step 03; heatmap shows scaled expression of top GO term across patients.  
- **K–M curves (IL10)** — run step 05; shows survival difference for IL10 high vs low overall and in the non‑malnourished subgroup.  
- **DMSR figures** — run step 07; `data/output/` contains NMDS and DMSR plots plus serialized objects.

---

## Tips 
- Ensure sample identifiers in `phenodata.csv` and counts align exactly; the discovery script will stop if they do not.  
- The discovery script drops a small, hard‑coded outlier list; adjust if your dataset differs.  
- For enrichment, make sure internet access for `biomaRt`/annotation is available if you recompute annotations.  
- The Random‑Forest in step 03 is for **exploration** (ntree=10,000 by default); treat its AUROC as internal validation only.

---



---

## License

- **Code**: MIT License (see `LICENSE`).  
- **Figures/Docs**: CC BY 4.0.  
- **Data**: follow the terms of the original data providers; do **not** re‑distribute restricted data in this repository.

---

## Acknowledgements

We thank the children and their caregivers for their participation and trust. We are indebted to the clinicians, nurses, and ward staff who supported enrolment, clinical characterisation, and follow-up, and to the laboratory teams who handled sample logistics, processing, sequencing/proteomics, and data curation. We acknowledge the data management and field teams for rigorous quality assurance throughout. We offer particular thanks to the CHAIN Network investigators and site teams for outstanding collaborative spirit - freely sharing ideas, protocols, and samples, and offering constructive critique at every stage. This work reflects the collective effort of many colleagues across sites and disciplines.
