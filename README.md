# pediatric-KFD-immune-states

Analysis code for **"Single-cell and spatial transcriptomics characterization of
interferon-associated immune states in pediatric Kikuchi–Fujimoto disease."**

## Data

| Modality | Accession |
|---|---|
| PBMC scRNA-seq, scTCR-seq (KFD) | GEO GSE345328 |
| GeoMx DSP, lymph node (KFD) | GEO GSE345380 |
| PBMC scRNA-seq (healthy pediatric controls) | as reported in the source study (ref. 23) |

## Repository layout

| Script | Analysis | Methods section |
|---|---|---|
| `01_scRNAseq/01_preprocessing_QC.R` | SoupX ambient-RNA correction, HTO demultiplexing, per-sample adaptive QC, scDblFinder | *Single-cell RNA-seq data processing, clustering, and annotation* |
| `01_scRNAseq/02_integration_annotation.R` | Harmony integration, Louvain clustering, lineage re-clustering and annotation | *Single-cell RNA-seq data processing, clustering, and annotation* |
| `01_scRNAseq/03_pseudobulk_DEG_pathway.R` | Pseudobulk KFD vs HC differential expression, GO/KEGG enrichment | *Differential expression analysis*; *Pathway enrichment analysis* |
| `01_scRNAseq/04_module_scores_signatures.R` | Hallmark IFN-α module score; severity-associated signature gene sets | *Type I IFN module scoring*; *Signature gene set construction* |
| `01_scRNAseq/05_cellphonedb.py` | CellPhoneDB ligand–receptor inference | *Cell–cell interaction analysis* |
| `02_TCR_repertoire/TCR_scirpy.py` | scirpy chain QC, TCRdist clonotype clusters, clonal expansion | *scTCR-seq analysis* |
| `03_GeoMx/geomx_DEG_AUCell.R` | ROI-level DE, pathway enrichment, AUCell projection of PBMC signatures | *Spatial transcriptomic profiling (GeoMx DSP) and analysis* |

## Software

- **Cell Ranger** v7.1.0 (`cellranger multi`; GRCh38-2020-A, V(D)J GRCh38-alts-ensembl-7.1.0)
- **R** (v4.4.0) — Seurat v5, SoupX, scDblFinder, harmony, clusterProfiler,
  msigdbr, escape (AUCell).
- **Python** — scirpy, muon, scanpy (v3.10); CellPhoneDB v5.0.1 (v3.9).
- **GeoMx DSP Data Analysis Suite** v3.1.0 (QC and Q3 normalization).

Statistical tests and multiple-testing correction follow the *Statistical
analysis* section of the manuscript.
