# scRNA-seq Wound Healing Reproducibility Study

Reproducibility study of Keiser et al. (2025): *Using R, Seurat, and CellChat to Analyze a Single-Cell Transcriptomics Dataset of Mouse Skin Wound Healing*. Published in the Journal of Visualized Experiments, DOI: [10.3791/67266](https://doi.org/10.3791/67266).

This repository contains the R code used to reproduce the single-cell RNA-seq analysis pipeline described in the paper, applied to a publicly available spatio-temporal mouse excisional skin wound healing dataset (Hu et al., 2023, GEO: GSE204777).

---

## Overview

The original paper presents a step-by-step bioinformatics workflow for analyzing single-cell transcriptomics data using **Seurat** and **CellChat** in RStudio. This study reproduces seven core methods from that protocol:

1. Environment setup and package installation
2. Data loading and quality control (QC)
3. Seurat clustering and cell type annotation
4. Fibroblast subtype analysis
5. Module scoring (wound healing phase-specific gene sets)
6. Cell-cell interaction analysis using CellChat (D1 vs D14)
7. Multi-batch dataset integration using RPCA

---

## Repository Structure

```
scRNA-wound-healing/
├── Final_Project.Rmd          # Main R Markdown file with full pipeline
├── dataset_cluster_markers.csv
├── fibroblast_cluster_markers.csv
├── doublet_scores.png
├── elbow_merged.png
├── umap_unintegrated.png
├── umap_integrated.png
├── umap_merged_clusters.png
├── umap_merged_annotated.png
├── cellchat_scatter_D1.png
├── cellchat_scatter_D14.png
├── collagen_circle_D1.png
├── collagen_circle_D14.png
├── cellchat_compare.png
└── README.md
```

---

## Data

The dataset used is **GSE204777** from NCBI GEO, a spatio-temporal single-cell RNA-seq dataset of mouse excisional skin wounds originally published by Hu et al. (2023).

- **Batch 1 (b1):** GSM6190913 — used for the primary single-batch analysis
- **Batch 3 (b3):** GSM6190915 — used for the multi-batch integration analysis

Each batch consists of three files:
- `*_barcodes.tsv.gz`
- `*_features.tsv.gz`
- `*_matrix.mtx.gz`

Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204777

Place batch files in the working directory before running the pipeline.

---

## Requirements

### R Version
R 4.4.1

### R Packages

**CRAN:**
```r
install.packages(c("devtools", "readxl", "openxlsx", "tidyverse", "scCustomize", "Seurat"))
```

**Bioconductor:**
```r
BiocManager::install(c("NMF", "ComplexHeatmap", "BiocNeighbors",
                        "SingleCellExperiment", "circlize",
                        "edgeR", "scDblFinder"))
```

**GitHub:**
```r
devtools::install_github("jinworks/CellChat")
```

---

## How to Run

1. Clone this repository:
```bash
git clone https://github.com/mayankgandhi13/scRNA-wound-healing.git
cd scRNA-wound-healing
```

2. Download the GEO data files for GSM6190913 (b1) and GSM6190915 (b3) and place them in the working directory.

3. Open `Final_Project.Rmd` in RStudio and set the working directory to the project folder:
```r
setwd("/path/to/scRNA-wound-healing")
```

4. Run the chunks sequentially, or knit the full document. Intermediate RDS checkpoints are saved at the end of each major step so you can resume without re-running the entire pipeline.

### Checkpoints

| File | Description |
|------|-------------|
| `dataset_post_Method2.rds` | Post-QC single batch (b1) |
| `dataset_post_Method3.rds` | Post-clustering and annotation |
| `dataset_fibroblast_post_Method4.rds` | Fibroblast subset |
| `dataset_post_Method5.rds` | Post-module scoring |
| `merged_dataset_post_Method7.rds` | Post-integration (b1 + b3) |
| `cellchat_D1.rds` | CellChat object for day 1 |
| `cellchat_D14.rds` | CellChat object for day 14 |
| `cellchat_merged.rds` | Merged CellChat comparison object |

---

## Pipeline Summary

### Quality Control
- Mitochondrial gene filtering: cells with > 25% mt genes removed
- Low-feature filtering: cells with < 200 genes removed
- Doublet detection and removal using `scDblFinder` (threshold: 0.25)
- HTO demultiplexing using `HTODemux` (positive.quantile = 0.99)

### Clustering and Annotation
- Normalization, variable feature selection, scaling, and PCA
- Clustering: 13 PCA dimensions, resolution = 0.1
- UMAP: seed = 123 for reproducibility
- Cell types identified: Macrophage, Neutrophil, Fibroblast, Epithelial cell, Endothelial cell, T cell, Smooth muscle cell

### Fibroblast Subtype Analysis
- Subset from main object and re-clustered using 9 PCA dimensions
- 3 fibroblast subtypes identified with distinct temporal expression patterns

### Module Scoring
- Three gene sets applied: Inflammatory, Proliferative, Resolution
- Phase-specific genes from Wietecha et al. (2023)
- Scored across both DPW groups and major cell types

### CellChat (D1 vs D14)
- Compared cell-cell communication at day 1 (inflammatory) vs day 14 (resolution)
- Focused on the COLLAGEN signaling pathway
- Fibroblasts show dramatically increased outgoing interactions at D14

### Dataset Integration
- Batches b1 and b3 merged and integrated using RPCA via `IntegrateLayers`
- Post-integration clustering resolved 10 clusters vs 8 pre-integration

---

## Key Results

- Immune cells (neutrophils and macrophages) dominate early wound time points (D1, D3)
- Fibroblasts expand significantly during resolution (D14)
- Three fibroblast subtypes show distinct temporal profiles across the healing continuum
- Collagen signaling from fibroblasts broadens from immune-restricted targets at D1 to multiple cell types at D14
- RPCA integration effectively removed batch effects between b1 and b3

---

## References

Keiser, S., Botello, N., Cruz, E., and Wietecha, M. S. (2025). Using R, Seurat, and CellChat to analyze a single-cell transcriptomics dataset of mouse skin wound healing. *Journal of Visualized Experiments, 222*, e67266. https://doi.org/10.3791/67266

Hu, K. H. et al. (2023). Transcriptional space-time mapping identifies concerted immune and stromal cell patterns and gene programs in wound healing and cancer. *Cell Stem Cell, 30*(6), 885-903. https://doi.org/10.1016/j.stem.2023.04.004

Hao, Y. et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology, 42*(2), 293-304.

Jin, S. et al. (2021). Inference and analysis of cell-cell communication using CellChat. *Nature Communications, 12*(1), 1088.

Jin, S., Plikus, M. V., and Nie, Q. (2024). CellChat for systematic analysis of cell-cell communication from single-cell transcriptomics. *Nature Protocols, 20*(1), 180-219.

Germain, P. L. et al. (2021). Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research, 10*, 979.

Wietecha, M. S. et al. (2023). Phase-specific signatures of wound fibroblasts and matrix patterns define cancer-associated fibroblast subtypes. *Matrix Biology, 119*, 19-56.

---

## Acknowledgments

This project was completed as part of a graduate-level bioinformatics reproducibility study at Northeastern University. Original data and protocol credit belongs to Keiser et al. (2025) and Hu et al. (2023).
