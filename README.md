# MOFA-MPR
## An Integrative Multi-Omics Framework Revealing Distinct Molecular Programmes of Breast Cancer Transformation and Progression


### Overview

This repository contains the computational workflow used for multi-omics biomarker discovery and validation in breast cancer. The analytical framework integrates molecular data from The Cancer Genome Atlas Breast Invasive Carcinoma cohort (TCGA-BRCA) using Multi-Omics Factor Analysis (MOFA+), followed by statistical evaluation of latent factors, feature-level convergence, external transcriptomic validation, and single-cell RNA sequencing (scRNA-seq) contextualization.

## Workflow

```text
TCGA-BRCA Multi-Omics Data
        │
        ▼
Data Download and Preprocessing
        │
        ▼
Cross-Omics Sample Alignment
        │
        ▼
MOFA+ Integration
        │
        ▼
Factor Statistical Analysis
        │
        ▼
Panel Selection
        │
        ├──────────────────────────────┐
        ▼                              ▼
External Validation              scRNA-seq Contextualization
        │                              │
        ├── Tumour vs Normal           ├── Harmony Integration
        ├── Early vs Late              ├── Cell-Type Annotation
        └── Stage I vs Stage II        ├── Gene Contextualization
                                       └── Compartment/Stage Analysis
```

## Datasets

The analysis uses publicly available datasets from TCGA and the Gene Expression Omnibus (GEO):

- TCGA-BRCA
- GSE109169
- GSE15852
- GSE20685
- GSE129551
- GSE96058
- GSE161529

## Requirements

The analyses were performed primarily in R.

Key packages include:

- MOFA2
- TCGAbiolinks
- limma
- GEOquery
- Seurat
- harmony
- dplyr
- ggplot2

## Citation

The associated manuscript citation will be added following publication.
