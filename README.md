# Cross-etiology transcriptomic conservation in hepatocellular carcinoma

This repository contains the full R workflow used to reproduce the analyses in the manuscript:

**“Cross-etiology transcriptomic conservation in hepatocellular carcinoma reveals opposing proliferation and hepatocyte-loss programs validated across cohorts.”**

## Overview

We integrate public transcriptomic cohorts to identify conserved HCC tumor-state programs across viral etiologies (HBV and HCV), distill these programs into compact gene modules, and validate module-based scores in independent cohorts. The workflow includes:

- Cohort acquisition and standardized metadata curation (GEO)
- Gene-level preprocessing (probe → gene symbol)
- Hallmark pathway scoring (GSVA/ssGSEA) and tumor vs non-tumor contrasts
- HBV injury-axis derivation and adjustment for proliferative state (E2F/G2M)
- Cross-cohort gene-level meta-analysis
- Module construction (Panel A/Panel B) and computation of:
  - **ProlifHubScore**
  - **HepLossScore**
  - **HCCStateScore = ProlifHubScore − HepLossScore**
- External validation in GSE14520 (ROC/AUC)
- Matched-size random gene-set null test
- Automated generation of figures and tables

## Datasets

Public datasets used:
- GEO: **GSE121248** (HBV-associated HCC; tumor vs paired non-tumor/adjacent)
- GEO: **GSE41804** (HCV-associated HCC; tumor vs paired non-tumor/adjacent; IL28B genotype)
- GEO: **GSE83148** (chronic HBV hepatitis with ALT/AST/HBV-DNA strata; non-tumor)
- GEO: **GSE38941** (HBV-associated acute liver failure; non-tumor)

External validation:
- **GSE14520** (independent HCC cohort). For reproducibility and stability, this workflow expects a pre-assembled expression matrix with a `#CLASS` row (tumor vs non-tumor) placed at:
  - `data/external/GSE14520_expression_matrix_reordered.csv`

Web-based cross-platform check (manual logging):
- **GEPIA3** (LIHC; median split 50%; OS/DFS in months)

## Repository structure

- `01_*.R` … `10_*.R`: pipeline scripts in execution order  
- `R/_shared.R`: shared helpers (directory creation, phenotype parsing, mapping/scoring utilities, GSVA API compatibility)
- `data/`: raw, curated, and external input files
- `results/`, `meta/`, `paper_package/`: intermediate and manuscript-ready outputs

## Requirements

- R (recommended ≥ 4.2)
- Bioconductor packages (limma, GSVA, msigdbr, annotation packages)
- PDF tools (optional): `qpdf` R package is used to combine PDFs

## Running the pipeline

From the repo root:

```bash
Rscript 01_data_import.R
Rscript 02_metadata_curation.R
Rscript 03_preprocessing_gene_mapping.R
Rscript 04_hallmark_scoring.R
Rscript 05_injury_axis_models.R
Rscript 06_meta_analysis.R
Rscript 07_module_construction.R
Rscript 08_validation_GSE14520.R
Rscript 09_random_null_test.R
Rscript 10_figures_tables.R

[![DOI](https://zenodo.org/badge/1172222258.svg)](https://doi.org/10.5281/zenodo.18918427)
