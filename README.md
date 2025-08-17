# DEG Analysis
The🧬 Differential Gene Expression Analysis from Microarray Data (GSE47960)
📌 Project Overview

This repository contains a complete workflow for identifying differentially expressed genes (DEGs) from the GSE47960 microarray dataset and visualizing the results with standard bioinformatics plots.

The data used for this study were extracted from the work of [Mitchell et al (2013)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0069374)

Data access on this [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47960)

The differential gene expression (DEG) analysis can be found...

The pipeline is split into two main parts:

DEG analysis – preprocessing, normalization, and statistical testing.

Visualization – PCA, volcano plots, heatmaps, Venn diagrams, DEG counts, and QC plots.

This repo is designed for reproducibility: anyone can clone it and reproduce the results with a single command.

📂 Repository Structure
microarray-DEG-analysis/
│
├── README.md               # This file (project documentation)
├── LICENSE                 # License (MIT/GPL-3.0 recommended)
├── requirements.txt        # List of R packages used
├── run_all.R               # Master script to run the whole pipeline
│
├── data/
│   ├── raw/                # Raw GEO data (not tracked in git)
│   ├── processed/          # Processed expression matrices
│   └── metadata/           # Sample annotation files
│
├── scripts/
│   ├── 01_DEG_analysis.R   # Preprocessing + DEG identification
│   └── 02_visualization.R  # Visualization of results
│
├── results/
│   ├── tables/             # DEG tables in CSV/TSV format
│   └── figures/            # All figures (PNG, PDF, TIFF)
│
└── docs/                   # Optional: RMarkdown reports or extra notes

⚙️ Installation & Requirements

You need R (≥ 4.0) and the following packages:

limma

GEOquery

edgeR

ggplot2

pheatmap

VennDiagram

cowplot

tidyverse

Install packages with:

install.packages(c("ggplot2", "pheatmap", "cowplot", "tidyverse"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "GEOquery", "edgeR"))


A full list is saved in requirements.txt.

🚀 Quick Start

Clone this repo and run the entire pipeline:

git clone https://github.com/YOURUSERNAME/microarray-DEG-analysis.git
cd microarray-DEG-analysis
Rscript run_all.R


This will:

Download/prepare the dataset

Identify DEGs in each condition

Generate plots and save them to results/figures/

Save DEG tables to results/tables/

📊 Key Outputs

QC Plots: Density plot, mean–variance trend

DEG Overlap: Venn diagram of shared DEGs

DEG Counts: Stacked bar plot of up/down DEGs per condition

Volcano Plots: For each condition

PCA Plot: Clustering of samples

Heatmap: Regulation patterns of common DEGs

Example (preview):




📑 Citation

If you use this workflow, please cite the original dataset:
Moldoveanu et al. (2013), GEO accession: GSE47960.

👩‍💻 Author

Your Name – Bioinformatics Researcher / Data Scientist

LinkedIn 