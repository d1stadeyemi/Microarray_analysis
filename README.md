# Differential Gene Expression Analysis from Microarray Data (GSE47960)<br>
## 📌 Project Overview <br>

This repository contains a complete workflow for identifying differentially expressed genes (DEGs) from the GSE47960 microarray dataset and visualizing the results with standard bioinformatics plots.

* The data used for this study were extracted from the work of [Mitchell et al (2013)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0069374)<br>

Data access on this [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47960)<br>

* The differential gene expression (DEG) analysis can be found...

The pipeline is split into two main parts:

(1) DEG analysis – preprocessing, normalization, and statistical testing.<br>
(2) Visualization – PCA, volcano plots, heatmaps, Venn diagrams, DEG counts, and QC plots.

This repo is designed for reproducibility: anyone can clone it and reproduce the results with a single command.

## 📂 Repository Structure<br>
microarray-DEG-analysis/<br>
│<br>
├── README.md               # This file (project documentation)<br>
├── LICENSE                 # License (MIT/GPL-3.0 recommended)<br>
├── requirements.txt        # List of R packages used<br>
├── run_all.R               # Master script to run the whole pipeline<br>
│<br>
├── data/<br>
│   ├── raw/                # Raw GEO data (not tracked in git)<br>
│   ├── processed/          # Processed expression matrices<br>
│   └── metadata/           # Sample annotation files<br>
│<br>
├── scripts/<br>
│   ├── 01_DEG_analysis.R   # Preprocessing + DEG identification<br>
│   └── 02_visualization.R  # Visualization of results<br>
│<br>
├── results/<br>
│   ├── tables/             # DEG tables in CSV/TSV format<br>
│   └── figures/            # All figures (PNG, PDF, TIFF)<br>
│<br>
└── docs/                   # Optional: RMarkdown reports or extra notes<br>

## ⚙️ Installation & Requirements

You need R (≥ 4.0) and the following packages:<br>
- limma<br>
- GEOquery<br>
- ggplot2<br>
- VennDiagram<br>
- cowplot<br>
- tidyverse<br>

Install packages with:

install.packages(c("ggplot2", "VennDiagram", "cowplot", "tidyverse"))<br>
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")<br>
BiocManager::install(c("limma", "GEOquery", "umap"))<br>

A full list is saved in [requirements.txt](https://github.com/d1stadeyemi/Microarray_analysis/blob/master/requirements.txt)

## 🚀 Quick Start

Clone this repo and run the entire pipeline:

git clone https://github.com/d1stadeyemi/Microarray_analysis<br>
cd microarray-DEG-analysis<br>
Rscript run_all.R<br>

This will:

Download/prepare the dataset<br>
Identify DEGs in each condition<br>
Generate plots and save them to results/figures/<br>
Save DEG tables to results/tables/<br>

## 📊 Key Outputs

QC Plots: Density plot, mean–variance trend<br>
DEG Overlap: Venn diagram of shared DEGs<br>
DEG Counts: Stacked bar plot of up/down DEGs per condition<br>
Volcano Plots: For each condition<br>

## 📑 Citation

If you use this workflow, please cite the original dataset:
Mitchell et al (2013), GEO accession: GSE47960.

## 👩‍💻 Author

Muhammad Ajagbe – Bioinformatics Researcher / Data Scientist <br>
[LinkedIn](www.linkedin.com/in/muhammadajagbe)