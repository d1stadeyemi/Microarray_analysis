# install_requirements.R
# Script to install all required R and Bioconductor packages

# ---- CRAN packages ----
cran_packages <- c(
    "ggplot2",
    "ggrepel",
    "VennDiagram",
    "patchwork",
    "cowplot",
    "dplyr",
    "tidyr",
    "grid",
    "gridExtra"
)

# ---- Bioconductor packages ----
bioc_packages <- c(
  "limma",
  "GEOquery",
  "umap" 
)

# ---- Install CRAN packages ----
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

invisible(lapply(cran_packages, install_if_missing))

# ---- Install Bioconductor packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
invisible(BiocManager::install(bioc_packages, ask = FALSE, update = TRUE))

cat("\n✅ All required packages are installed!\n")