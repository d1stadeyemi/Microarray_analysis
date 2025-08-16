if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("umap")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("VennDiagram")
install.packages("patchwork")
install.packages("cowplot")
install.packages("dplyr")
install.packages("tidyr")
install.packages("grid")
install.packages("gridExtra")

# Load required packages for Differential Gene Expression Analysis
library(GEOquery)
library(limma)
library(umap)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(patchwork)
library(cowplot)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)

# Download expression dataset (in matrix format), metadata and the annotation if available
gset <- getGEO("GSE47960", GSEMatrix =TRUE, AnnotGPL=TRUE)

# Since there are only 1 platform in the GEO Dataset,
# then the code below is not needed
# As it selects, if there are more than 1 platform, the dataset of the first platform
# if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1

# Extract the ExpressionSet of the only available platform in the GEO dataset
gset <- gset[[1]]

# Convert the feature variable labels (column names) into a standard R compatible variable names
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Create groups for samples
# Group 0 is for the control (mock)
# while group 1 to 4 for treatment samples

# Define a function that uses number patterns to classify the samples into groups
classify_samples <- function(sample_names, pattern_map) {

  # Start with the default group for samples
  sml <- rep("X", length(sample_names))

  # Create groups based on the pattern map
  for (pat in names(pattern_map)) {
    sml[grepl(pat, sample_names, ignore.case=TRUE)] <- pattern_map[[pat]]
  }

  return(sml)
}

# Assign the samples' column to a variable
sample_names <- gset$`infection code:ch1`

# Define patterns and group label
pattern_map <- list(
  "mock" = "0",
  "H1N1" = "1",
  "dORF6" = "2",
  "BatSRBD" = "3",
  "icSARS" = "4"
)

# Apply function to classify samples into groups
sml <- classify_samples(sample_names, pattern_map)

# Prepare expression dataset for log2 transformation if needed
# by checking if the expression data looks like raw data or not
# If yes, then conduct a log2 transformation
# Extract the expression matrix from the gset object
ex <- exprs(gset)

# Calculate different percentiles of the expression values and omit NA
# Then convert the result into numeric values
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

# Check the percentiles of the expression data to determine if raw or not
LogC <- (qx[5] > 100) || # Check if the 99th percentile is very high
          (qx[6]-qx[1] > 50 && qx[2] > 0) # Check if the range is very high and most values are positive

# Applied log2 transformation if the above conditions are true
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# Prepare to assign samples to groups
# Convert vector (suppose group labels) to factor which is used for categorical variables
gs <- factor(sml)

# Create an R friendly variable names to facilitate downstream usage
groups <- make.names(c('mock', 'H1N1', 'dORF6', 'BatSRBD', 'icSARS'))

# Convert the group labels to R valid variable names
levels(gs) <- groups

# Add a new "group" column to define sample groups
gset$group <- gs

# Create a no-intercept design matrix with one column per group
design <- model.matrix(~ group + 0, gset)

# Rename design matrix columns for consistency
colnames(design) <- levels(gs)

# Data Cleaning: Keep only rows (probes/genes) without missing (NA) values in all samples
gset <- gset[complete.cases(exprs(gset)), ]

# Fit a linear model
fit <- lmFit(gset, design)

# Create contrast: Find the difference in expression between mock and treatment groups
cont.matrix <- makeContrasts(
  H1N1_vs_mock = H1N1 - mock,
  dORF6_vs_mock = dORF6 - mock,
  BatSRBD_vs_mock = BatSRBD - mock,
  icSARS_vs_mock = icSARS - mock,
  levels = design
)

# Apply the difference (contrast) to the fitted model
# to isolate only the difference between the two groups (control and treatment) for each gene
fit2 <- contrasts.fit(fit, cont.matrix)

# Apply empirical bayes moderation to the linear model results in fit2
# To reduce linear model estimation noise
fit2 <- eBayes(fit2)

# Compute the main output table for each sample
tT_H1N1 <- topTable(
          fit2,                   # The moderated model
          coef = "H1N1_vs_mock",  # define which sample to select from
          adjust = "fdr",         # Adjust p-values for multiple testing using False Discovery Rate
          sort.by = "B",          # Sort genes by B-statistic (probability of being Differentially Expressed)
          number = Inf            # Get all genes
          )

tT_dORF6 <- topTable(fit2, coef="dORF6_vs_mock" ,adjust="fdr", sort.by="B", number=Inf)
tT_BatSRBD <- topTable(fit2, coef="BatSRBD_vs_mock", adjust="fdr", sort.by="B", number=Inf)
tT_icSARS <- topTable(fit2, coef="icSARS_vs_mock", adjust="fdr", sort.by="B", number=Inf)

# Extract only the important columns from the main output (final result; tT)
tT_H1N1 <- subset(tT_H1N1, select = c("ID","adj.P.Val","P.Value","t","B","logFC","GenBank.Accession","Gene.symbol","Gene.title"))
tT_dORF6 <- subset(tT_dORF6, select = c("ID","adj.P.Val","P.Value","t","B","logFC","GenBank.Accession","Gene.symbol","Gene.title"))
tT_BatSRBD <- subset(tT_BatSRBD, select = c("ID","adj.P.Val","P.Value","t","B","logFC","GenBank.Accession","Gene.symbol","Gene.title"))
tT_icSARS <- subset(tT_icSARS, select = c("ID","adj.P.Val","P.Value","t","B","logFC","GenBank.Accession","Gene.symbol","Gene.title"))

# Extract the DEGs (adjusted p-value <0.05 and absolute value of logFC > 1)
degs_H1N1 <- subset(tT_H1N1, adj.P.Val < 0.05 & abs(logFC) > 1)
degs_dORF6 <- subset(tT_dORF6, adj.P.Val < 0.05 & abs(logFC) > 1)
degs_BatSRBD <- subset(tT_BatSRBD, adj.P.Val < 0.05 & abs(logFC) > 1)
degs_icSARS <- subset(tT_icSARS, adj.P.Val < 0.05 & abs(logFC) > 1)

# Another way of extracting upregulated and downregulated genes using limma
# summarize test results as "up", "down" or "not expressed"
# dT <- decideTests(fit2, adjust.method = "fdr", p.value = 0.01, lfc = 0)

write.table(degs_H1N1, file = "DEGs_H1N1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(degs_dORF6, file = "DEGs_dORF6.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(degs_BatSRBD, file = "DEGs_BatSRBD.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(degs_icSARS, file = "DEGs_icSARS.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Get common DEGs in paired conditions
common_degs_H1N1_dORF6 <- intersect(degs_H1N1$ID, degs_dORF6$ID)
common_degs_H1N1_BatSRBD <- intersect(degs_H1N1$ID, degs_BatSRBD$ID)
common_degs_H1N1_icSARS <- intersect(degs_H1N1$ID, degs_icSARS$ID)
common_degs_dORF6_BatSRBD <- intersect(degs_dORF6$ID, degs_BatSRBD$ID)
common_degs_dORF6_icSARS <- intersect(degs_dORF6$ID, degs_icSARS$ID)
common_degs_BatSRBD_icSARS <- intersect(degs_BatSRBD$ID, degs_icSARS$ID)

# Get DEGs that are common in all conditions
common_degs <- Reduce(intersect, list(degs_H1N1$ID, degs_dORF6$ID, degs_BatSRBD$ID, degs_icSARS$ID))