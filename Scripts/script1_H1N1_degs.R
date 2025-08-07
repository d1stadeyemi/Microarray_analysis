# If required install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the required packages
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("umap")
install.packages("ggplot2")
install.packages("ggrepel")

# Load required packages for Differential Gene Expression Analysis
library(GEOquery)
library(limma)
library(umap)
library(ggplot2)
library(ggrepel)

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

# split samples into groups
# Control samples ("mock") are assigned to group 0
# while H1N1 treatment samples are assigned to group 1
# and others are assigned to group X
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX11111111111111111",
        "111000000000000000000000000000000000000XXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        "XXXXXXXXXXXXX") # NB: paste0 combines several long stings into a single long string

# Slipt the string "gsms" into individual character and extract them from the list
sml <- strsplit(gsms, split = "")[[1]]

# Find the indices of non-X samples
sel <- which(sml != "X")

# Filter group labels
sml <- sml[sel]

# Filter columns in the expression dataset
gset <- gset[ ,sel]

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
groups <- make.names(c("0","1"))

# Convert the group labels to R valid variable names
levels(gs) <- groups

# Add a new "group" column to define sample groups
gset$group <- gs

# Create a no-intercept design matrix with one column per group
design <- model.matrix(~group + 0, gset)

# Rename design matrix columns for consistency
colnames(design) <- levels(gs)

# Data Cleaning: Keep only rows (probes/genes) without missing (NA) values in all samples
gset <- gset[complete.cases(exprs(gset)), ]

# Fit a linear model
fit <- lmFit(gset, design)

# Compare group X0 and group X1
cts <- paste(groups[1], groups[2], sep = "-")

# Find the difference in expression between the two groups
cont.matrix <- makeContrasts(contrasts=cts, levels=design)

# Apply the difference (contrast) to the fitted model
# to isolate only the difference between the two groups (control and treatment) for each gene
fit2 <- contrasts.fit(fit, cont.matrix)

# Apply empirical bayes moderation to the linear model results in fit2
# To reduce linear model estimation noise
fit2 <- eBayes(fit2, 0.01)

# Compute the main output table
tT <- topTable(
          fit2,               # The moderated model
          adjust = "fdr",     # Adjust p-values for multiple testing using False Discovery Rate
          sort.by = "B",      # Sort genes by B-statistic (probability of being Differentially Expressed)
          number = Inf        # Get all genes
          )

# Extract only the important columns from the main output (final result; tT)
tT <- subset(tT, select = c("ID","adj.P.Val","P.Value","t","B","logFC","GenBank.Accession","Gene.symbol","Gene.title"))

# Extract the DEGs (adjusted p-value <0.05 and absolute value of logFC > 1)
degs <- subset(tT, adj.P.Val < 0.05 & abs(logFC) > 1)

# Extract genes that are upregulated and downregulated
upregulated_genes <- subset(degs, logFC > 1)
downregulated_genes <- subset(degs, logFC < -1)

# Another way of extracting upregulated and downregulated genes using limma
# summarize test results as "up", "down" or "not expressed"
# dT <- decideTests(fit2, adjust.method = "fdr", p.value = 0.01, lfc = 0)

write.table(degs, file = "DEGs_H1N1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Visualizations and Quality Control

# Plot a histogram to visualize the distribution of Adjusted p-values
hist(tT$adj.P.Val, col = "grey", border = "white",
      xlab = "P-adj value", ylab = "Number of genes",
      main = "P-adj value distribution")

# Filters out genes that failed model fitting (no F-statistic)
t.good <- which(!is.na(fit2$F))

# Create a quantile-quantile (Q-Q) plot of the moderated t-statistics calculated by limma
# to assess the distribution and reliability of the test statistics
# In other words, draw a Q-Q plot to assess if the t-statistics
# behaved as expected under the null hypothesis
qqt(fit2$t[t.good], fit2$df.total[t.good], main = "Moderated t statistic")

# Create a volcano plot
# Ensure required columns exist for the plot
# by taking the -log10 of the adjusted pvalues
tT$logP <- -log10(tT$adj.P.Val)

# Plot
plot(
     tT$logFC, tT$logP,
     col="grey", pch=20,
     xlab="log fold change",
     ylab="-log10(Adj-pvalue)",
     main="Volcano Plot of DEGs",
     )

# Add threshold lines
abline(h = -log10(0.05), col = "blue", lty = 2)  # p-value threshold
abline(v = c(-1, 1), col = "red", lty = 2)       # fold change thresholds

# Subset upregulated and downregulated genes
up <- degs[degs$logFC > 1, ]
down <- degs[ degs$logFC < -1, ]

# Highlight DEGs
points(up$logFC, -log10(up$adj.P.Val), col="green", pch=20)
points(down$logFC, -log10(down$adj.P.Val), col="red", pch=20)

# UMAP plot (dimensionality reduction)
# Clean the expression data
ex <- na.omit(ex)                   # Remove missing values
ex <- ex[!duplicated(ex), ]         # Remove duplicate genes

# Run UMAP (transpose so samples are rows)
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)

df <- data.frame(UMAP1 = ump$layout[,1],
                 UMAP2 = ump$layout[,2],
                 Sample = rownames(ump$layout),
                 Group = gs)

# Customize UMAP with ggplot2
ggplot(df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 100) +
  labs(title = "UMAP Plot (n_neighbors = 15)",
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE47960")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE47960", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")