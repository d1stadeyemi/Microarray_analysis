#!/usr/bin/env Rscript
# script2_visualization.R
# ---------

# Install required packages
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("VennDiagram")
install.packages("patchwork")
install.packages("cowplot")
install.packages("dplyr")
install.packages("tidyr")
install.packages("grid")
install.packages("gridExtra")

# Load required packages for visualization
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(patchwork)
library(cowplot)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)


# Load key results
results <- readRDS("results/deg_results.rds")

degs_H1N1   <- results$degs_H1N1
degs_dORF6  <- results$degs_dORF6
degs_BatSRBD <- results$degs_BatSRBD
degs_icSARS <- results$degs_icSARS
fit2        <- results$fit2
ex          <- results$ex
gs          <- results$gs
gset        <- results$gset

# Create a venn diagram to visualize common DEGs
venn.plot <- venn.diagram(
  x = list(
    H1N1 = degs_H1N1$ID,
    dORF6 = degs_dORF6$ID,
    BatSRBD = degs_BatSRBD$ID,
    icSARS = degs_icSARS$ID
  ),
  filename = NULL,
  main = "Common DEGs across conditions",
  main.cex = 1.5,
  fill = c("red", "blue", "green", "purple"),
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "Helvetica",
  cat.cex = 1.2,
  cat.fontfamily = "Helvetica",
  cat.col = c("red", "blue", "green", "purple"),
)

png("venn_plot_common_DEGs.png", width = 2000, height = 2000, res = 300)
file.remove(list.files(pattern = "*.log")) # remove auto generated log file

grid.draw(venn.plot)
dev.off()

# Create histograms of Adjusted p-values for each condition
p1 <- ggplot(tT_H1N1, aes(x = adj.P.Val)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue", color = "white") +
  labs(title = "H1N1", x = "Adjusted P-value", y = "Number of genes") +
  theme_minimal(base_size = 14)

p2 <- ggplot(tT_dORF6, aes(x = adj.P.Val)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue", color = "white") +
  labs(title = "dORF6", x = "Adjusted P-value", y = "Number of genes") +
  theme_minimal(base_size = 14)

p3 <- ggplot(tT_BatSRBD, aes(x = adj.P.Val)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue", color = "white") +
  labs(title = "BatSRBD", x = "Adjusted P-value", y = "Number of genes") +
  theme_minimal(base_size = 14)

p4 <- ggplot(tT_icSARS, aes(x = adj.P.Val)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue", color = "white") +
  labs(title = "icSARS", x = "Adjusted P-value", y = "Number of genes") +
  theme_minimal(base_size = 14)

# Combine to a 2x2 grid and label A-D
final_plot <- plot_grid(
  p1, p2, p3, p4,
  labels = c("A", "B", "C", "D"),
  label_size = 18,
  ncol = 2,
  label_fontface = "bold"
)

# Export as png
ggsave("pval_histograms.png", final_plot, width = 12, height = 10, dpi = 300)

# Create a volcano plot
# Ensure required columns exist for the plot
# by taking the -log10 of the adjusted pvalues

# Function to create a single volcano plot without labels inside the plot
make_volcano <- function(tT, degs, title) {
  tT$logP <- -log10(tT$adj.P.Val)
  tT$category <- "NS"
  tT$category[tT$logFC > 1 & tT$adj.P.Val < 0.05] <- "Up"
  tT$category[tT$logFC < -1 & tT$adj.P.Val < 0.05] <- "Down"

  ggplot(tT, aes(x = logFC, y = logP, color = category)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Up" = "green", "Down" = "red", "NS" = "grey")) +
    geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "brown", linetype = "dashed") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(Adj p-value)") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
}

# Create each plot
p1 <- make_volcano(tT_H1N1, degs_H1N1, "H1N1")
p2 <- make_volcano(tT_dORF6, degs_dORF6, "dORF6")
p3 <- make_volcano(tT_BatSRBD, degs_BatSRBD, "BatSRBD")
p4 <- make_volcano(tT_icSARS, degs_icSARS, "icSARS")

# Combine into 2x2 grid with labels A, B, C, D
final_plot <- plot_grid(
  p1, p2, p3, p4,
  labels = c("A", "B", "C", "D"),
  label_size = 18,
  ncol = 2,
  label_fontface = "bold"
)

# Export as png
ggsave("volcano_plot_DEGs.png", final_plot, width = 12,
      height = 10, dpi = 300)

# Count upregulated and downregulated DEGs in each condition
deg_counts <- data.frame(
  Condition = c("H1N1", "dORF6", "BatSRBD", "icSARS"),
  Up = c(
    sum(degs_H1N1$logFC > 1),
    sum(degs_dORF6$logFC > 1),
    sum(degs_BatSRBD$logFC > 1),
    sum(degs_icSARS$logFC > 1)
  ),
  Down = c(
    sum(degs_H1N1$logFC < -1),
    sum(degs_dORF6$logFC < -1),
    sum(degs_BatSRBD$logFC < -1),
    sum(degs_icSARS$logFC < -1)
  )
)

# Add total DEGs column for sorting
deg_counts <- deg_counts %>%
  mutate(Total = Up + Down)

# Convert to long format for ggplot
deg_counts_long <- deg_counts %>%
  pivot_longer(cols = c("Down", "Up"), names_to = "Direction", values_to = "Count")

# Stacked bar plot with labels
final_plot_bar <- ggplot(deg_counts_long, aes(x = reorder(Condition, -Total), y = Count, fill = Direction)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(Count > 0, Count, "")),
            position = position_stack(vjust = 0.5), # center labels inside each bar segment
            color = "black",
            size = 4) +
  scale_fill_manual(values = c("Down" = "orange", "Up" = "green")) +
  labs(
    title = "Differentially Expressed Genes per Condition",
    x = "Condition",
    y = "Number of DEGs",
    fill = "Regulation"
  ) +
  theme_minimal(base_size = 12)

# Export as ong
ggsave("bar_plot_DEG_count.png", final_plot_bar, width = 6,
      height = 8, dpi = 300)


# Create a expression value distribution plot and Density plot
png("plotSA_Plotdensity.png", width = 3500, height = 2000, res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE47960")
mtext("A", side = 3, adj = 0, line = 0.5, font = 2, cex = 1.5)

# expression value distribution
title <- paste ("GSE47960", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")
mtext("B", side = 3, adj = 0, line = 0.5, font = 2, cex = 1.5)

dev.off()