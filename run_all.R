#!/usr/bin/env Rscript
# run_all.R
# ---------
# Run DEG analysis & visualization in sequence

cat("🔹 Running DEG analysis...\n")
source("Scripts/script1_DEG_analysis.R")

cat("🔹 Running visualizations...\n")
source("Scripts/script2_visualization.R")

cat("✅ All analysis completed. Results saved in results/ folder.\n")