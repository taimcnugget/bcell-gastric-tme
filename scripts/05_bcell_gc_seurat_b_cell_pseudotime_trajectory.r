# =============================================================
# Script: 06_bcell_gc_seurat_b_cell_trajectory_analysis.R
# Purpose: Investigate the trajectory of b cell differentiation in
#          healthy, primary TME, and metastatic TME to investigate
#          if there are certain subsets that remain longer in the TMEs
#          and if so, what are the genetic signatures of these cells.
#          These cells can then be further investigated in their role
#          in tumor progression. 
#
# Dataset: GSE163558 — gastric cancer TME scRNA-seq
#          10 samples: primary tumor (3), healthy (1), and
#          metastatic sites (liver (2), lymph node (2), ovary,
#          peritoneum)
#
# Input:   B cell subset annotated Seurat Object (labeled as b_cells)
# Output:  06_bcell_gc_seurat_b_cell_trajectory_analysis.rds
#          06_b_cell_trajectory.png
#          06_lineage_b_cell_type.png
#          06_pseudotime_umap.png
#
# Author:  Tailynn McCarty
# Date:    April 2026
# =============================================================

# --- Dependencies
# R packages required to run this script:
#
# CRAN: Seurat, Matrix, dplyr, ggplot2, patchwork, slingshot, 
#       tradeSeq, SingleCellExperiment, IRanges, DelayedMatrixStates
#
# Note: On Kaggle, Matrix, dplyr, ggplot2, and patchwork
# are pre-installed. Others require installation
# each session (see install block below).

# --- Directories
output_directory <- "/kaggle/working"
figure_directory <- file.path(output_directory, "figures")

dir.create("figures")

# ---- Packages
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install(c("slingshot", 
                       "tradeSeq", 
                       "SingleCellExperiment",
                       "IRanges",
                       "Seurat", 
                       "DelayedMatrixStats"))

# --- Libraries
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(Seurat)
library(BiocParallel)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# Use Seurat v3 assay format for compatibility
options(Seurat.object.assay.version = "v3")


# --- Import data
b_cells <- readRDS("/kaggle/input/datasets/taimcnugget/bcell-gc-tme-rds-updated/05_bcell_gc_seurat_b_cell_deg_analysis.rds")

# --- Convert object to singlecellexperiment object
b_cell_sce <- Seurat:::as.SingleCellExperiment.Seurat(b_cells)

#check
b_cell_sce

#check dimension reductions
reducedDimNames(b_cell_sce)
table(b_cell_sce$b_cell_type)

# --- Slingshot
set.seed(3)
b_cell_sce <- slingshot(b_cell_sce, 
                        clusterLabels = "cell_type", 
                        reducedDim = "HARMONY", 
                        start.clus = "Naive B cells", end.clus = c("Memory B cells", "Plasma cells"))

slingLineages(b_cell_sce)

# --- Replace N/A values with 0 for plotting
pt <- slingPseudotime(b_cell_sce)
pt[is.na(pt)] <- 0


# --- Add lineages to Seurat object
b_cells$lineage_1 <- pt[, 1]
b_cells$lineage_2 <- pt[, 2]
b_cells$lineage_3 <- pt[, 3]
b_cells$lineage_4 <- pt[, 4]


# --- Slingshot curves
umap_coords <- Embeddings(b_cells, "harmony")

png(filename = file.path(figure_directory, "06_b_cell_trajectory.png"),
    width = 800,
    height = 800,
    res = 150
)

plot(umap_coords,
     col = RColorBrewer::brewer.pal(9, "Set1")[factor(b_cells$cell_type)], 
     pch = 16, 
     cex = 0.3, main = "B Cell Trajectory")

lines(SlingshotDataSet(b_cell_sce), 
      lwd = 2, 
      col = "black")

dev.off()

#comparing pseudotime across conditions

b_cells$condition_grouped <- factor(b_cells$condition_grouped, 
                                    levels = c("Healthy", "Primary", "Metastatic")
)

lineage_b_cell_type <- VlnPlot(b_cells,
        features = c("lineage_1", "lineage_2", "lineage_3", "lineage_4"),
        group.by = "condition_grouped",
        pt.size = 0, 
        combine = TRUE, 
        ncol = 2)

lineage_b_cell_type

ggsave(file.path(figure_directory, "06_lineage_b_cell_type.png"),
       lineage_b_cell_type,
       width = 8,
       height = 8,
       dpi = 150)

# --- Retrieve data
counts_matrix <- GetAssayData(b_cells,
                              assay = "RNA", 
                              layer = "counts")

var_genes <- VariableFeatures(b_cells)
counts_matrix <- counts_matrix[var_genes, ]

pseudotime <- slingPseudotime(b_cell_sce)

weights <-slingCurveWeights(b_cell_sce)

head(pseudotime)

# --- Remove cells that do not fit into any lineage
pseudotime_clean <- slingPseudotime(b_cell_sce)
pseudotime_clean[is.na(pseudotime_clean)] <- 0

weights_clean <- slingCurveWeights(b_cell_sce)

all_zero <- apply(pseudotime_clean, 1, function(x) all(x == 0))
cat("Cells removed (no lineage):", sum(all_zero), "\n")

pseudotime_clean <- pseudotime_clean[!all_zero, ]
weights_clean    <- weights_clean[!all_zero, ]
counts_clean     <- counts_matrix[, !all_zero]

# check dimensions
cat("Pseudotime dims:", dim(pseudotime_clean), "\n")
cat("Weights dims:",    dim(weights_clean), "\n")
cat("Counts dims:",     dim(counts_clean), "\n")

#check removal
cat("NAs remaining:", sum(is.na(pseudotime_clean)), "\n")

# --- Note: subsetting random cells to make this run faster
set.seed(3)
top_genes <- VariableFeatures(b_cells)[1:1000]
counts_clean_filtered <- counts_clean[rownames(counts_clean) %in% top_genes, ]

b_cell_sample <- sample(ncol(counts_clean_filtered), 3000)
counts_clean_small <- counts_clean_filtered[, b_cell_sample]
pseudotime_clean   <- pseudotime_clean[b_cell_sample, ]
weights_clean      <- weights_clean[b_cell_sample, ]

# Verify changes
cat("Counts:", dim(counts_clean_small), "\n")
cat("Pseudotime:", dim(pseudotime_clean), "\n")
cat("Weights:", dim(weights_clean), "\n")

# --- Set up parallel processing to speed up fitGAM
BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))

# --- Fit the model (may take time to run)
b_cell_sce_tradeseq <- fitGAM(counts = counts_clean_small, 
                              pseudotime = pseudotime_clean, 
                              cellWeights = weights_clean, 
                              nknots = 4, 
                              verbose = TRUE, 
                              parallel = TRUE)


# --- Save immediately in case of a crash 
saveRDS(b_cell_sce_tradeseq, file.path(output_directory,"06_bcell_gc_seurat_b_cell_trajectory_analysis.rds"))
message("Saved: 06_bcell_gc_seurat_b_cell_trajectory_analysis.rds")


# --- Pseudotime UMAP - one plot per lineage
lineage_plots <- list()

for (i in 1:4) {
  lineage_col <- paste0("lineage_", i)
  
  lineage_plots[[i]] <- FeaturePlot(b_cells,
                                     features = lineage_col,
                                     reduction = "umap",
                                     pt.size = 0.3) +
    scale_color_viridis_c(option = "inferno") +
    ggtitle(paste("Pseudotime - Lineage", i)) +
    theme_classic()
}

# --- Combine into one figure
pseudotime_umap <- wrap_plots(lineage_plots, ncol = 2)

pseudotime_umap

ggsave(file.path(figure_directory, "06_pseudotime_umap.png"),
       pseudotime_umap,
       width = 12,
       height = 10,
       dpi = 150)
message("Saved: 06_pseudotime_umap.png")

# --- Association test - genes that differ along any lineage
assoc_test <- associationTest(b_cell_sce_tradeseq)
write.csv(assoc_test, 
          file.path(output_directory, "06_tradeseq_association.csv"))

# --- Diff end test - genes that differ between lineage endpoints
end_test <- diffEndTest(b_cell_sce_tradeseq, pairwise = TRUE)
write.csv(end_test,
          file.path(output_directory, "06_tradeseq_diffend.csv"))


# --- Save slingshot analysis
saveRDS(b_cell_sce, file.path(output_directory,"06_bcell_gc_seurat_b_cell_sce_pseudotime.rds"))
message("Saved: 06_bcell_gc_seurat_b_cell_sce_pseudotime.rds")


# --- Save pseudotime stats across all 4 lineages
sink(file.path(output_directory, "pseudotime_stats.txt"))

for (i in 1:4) {
  lineage_col <- paste0("lineage_", i)
  cat(paste0("=== LINEAGE ", i, " ===\n"))
  
  pseudotime_vals <- b_cells@meta.data[[lineage_col]]
  condition_vals  <- b_cells$condition_grouped
  
  # remove NAs
  keep <- !is.na(pseudotime_vals)
  pseudotime_vals <- pseudotime_vals[keep]
  condition_vals  <- condition_vals[keep]
  
  print(kruskal.test(pseudotime_vals ~ condition_vals))
  print(pairwise.wilcox.test(pseudotime_vals,
                              condition_vals,
                              p.adjust.method = "BH"))
  cat("\n")
}

sink()
message("Saved: 06_pseudotime_stats.txt")

# --- Save pseudotime values
pseudotime_df <- data.frame(
  cell      = colnames(b_cells),
  lineage_1 = b_cells$lineage_1,
  lineage_2 = b_cells$lineage_2,
  lineage_3 = b_cells$lineage_3,
  lineage_4 = b_cells$lineage_4,
  condition = b_cells$condition_grouped,
  cell_type = b_cells$cell_type
)

write.csv(pseudotime_df,
          file.path(output_directory, "06_pseudotime_values.csv"),
          row.names = FALSE)
message("Saved: 06_pseudotime_values.csv")