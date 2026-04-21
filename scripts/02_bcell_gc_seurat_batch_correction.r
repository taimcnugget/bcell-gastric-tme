# =============================================================
# Script: 02_bcell_gc_seurat_batch_correction.R
# Purpose: Batch correct cells that passed QC from the GSE163558 
#          scRNA-seq data from GEO. 
#
# Dataset: GSE163558 — gastric cancer TME scRNA-seq
#          10 samples: primary tumor (3), healthy (1), and
#          metastatic sites (liver (2), lymph node (2), ovary,
#          peritoneum)
#
# Input:   QC filtered Seurat object
# Output:  02_bcell_gc_seurat_correction.rds
#          figures/02_top_20_var_features.png
#          figures/02_pca_plot.png
#          figures/02_dim_red_genes.png
#          figures/02_UMAP_no_correction.png
#          figures/02_harmony_correction.png
#          figures/02_UMAP_harmony_correction.png
#
# Author:  Tailynn McCarty
# Date:    April 2026
# =============================================================

# --- Dependencies
# R packages required to run this script:
#
# CRAN: Seurat, Matrix, dplyr, ggplot2, patchwork
#
# Note: On Kaggle, Matrix, dplyr, ggplot2, and patchwork
# are pre-installed. Seurat and Harmony require installation
# each session (see install block below).

# --- Libraries
library(Seurat)
library(ggplot2)
library(harmony)

# Use Seurat v3 assay format for compatibility
options(Seurat.object.assay.version = "v3")

# Set seed for reproducibility
set.seed(3)

# --- Directories
output_directory <- "/kaggle/working"
figure_directory <- file.path(output_directory, "figures")

dir.create("figures")

#import qc file
results <- readRDS("/kaggle/input/datasets/taimcnugget/gastric-cancer-rds-files/b_cell_GC_TME_results_qc.rds")


# --- Parameters
feature_number <- 2000
n_pcs <- 20 #matches author's number of principle components

# --- Normalize Data
results <- NormalizeData(results)

results <- FindVariableFeatures(results, selection.method = "vst", nfeatures = feature_number)

# --- Investigate the top 20
top20 <- head(VariableFeatures(results), 20)

plot_1 <- VariableFeaturePlot(results)
top_20_viz <- LabelPoints(plot = plot_1, points = top20, repel = TRUE)

top_20_viz

ggsave(file.path(figure_directory, "02_top_20_var_features.png"),
       top_20_viz,
       width = 12,
       height = 5,
       dpi = 150)

# --- Scale Variable Genes
all_genes <- rownames (results)

results <- ScaleData(results, 
                     features = VariableFeatures(results),  # not all genes
                     vars.to.regress = c("percent.mt", "nCount_RNA"))

# --- PCA
results <- RunPCA(results, 
                 features = VariableFeatures(results),
                 npcs = 50, 
                 verbose = FALSE)

# Elbow plot — look for plateau in variance explained
ElbowPlot(results, ndims = 50)

# paper chose 20 PC, but it looks like the elbow is at 15 PCs

#vizualize
pca_plot <- DimPlot(results, reduction = 'pca')
dim_red_genes <- VizDimLoadings(results, dims = 1:n_pcs, reduction = 'pca')

ggsave(file.path(figure_directory, "02_pca_plot.png"),
       pca_plot,
       width = 12,
       height = 5,
       dpi = 150)

ggsave(file.path(figure_directory, "02_dim_red_genes.png"),
       plot = dim_red_genes,
       width = 12,
       height = 5,
       dpi = 150)


# --- UMAP reduction
# Before correction

results <- RunUMAP(results,
                   dims = 1:n_pcs,
                   verbose = FALSE)

# No correction UMAP
UMAP_no_correction <- DimPlot(results, 
        reduction = 'umap',
        group.by = "condition") + ggtitle("UMAP - no batch correction")


ggsave(file.path(figure_directory, "02_UMAP_no_correction.png"),
       UMAP_no_correction,
       width = 12,
       height = 5,
       dpi = 150)

# --- Harmony Correction
results_corrected <- RunHarmony(results,
                      group.by.vars = "orig.ident",
                      lambda = 1.5,
                      dims.use = 1:n_pcs,
                      verbose = FALSE)

# --- Investigate sample mixing after correction
harmony_correction <- DimPlot(results_corrected,
        reduction = 'harmony',
        group.by = "orig.ident",
        pt.size = 0.2) + ggtitle("Harmony — colored by sample")


ggsave(file.path(figure_directory, "02_harmony_correction.png"),
       harmony_correction,
       width = 12,
       height = 5,
       dpi = 150)


# --- UMAP reduction
# After correction
results_corrected <- RunUMAP(results_corrected,
                   reduction = "harmony",
                   dims = 1:n_pcs,
                   verbose = FALSE)

UMAP_harmony_correction <- DimPlot(results_corrected,
        reduction = 'umap',
        group.by = "condition") + ggtitle("UMAP - Harmony correction")


ggsave(file.path(figure_directory, "02_UMAP_harmony_correction.png"),
       UMAP_harmony_correction,
       width = 12,
       height = 5,
       dpi = 150)


# --- Compare correction to no correction
UMAP_no_correction/
UMAP_harmony_correction


# --- Save output
saveRDS(results_corrected, file.path(output_directory, "02_bcell_gc_seurat_correction.rds"))
message("Saved: 02_bcell_gc_seurat_correction.rds")