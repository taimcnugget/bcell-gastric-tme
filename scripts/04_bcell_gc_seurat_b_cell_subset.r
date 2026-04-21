# =============================================================
# Script: 04_bcell_gc_seurat_b_cell_subset.R
# Purpose: Subset B cells from the gastric cancer dataset to
#          investigate whether there are trends in which B cells
#          within a TME, which can be used for immunotherapeutic
#          treatment research. 
#
# Dataset: GSE163558 — gastric cancer TME scRNA-seq
#          10 samples: primary tumor (3), healthy (1), and
#          metastatic sites (liver (2), lymph node (2), ovary,
#          peritoneum)
#
# Input:   Annotated Seurat Object (labeled as results)
# Output:  03_bcell_gc_seurat_clustering_annotation.rds
#          figures/03_annotated_UMAP.png
#          figures/03_knn_cluster_UMAP.png
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
# are pre-installed. Seurat requires installation
# each session (see install block below).

# --- Directories
output_directory <- "/kaggle/working"
figure_directory <- file.path(output_directory, "figures")

dir.create("figures")

# --- Libraries
library(Seurat)
library(ggplot2)
library(harmony)
library(patchwork)
library(dplyr)

# Use Seurat v3 assay format for compatibility
options(Seurat.object.assay.version = "v3")

# Set seed for reproducibility
set.seed(3)

# --- Load data
results <- readRDS("/kaggle/input/datasets/taimcnugget/bcell-gc-tme-rds-updated/03_bcell_gc_seurat_clustering_annotation.rds")

# --- Parameters
feature_number <- 2000
n_pcs <- 10
knn_param <- 20

# --- Subset B cells
b_cells <- subset(results, idents = c("Naive B cells",
                                       "GC B cells", 
                                       "Transitional B cells",
                                       "Plasma cells"))

#check
message(paste("There are", ncol(b_cells), "cells in the subset"))
table(Idents(b_cells))

top_features <- VlnPlot(b_cells, features = c("CD19",
                              "IGHD", 
                              "SERPINA9", 
                              "TCL1A", 
                              "IGHG1",
                              "CD38"),
        group.by = 'cell_type', pt.size = 0) + plot_annotation(title = "Top B Cell Features")

ggsave(file.path(figure_directory, "04_top_features.png"),
       top_features,
       width = 10,
       height = 10,
       dpi = 150)

# --- Rerun pipeline on B cells only
# --- Normalize
b_cells <- NormalizeData(b_cells)

# -- Find variable features
b_cells <- FindVariableFeatures(b_cells, 
                                selection.method = 'vst', 
                                nfeatures = feature_number)
# -- Scale B cells
b_cells <- ScaleData(b_cells,
                     features = VariableFeatures(b_cells),
                     vars.to.regress = c("percent.mt",
                                         "nCount_RNA"))

# --- Dimensional reduction (PCA)
b_cells <- RunPCA(b_cells, 
                  npcs = 50,
                  features = VariableFeatures(b_cells),
                  verbose = FALSE)

ElbowPlot(b_cells, ndims = 50)

pca_plot <- DimPlot(b_cells, reduction = 'pca')
dim_red_genes <- VizDimLoadings(b_cells, dims = 1:n_pcs, reduction = 'pca')

ggsave(file.path(figure_directory, "04_pca_plot.png"),
       pca_plot,
       width = 8,
       height = 5,
       dpi = 150)

ggsave(file.path(figure_directory, "04_dim_red_genes.png"),
       plot = dim_red_genes,
       width = 8,
       height = 5,
       dpi = 150)

# --- Harmony Correction
b_cells <- RunHarmony(b_cells,
                      group.by.vars = "orig.ident",
                      lambda = 1.5,
                      dims.use = 1:n_pcs,
                      verbose = FALSE)

b_cell_harmony_correction <- DimPlot(b_cells,
        reduction = 'harmony',
        group.by = "condition",
        pt.size = 0.2) + ggtitle("Harmony — colored by sample")


ggsave(file.path(figure_directory, "04_b_cell_harmony_correction.png"),
       b_cell_harmony_correction,
       width = 8,
       height = 5,
       dpi = 150)

# --- Dimensional reduction (UMAP)
b_cells <- RunUMAP(b_cells, 
                   reduction = "harmony", 
                   dims = 1:n_pcs)

# --- Cluster

b_cells <- FindNeighbors(b_cells, 
                         reduction = "harmony", 
                         dims = 1:n_pcs, 
                         k.param = knn_param)

b_cells <- FindClusters(b_cells, resolution = 0.3)


b_cell_UMAP <- DimPlot(b_cells, 
                       reduction = "umap",
                       label = TRUE, repel = TRUE) + ggtitle("B Cell Subtypes")

ggsave(file.path(figure_directory, "04_b_cell_UMAP.png"),
       b_cell_UMAP,
       width = 8,
       height = 5,
       dpi = 150)

# --- Find markers
b_cells <- JoinLayers(b_cells)

b_cell_markers <- FindAllMarkers(b_cells, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

top5 <- b_cell_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  summarise(top_genes = paste(gene, collapse = ", "),
            top_logFC = paste(round(avg_log2FC, 2), collapse = ", "))

print(top5)

# --- Cluster investigation (clusters 2 and 5 seem similar and 0 and 4 seem similar)
markers_2v5 <- FindMarkers(b_cells,
                           ident.1 = "2",
                           ident.2 = "5",
                           min.pct = 0.25)
head(markers_2v5, 10)

# compare cluster 0 vs 4
markers_0v4 <- FindMarkers(b_cells,
                           ident.1 = "0",
                           ident.2 = "4",
                           min.pct = 0.25)
head(markers_0v4, 10)

# --- Remove contaminating clusters (3, 7, 8)

b_cells <- subset(b_cells, idents = c("3","7", "8"), invert = TRUE)
b_cells@active.ident <- droplevels(b_cells@active.ident)

# confirm they're gone
levels(b_cells)

# --- Annotate clusters
cell_type_labels <- c("0" = "Memory B cells",
                      "1" = "Plasma cells",
                      "2" = "Activated Naive B cells",
                      "4" = "Activated Memory B cells",
                      "5" = "Naive B cells",
                      "6" = "GC B cells")

b_cells <- RenameIdents(b_cells, cell_type_labels)
b_cells$cell_type <- Idents(b_cells)

b_cells$cell_type <- Idents(b_cells)


b_cell_annotated_UMAP <- DimPlot(b_cells, 
        reduction = "umap",
        label = TRUE, 
        repel = TRUE,
        pt.size = 0.2) + ggtitle("") + NoLegend()

b_cell_annotated_UMAP


ggsave(file.path(figure_directory, "04_b_cell_annotated_UMAP.png"),
       b_cell_annotated_UMAP,
       width = 18,
       height = 5,
       dpi = 150)