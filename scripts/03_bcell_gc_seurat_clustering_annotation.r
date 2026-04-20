# =============================================================
# Script: 03_bcell_gc_seurat_clustering_annotation.R
# Purpose: Clustering and annotating cell type clusters from
#          batch corrected variable gene counts. This step is
#          important for down stream analysis to ensure cell
#          types are accurate.
#
# Dataset: GSE163558 — gastric cancer TME scRNA-seq
#          10 samples: primary tumor (3), healthy (1), and
#          metastatic sites (liver (2), lymph node (2), ovary,
#          peritoneum)
#
# Input:   Harmony Batch corrected Seurat Object (labeled as results)
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
# CRAN: Seurat, Matrix, dplyr, ggplot2, patchwork, glue
#
# Note: On Kaggle, Matrix, dplyr, ggplot2, and patchwork
# are pre-installed. Seurat and glue requires installation
# each session (see install block below).

# --- Libraries
library(Seurat)
library(ggplot2)
library(harmony)
library(glue)

# Use Seurat v3 assay format for compatibility
options(Seurat.object.assay.version = "v3")

# Set seed for reproducibility
set.seed(3)

# --- Directories
output_directory <- "/kaggle/working"
figure_directory <- file.path(output_directory, "figures")

dir.create("figures")

# --- Import corrected object
results <- readRDS("/kaggle/input/datasets/taimcnugget/bcell-gc-tme-rds-updated/02_bcell_gc_seurat_correction.rds")


# --- Parameters
n_pcs <- 20 #matches author's number of principle components
knn_param <- 20
louvain_res <- 0.3 #resolution for clustering (Louvain algo.)

# --- Clustering (KNN)
results <- FindNeighbors(results,
                         reduction = "harmony",
                         dims = 1:n_pcs,
                         k.param = knn_param)
results <- FindClusters(results,
                        resolution = louvain_res)
table(Idents(results))

# --- Cluster UMAP
knn_cluster_UMAP <- DimPlot(results, 
                              reduction = 'umap', 
                              label = TRUE,
                              repel = TRUE, 
                              pt.size = 0.3) + ggtitle(glue("Louvain Clusters (res = {louvain_res})"))


ggsave(file.path(figure_directory, "03_knn_cluster_UMAP.png"),
       knn_cluster_UMAP,
       width = 12,
       height = 5,
       dpi = 150)

# --- Join layers
results <- JoinLayers(results)

# --- Find the marker genes for each cluster
markers <- FindAllMarkers(results, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              test.use = "wilcox")


top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  summarise(top_genes = paste(gene, collapse = ", "),
            top_logFC = paste(round(avg_log2FC, 2), collapse = ", "))
print(top5)

# --- Feature Plot
FeaturePlot(results,
            features = c("MAL",       # 0 - CD4+ T cells
                         "GZMH",      # 1 - CD8+ T cells
                         "IGHD",      # 2 - Naive B cells
                         "CXCR2",     # 3 - Neutrophils
                         "FOXP3",     # 4 - Tregs
                         "IGHG1",     # 5 - Plasma cells
                         "TFF1",      # 6 - Epithelial cells
                         "CD163",     # 7 - Macrophages
                         "SH2D1B",    # 8 - NK cells
                         "INTS6",     # 9 - Unknown
                         "COL14A1",   # 10 - Fibroblasts
                         "SPC25",     # 11 - Proliferating cells
                         "KDR",       # 12 - Endothelial cells
                         "TCL1A",     # 13 - Transitional B cells
                         "TPSAB1",    # 14 - Mast cells
                         "SERPINA9"), # 15 - Germinal center B cells
            ncol = 3, pt.size = 0.1)


# --- Investigate the "Unknown" cluster

VlnPlot(results, 
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        idents = "9")

# --- High number of low quality cells - must have slipped through QC (technially they hit the
# --- minimum threshold we set, but they are still skewing results so they're getting removed)

results <- subset(results, idents = "9", invert = TRUE)
results@active.ident <- droplevels(results@active.ident)

# confirm 9 is gone
levels(results)

# --- Annotate UMAP
cell_type_labels <- c("0"  = "CD4+ T cells",
                      "1"  = "CD8+ T cells",
                      "2"  = "Naive B cells",
                      "3"  = "Neutrophils",
                      "4"  = "Tregs",
                      "5"  = "Plasma cells",
                      "6"  = "Epithelial cells",
                      "7"  = "Macrophages",
                      "8"  = "NK cells",
                      "10" = "Fibroblasts",
                      "11" = "Proliferating cells",
                      "12" = "Endothelial cells",
                      "13" = "Transitional B cells",
                      "14" = "Mast cells",
                      "15" = "GC B cells")

results <- RenameIdents(results, cell_type_labels)

results$cell_type <- Idents(results)


annotated_UMAP <- DimPlot(results,
                          reduction = "umap",
                          label = TRUE,
                          repel = TRUE,
                          pt.size = 0.2) + ggtitle("") + NoLegend()

annotated_UMAP


ggsave(file.path(figure_directory, "03_annotated_UMAP.png"),
       annotated_UMAP,
       width = 12,
       height = 5,
       dpi = 150)

# --- Save file
saveRDS(results, file.path(output_directory, "03_bcell_gc_seurat_clustering_annotation.rds"))
message("Saved: 03_bcell_gc_seurat_clustering_annotation.rds")