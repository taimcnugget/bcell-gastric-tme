# =============================================================
# Script: 02_batch_correction.R
# Purpose: Load qc rds file, normalize and scale the data,
#          reduce dimensions using PCA, batch correct, and
#          save for further downstream analysis.
#
# Dataset: GSE163558 — gastric cancer TME scRNA-seq
#          10 samples: primary tumor (3), healthy (1), and
#           metastatic
#          sites (liver (2), lymph node (2), ovary, peritoneum)
#
# Input:   Raw data downloaded from GEO (GSE163558_RAW.tar)
# Output:  02_bcell_gc_seurat_batch_correction.rds
#          figures/02_top20_variable_features.png
#          figures/02_pca_elbow.png
#          figures/02_pca.png
#          figures/02_harmony_correction.png
#
# Author:  Tailynn McCarty
# Date:    March 2026
# ============================================================

# --- Dependencies
# R packages required to run this script:
#
# CRAN: Seurat, ggplot2, harmony
#
# Note: On Kaggle, ggplot2 is pre-installed.
# Seurat and harmony require installation
# each session (see install block below).
#
# Local install:
#   install.packages(c("Seurat",
#                      "ggplot2",
#                      "harmony"))

# --- Libraries
library(Seurat)
library(ggplot2)
library(harmony)

set.seed(3)
options(Seurat.object.assay.version = "v3")


# --- Paths
input <- "/kaggle/input/datasets/taimcnugget/gastric-cancer-rds-files/"
output_directory <- "/kaggle/working"

qc_file <- file.path(input, "01_bcell_gc_seurat_qc.rds")
figure_directory <- file.path(output_directory, "figures")

dir.create(figure_directory, showWarnings = FALSE, recursive = TRUE)

# load with existence check
if (file.exists(qc_file)) {
  results <- readRDS(qc_file)
} else {
  stop("File not found: ", qc_file)
}


# --- Normalize
results <- NormalizeData(qc_file)

results <- FindVariableFeatures(results,
                                selection.method = "vst", nfeatures = 2000)

# check top 20 variable features
top20 <- head(VariableFeatures(results), 20)
print(top20)



#vizualize top 20
plot_1 <- VariableFeaturePlot(results)
top_20_viz <- LabelPoints(plot = plot_1,
                          points = top20,
                          repel = TRUE)

top_20_viz

ggsave(file.path(figure_directory,
                 "02_top20_variable_features.png"),
       top_20_viz, width = 8, height = 6, dpi = 150)


# ---- Scaling (variable genes only)
message("Initiating Scaling...")
results <- ScaleData(results,
                     features = VariableFeatures(results),
                     vars.to.regress = c("percent.mt", "nCount_RNA"))
message("Scaling complete")


# --- PCA
message("Initiating PCA...")
results <- RunPCA(results,
                  features = VariableFeatures(results),
                  npcs = 50,
                  verbose = FALSE)
message("PCA complete")

# Elbow plot — look for plateau in variance explained
elbow <- ElbowPlot(results, ndims = 50)

ggsave(file.path(figure_directory,
                 "02_pca_elbow.png"),
       elbow, width = 8, height = 6, dpi = 150)

# paper chose 20 PC, elbow looks closer to 15 but staying consistent with paper

n_pcs <- 20

# vizualize PCA
pca_plot <- DimPlot(results, reduction = "pca")
VizDimLoadings(results, dims = 1:n_pcs, reduction = "pca")

ggsave(file.path(figure_directory,
                 "02_pca.png"),
       pca_plot, width = 8, height = 6, dpi = 150)

# --- Batch correction, Harmony
# Correcting for donor/organ heterogeneity across patients and sample sites
message("Initiating Harmony batch correction...")
results <- RunHarmony(results,
                      group.by.vars = "orig.ident",
                      dims.use = 1:n_pcs,
                      verbose = FALSE)
message("Harmony batch correction complete...")

# Quick check — sample mixing after correction
DimPlot(results,
        reduction = "harmony",
        group.by = "orig.ident",
        pt.size = 0.2) + ggtitle("Harmony — colored by sample")

ggsave(file.path(figure_directory,
                 "02_harmony_correction.png"),
       pca_plot, width = 8, height = 6, dpi = 150)


# --- Save output
saveRDS(results, file.path("results/02_bcell_gc_seurat_batch_correction.rds"))