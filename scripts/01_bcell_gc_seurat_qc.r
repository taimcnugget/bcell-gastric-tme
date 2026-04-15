# =============================================================
# Script: 01_load_and_qc.R
# Purpose: Load GSE163558 scRNA-seq data from GEO, merge into
#          a single Seurat object, assign condition metadata,
#          perform QC filtering, and save for downstream use.
#
# Dataset: GSE163558 — gastric cancer TME scRNA-seq
#          10 samples: primary tumor (3), healthy (1), and
#          metastatic sites (liver (2), lymph node (2), ovary,
#          peritoneum)
#
# Input:   Raw data downloaded from GEO (GSE163558_RAW.tar)
# Output:  01_bcell_gc_seurat_qc.rds
#          figures/01_vln_qc_metrics.png
#          figures/01_scatter_umi_mt.png
#          figures/01_scatter_umi_genes.png
#
# Author:  Tailynn McCarty
# Date:    March 2026
# =============================================================

# --- Dependencies
# R packages required to run this script:
#
# CRAN: Seurat, Matrix, dplyr, ggplot2, patchwork
# Bioconductor: GEOquery
#
# Note: On Kaggle, Matrix, dplyr, ggplot2, and patchwork
# are pre-installed. Seurat and GEOquery require installation
# each session (see install block below).
#
# Local install:
#   install.packages(c("Seurat", "Matrix", "dplyr",
#                      "ggplot2", "patchwork"))
#   BiocManager::install("GEOquery")

# --- Libraries

library(Seurat)
library(GEOquery)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)

# Use Seurat v3 assay format for compatibility
options(Seurat.object.assay.version = "v3")

# Set seed for reproducibility
set.seed(3)

# --- Parameters

data_directory <- "/kaggle/working/GSE163558"
extracted_directory <- file.path(data_directory, "extracted")
output_directory <- "/kaggle/working"
figure_directory <- file.path(output_directory, "figures")

# QC thresholds — matched to original authors' parameters
min_features <- 200
max_features <- 5000
max_mt_pct <- 20

# Seurat object creation thresholds
min_cells   <- 3
dir.create(figure_directory, showWarnings = FALSE, recursive = TRUE)



# --- Condition mapping
# Maps GEO sample IDs to human-readable condition labels
condition_map <- c(
  "GSM5004180" = "Primary_1",
  "GSM5004181" = "Primary_2",
  "GSM5004182" = "Primary_3",
  "GSM5004183" = "Healthy",
  "GSM5004184" = "Lymph_1",
  "GSM5004185" = "Lymph_2",
  "GSM5004186" = "Ovary",
  "GSM5004187" = "Peritoneum",
  "GSM5004188" = "Liver_1",
  "GSM5004189" = "Liver_2"
)

# --- Download & extract data
dir.create(data_directory, showWarnings = FALSE)

getGEOSuppFiles("GSE163558", makeDirectory = TRUE,
                baseDir = "/kaggle/working/")

untar(file.path(data_directory, "GSE163558_RAW.tar"),
      exdir = extracted_directory)



# --- Update all list.files() calls to use new path
all_files <- list.files(extracted_directory,
                        full.names = FALSE)


# --- load data into Seurat
gsm_ids <- unique(gsub("_.*", "", all_files))
print(gsm_ids)


seurat_list <- lapply(gsm_ids, function(gsm) {
  mat_file <- list.files(extracted_directory,
                         pattern = paste0(gsm, ".*matrix\\.mtx\\.gz$"),
                         full.names = TRUE)
  bar_file <- list.files(extracted_directory,
                         pattern = paste0(gsm, ".*barcodes\\.tsv\\.gz$"),
                         full.names = TRUE)
  feat_file <- list.files(extracted_directory,
                          pattern = paste0(gsm, ".*features\\.tsv\\.gz|genes\\.tsv\\.gz$"),
                          full.names = TRUE)

  counts   <- readMM(mat_file)
  barcodes <- read.delim(bar_file, header = FALSE)
  features <- read.delim(feat_file, header = FALSE)

  rownames(counts) <- make.unique(features$V2)  # gene symbols
  colnames(counts) <- barcodes$V1

  CreateSeuratObject(counts, 
                     project = gsm,
                     min.cells = min_cells, 
                     min.features = min_features)
})


# Name the list by GSM ID
names(seurat_list) <- gsm_ids


# --- Merge into a single object
results <- Reduce(function(a,b) merge(a,b), seurat_list)
message("Cells before QC: ", ncol(results))
print(table(results$orig.ident))



# --- Create a named vector of conditions matched to cell names

condition_vector <- condition_map[as.character(results$orig.ident)]
names(condition_vector) <- colnames(results)


results <- AddMetaData(results, 
                       metadata = condition_vector,
                       col.name = "condition")

# Verify
print(table(results$condition))


# --- QC metrics
# Mitochondrial percentage - high MT% indicates dying/damaged cells
# (will have MT- as prefix)
results[["percent.mt"]] <- PercentageFeatureSet(results, pattern = "^MT-")

message("QC metric summary (pre-filtering):")
print(summary(results@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]))


# --- QC Vizualizations

# Violin plot
plt_vln <- VlnPlot(results,
                   features = c("nFeature_RNA",
                                "nCount_RNA",
                                "percent.mt"),
                   group.by = "condition",
                   ncol = 3,
                   pt.size = 0)

ggsave(file.path(figure_directory, "01_vln_qc_metrics.png"),
       plt_vln,
       width = 12,
       height = 5,
       dpi = 150)

# Scatter plots
# UMI count vs MT% — dying cells: low UMI + high MT%
plt_umi_mt <- FeatureScatter(results,
                             feature1 = "nCount_RNA",
                             feature2 = "percent.mt",
                             group.by = "condition") +
  ggtitle("UMI Count vs Mitochondrial %") +
  xlab("UMI Count") +
  ylab("Mitochondrial %")

ggsave(file.path(figure_directory, "01_scatter_umi_mt.png"),
       plt_umi_mt, width = 8, height = 6, dpi = 150)

#UMI count vs gene count - doublets: high UMI + high gene count
plt_umi_genes <- FeatureScatter(results,
                                feature1 = "nCount_RNA",
                                feature2 = "nFeature_RNA",
                                group.by = "condition") +
  ggtitle("UMI Count vs Gene Count") +
  xlab("UMI Count") +
  ylab("Gene Count")

ggsave(file.path(figure_directory, "01_scatter_umi_genes.png"),
       plt_umi_genes, width = 8, height = 6, dpi = 150)

# --- QC filtering
# Thresholds matched to original authors' QC parameters
# Removes: low-quality cells, potential doublets,and dying cells (high MT%)

results <- subset(results,
                  subset = nFeature_RNA > min_features &
                    nFeature_RNA < max_features &
                    percent.mt < max_mt_pct)


message("Cells after QC filtering:", ncol(results))
print((table(results$condition)))


# --- Save output

saveRDS(results, file.path(output_directory, "01_bcell_gc_seurat_qc.rds"))

message("Saved: results/01_bcell_gc_seurat_qc.rds")