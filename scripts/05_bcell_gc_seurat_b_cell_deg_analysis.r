# =============================================================
# Script: 05_bcell_gc_seurat_b_cell_deg_analysis.R
# Purpose: Investigate DEG within the B cell subset B within
#          healthy, primary TME, and metastatic TME to investigate
#          differentially expressed genes, which can be used for
#          immunotherapeutic treatment research.
#
# Dataset: GSE163558 — gastric cancer TME scRNA-seq
#          10 samples: primary tumor (3), healthy (1), and
#          metastatic sites (liver (2), lymph node (2), ovary,
#          peritoneum)
#
# Input:   B cell subset annotated Seurat Object (labeled as b_cell)
# Output:  05_bcell_gc_seurat_b_cell_deg_analysis.rds
#          05_b_cell_TME_proportions.png
#          05_sig_DEG_comparison.png
#          05_memory_healthy_mets_deg.png
#          05_memory_healthy_primary_deg.png
#          05_memory_primary_mets_deg.png
#          05_plasma_healthy_mets_deg.png
#          05_plasma_healthy_primary_deg.png
#          05_plasma_primary_mets_deg.png
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
library(patchwork)
library(ggrepel)
library(dplyr)

# Use Seurat v3 assay format for compatibility
options(Seurat.object.assay.version = "v3")

# --- Set seed for reproducibility
set.seed(3)

# --- Load data
b_cells <- readRDS("/kaggle/input/datasets/taimcnugget/bcell-gc-tme-rds-updated/04_bcell_gc_seurat_b_cell_subset.rds")

# --- Parameters
p_adj_threshold <- 0.05
logfc_threshold <- 0.5
min_pct <- 0.1
min_cells <- 10

# ---- Calculate proportions of B cell types within each condition
b_cell_proportions <- as.data.frame(
  prop.table(table(b_cells$cell_type, 
                   b_cells$condition), margin = 2)
)

# Rename columns
colnames(b_cell_proportions) <- c("B_Cell_Type", "Condition", "Proportion")

# Plot
b_cell_TME_proportions <- ggplot(b_cell_proportions, 
       aes(x = Condition, y = Proportion, fill = B_Cell_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) +
  labs(title = "B Cell Subtype Proportions by Condition",
       x = "Condition",
       y = "Proportion",
       fill = "B Cell Type")

b_cell_TME_proportions

ggsave(file.path(figure_directory, "05_b_cell_TME_proportions.png"),
       b_cell_TME_proportions,
       width = 8,
       height = 5,
       dpi = 150)


# --- Group conditions to get a high level comparison between B cell types in healthy, primary, and metastatic TME
b_cells$condition_grouped <- ifelse(
  grepl("Primary", b_cells$condition), "Primary",
  ifelse(grepl("Lymph|Ovary|Peritoneum|Liver", b_cells$condition), 
         "Metastatic", "Healthy"))

table(b_cells$condition_grouped)


subtypes <- unique(b_cells$cell_type)

comparisons <- list(c("Primary", "Healthy"),
                    c("Metastatic", "Healthy"),
                    c("Metastatic", "Primary"))


# --- Differentially Expressed Genes loop
library(dplyr)
de_results <- list()

for (subtype in subtypes) {
  b_sub <- subset(b_cells, cell_type == subtype)
  Idents(b_sub) <- b_sub$condition_grouped
  
  for (comp in comparisons) {
    cond1 <- comp[1]
    cond2 <- comp[2]
    
    n1 <- sum(b_sub$condition_grouped == cond1)
    n2 <- sum(b_sub$condition_grouped == cond2)
    
    if (n1 < min_cells | n2 < min_cells) {
      cat("Skipping:", subtype, cond1, "vs", cond2,
          "- not enough cells\n")
      next
    }
    
    de <- tryCatch(
      expr = {
        FindMarkers(b_sub,
                    ident.1 = cond1,
                    ident.2 = cond2,
                    min.pct = min_pct,
                    test.use = "wilcox")
      },
      error = function(e) {
        cat("Error in:", subtype, cond1, "vs", cond2, "\n")
        return(NULL)
      }
    )
    
    if (!is.null(de)) {
      de$subtype    <- subtype
      de$condition1 <- cond1
      de$condition2 <- cond2
      de$gene       <- rownames(de)
      key <- paste(subtype, cond1, "vs", cond2, sep = "_")
      de_results[[key]] <- de
    }
  }
}

all_de <- bind_rows(de_results)

# --- FDR correction for multiple comparisons
all_de$p_val_adj_global <- p.adjust(all_de$p_val_adj, method = "BH")

# --- Filter for only significant DE genes
all_de_sig <- all_de %>%
  filter(p_val_adj_global < p_adj_threshold,
         abs(avg_log2FC) > logfc_threshold)

# Summary
message(paste("Total DE genes:", nrow(all_de_sig), "\n"))
table(all_de_sig$subtype, all_de_sig$condition1)

# --- Count significant genes per subtype per comparison
summary_de <- all_de_sig %>%
  group_by(subtype, condition1) %>%
  summarise(n_genes = n(), .groups = "drop")


# --- Vizualize significant genes per subtype per comparison
sig_DEG_comparison <- ggplot(summary_de, 
       aes(x = condition1, y = subtype, fill = n_genes)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_genes), size = 4) +
  scale_fill_gradient(low = "lightblue", high = "#1F3864") +
  theme_classic() +
  labs(title = "Number of Significant DE Genes",
       x = "Comparison", y = "B Cell Subtype",
       fill = "# DE Genes")

sig_DEG_comparison

ggsave(file.path(figure_directory, "05_sig_DEG_comparison.png"),
       sig_DEG_comparison,
       width = 8,
       height = 5,
       dpi = 150)

# --- Function to make volcano plot for every subtype/comparison
make_volcano <- function(subtype_name, cond1, cond2) {
  df <- all_de %>%
    filter(subtype == subtype_name,
           condition1 == cond1,
           condition2 == cond2)
  
  # guard for empty results
  if (nrow(df) == 0) {
    message("No DE results found for: ", subtype_name, " ", cond1, " vs ", cond2)
    return(NULL)
  }
  
  # guard for fewer than 10 significant genes in either direction
  up_genes   <- df$avg_log2FC[df$avg_log2FC > logfc_threshold & 
                               df$p_val_adj_global < p_adj_threshold]
  down_genes <- df$avg_log2FC[df$avg_log2FC < -logfc_threshold & 
                               df$p_val_adj_global < p_adj_threshold]
  
  top_up_threshold   <- ifelse(length(up_genes) >= 10,
                               sort(up_genes, decreasing = TRUE)[10],
                               ifelse(length(up_genes) > 0, 
                                      min(up_genes, na.rm = TRUE), NA))
  
  top_down_threshold <- ifelse(length(down_genes) >= 10,
                               sort(down_genes, decreasing = FALSE)[10],
                               ifelse(length(down_genes) > 0,
                                      max(down_genes, na.rm = TRUE), NA))
  
  df <- df %>%
    mutate(
      direction = case_when(
        avg_log2FC > logfc_threshold  & p_val_adj_global < p_adj_threshold ~ "Up",
        avg_log2FC < -logfc_threshold & p_val_adj_global < p_adj_threshold ~ "Down",
        TRUE ~ "NS"
      ),
      label = case_when(
        direction == "Up"   & !is.na(top_up_threshold) & 
          avg_log2FC >= top_up_threshold   ~ gene,
        direction == "Down" & !is.na(top_down_threshold) & 
          avg_log2FC <= top_down_threshold ~ gene,
        TRUE ~ ""
      )
    )
  
  ggplot(df, aes(x = avg_log2FC,
                  y = -log10(p_val_adj_global),
                  color = direction,
                  label = label)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_text_repel(size = 3, max.overlaps = 20,
                    box.padding = 0.5, show.legend = FALSE) +
    scale_color_manual(values = c("Up"   = "#C00000",
                                   "Down" = "#2E75B6",
                                   "NS"   = "grey")) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dashed") +
    geom_hline(yintercept = -log10(p_adj_threshold), 
               linetype = "dashed") +
    theme_classic() +
    labs(title = paste(subtype_name, ":", cond1, "vs", cond2),
         x = "log2 Fold Change",
         y = "-log10 adjusted p-value",
         color = "Direction")
}


# --- For this project I'm only going to focus on Plasma & Memory B cells 
plasma_healthy_primary <- make_volcano("Plasma cells", "Primary", "Healthy")
plasma_healthy_mets <- make_volcano("Plasma cells", "Metastatic", "Healthy")
plasma_primary_mets <- make_volcano("Plasma cells", "Metastatic", "Primary")
memory_healthy_primary <- make_volcano("Memory B cells", "Primary", "Healthy")
memory_healthy_mets <- make_volcano("Memory B cells", "Metastatic", "Healthy")
memory_primary_mets <- make_volcano("Memory B cells", "Metastatic", "Primary")


# --- find a way to save these as well all at once? 
save_multiple_plots <- function(plots, filenames, folder = figure_directory) {
  dir.create(folder, showWarnings = FALSE)
  for (i in seq_along(plots)) {
    ggsave(
      filename = file.path(folder, filenames[i]),
      plot = plots[[i]]
    )
  }
}


save_multiple_plots(
  plots = list(plasma_healthy_primary,
               plasma_healthy_mets,
               plasma_primary_mets,
               memory_healthy_primary,
               memory_healthy_mets,
               memory_primary_mets),
  filenames = c("05_plasma_healthy_primary_deg.png", 
                "05_plasma_healthy_mets_deg.png", 
                "05_plasma_primary_mets_deg.png",
                "05_memory_healthy_primary_deg.png", 
                "05_memory_healthy_mets_deg.png", 
                "05_memory_primary_mets_deg.png")
)

# --- Save DE results as CSV
write.csv(all_de, file.path(output_directory,"05_bcell_DEG_all.csv"))
message("Saved: 05_bcell_DEG_all.csv")

# --- Save Signifcant DE results as CSV
write.csv(all_de_sig, file.path(output_directory,"05_bcell_DEG_significant.csv"))
message("Saved: 05_bcell_DE_significant.csv")

# --- Save Seurat Object
saveRDS(b_cells, file.path(output_directory,"05_bcell_gc_seurat_b_cell_deg_analysis.rds"))
message("Saved: 05_bcell_gc_seurat_b_cell_deg_analysis.rds")