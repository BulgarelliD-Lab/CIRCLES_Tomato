# Required libraries
library(qiime2R)
library(tidyverse)
library(vegan)
library(pheatmap)
library(readr)
library(dplyr)

################################################################################
############################## HEATMAP: AMR GENES ##############################
################################################################################

# Set working directory (relative to project root if needed)
setwd("Analysis/R_scripts/NMDS")

# Load FPKM abundance table (AMR genes)
abundance_table <- read_tsv("../../abundance_tables/ResFinder_genes_FPKM_genomic_20220524.txt")

# Load metadata
metadata <- read_tsv("../../Metadata/CN04_metadata.txt")

# Set sampleID as rownames
abundance_table <- abundance_table %>% column_to_rownames(var = "sampleID")

# Transpose to have features (genes) as rows
fpkm_ab <- abundance_table %>% t()

# Compute Bray-Curtis distance matrix between samples
distmatrix <- vegdist(t(fpkm_ab), method = "bray", na.rm = TRUE)

# Subset to top 25 most abundant features
fpkm_ab <- if (nrow(fpkm_ab) > 25) fpkm_ab[order(rowSums(fpkm_ab), decreasing = TRUE)[1:25], ] else fpkm_ab

# Log-transform data
fpkm_ab <- log(fpkm_ab + 1)

# Prepare annotation metadata
metadata_annotation <- metadata %>%
  select(sampleID, Field, Microenvironment, Genotype) %>%
  column_to_rownames("sampleID") %>%
  as.data.frame()

# Define color palettes for annotation
field_palette <- setNames(c("#ADD8E6", "#FFC0CB"), c("Larino-3", "Larino-4"))
microenv_palette <- setNames(c("#90EE90", "#FA8072"), c("Rhizosphere", "Bulk"))
genotype_palette <- setNames(c("#FFD700", "#87CEEB"), c("Abundo", "SV5197TP"))

# Combine color palettes
ann_colors <- list(Field = field_palette, 
                   Microenvironment = microenv_palette, 
                   Genotype = genotype_palette)

# Heatmap color palette
general <- colorRampPalette(c("black", "#19375D", "#00214B", "#EE1652", "#F68025", "#FFFFB2"))

# Output directory and filename
outdir <- "../../Plots/WP03_tomatoes_obs_phase/Resistome/Heatmaps/"
prefix <- "heatmap_AMR_genes_FPKM_"

# Generate heatmap and save to PDF
pdf(paste0(outdir, prefix, "WP03_tomatoes.pdf"), width = 13, height = 7)
plt1 <- pheatmap(fpkm_ab, 
                 margins = c(3, 3),
                 fontsize_row = 8, 
                 angle_col = 45, 
                 treeheight_row = 100, 
                 treeheight_col = 100,
                 scale = "none", 
                 col = general(100), 
                 clustering_method = "complete", 
                 clustering_distance_cols = distmatrix, 
                 clustering_distance_rows = "correlation", 
                 annotation = metadata_annotation,
                 annotation_colors = ann_colors,
                 show_colnames = TRUE)
dev.off()

################################################################################
########################### HEATMAP: AMR CLASSES ###############################
################################################################################

# Load FPKM abundance table (AMR classes)
abundance_table <- read_tsv("../../abundance_tables/ResFinder_AMR_class_FPKM_genomic_20220524.txt")

# Load metadata
metadata <- read_tsv("../../Metadata/CN04_metadata.txt")

# Set sampleID as rownames
abundance_table <- abundance_table %>% column_to_rownames(var = "sampleID")

# Transpose to have features (classes) as rows
fpkm_ab <- abundance_table %>% t()

# Compute Bray-Curtis distance matrix between samples
distmatrix <- vegdist(t(fpkm_ab), method = "bray", na.rm = TRUE)

# Subset to top 40 most abundant features
fpkm_ab <- if (nrow(fpkm_ab) > 40) fpkm_ab[order(rowSums(fpkm_ab), decreasing = TRUE)[1:40], ] else fpkm_ab

# Log-transform data
fpkm_ab <- log(fpkm_ab + 1)

# Prepare annotation metadata
metadata_annotation <- metadata %>%
  select(sampleID, Field, Microenvironment, Genotype) %>%
  column_to_rownames("sampleID") %>%
  as.data.frame()

# Define color palettes for annotation
field_palette <- setNames(c("#ADD8E6", "#FFC0CB"), c("Larino-3", "Larino-4"))
microenv_palette <- setNames(c("#90EE90", "#FA8072"), c("Rhizosphere", "Bulk"))
genotype_palette <- setNames(c("#FFD700", "#87CEEB"), c("Abundo", "SV5197TP"))

# Combine color palettes
ann_colors <- list(Field = field_palette, 
                   Microenvironment = microenv_palette, 
                   Genotype = genotype_palette)

# Heatmap color palette
general <- colorRampPalette(c("black", "#19375D", "#00214B", "#EE1652", "#F68025", "#FFFFB2"))

# Output directory and filename
outdir <- "../../Plots/WP03_tomatoes_obs_phase/Resistome/Heatmaps/"
prefix <- "heatmap_AMR_class_FPKM_"

# Generate heatmap and save to PDF
pdf(paste0(outdir, prefix, "WP03_tomatoes.pdf"), width = 13, height = 7)
plt1 <- pheatmap(fpkm_ab, 
                 margins = c(3, 3),
                 fontsize_row = 8, 
                 angle_col = 45, 
                 treeheight_row = 100, 
                 treeheight_col = 100,
                 scale = "none", 
                 col = general(100), 
                 clustering_method = "complete", 
                 clustering_distance_cols = distmatrix, 
                 clustering_distance_rows = "correlation", 
                 annotation = metadata_annotation,
                 annotation_colors = ann_colors,
                 show_colnames = TRUE)
dev.off()
