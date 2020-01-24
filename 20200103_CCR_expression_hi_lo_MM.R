#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Collate all data into SingleCellExperiment object, and filter cells and genes based on available data. Mouse cells and genes are removed.
# --------------------------------------------------------------------------

# Set working directory, load data and libraries
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  location <- "/Users/mac/cloudstor/"
} else {
  location <- "/share/ScratchGeneral/scoyou/"
}

setwd(paste0(location, "projects/myeloma_osteocytes_cohort1/project_results/tumor_SC/"))

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('EnsDb.Hsapiens.v75')
library('grid')
library('ggplot2')
library('ggridges')
library('viridis')
library('Matrix')
library('reshape2')

# Read CSV with Chemokine Ligand-Receptor detail
CCRs <- read.csv("Myeloma male count file/data/20200124_Chemokine_receptor_targets_curated.csv", header = TRUE)

# Make single cell object using counts matrix
counts_matrix <- read.table("Myeloma male count file/RSEM_MM.gene.counts.matrix")
raw_experiment <- SingleCellExperiment(assays = list(counts = counts_matrix))

# Extract cell info from cell ids
cell_ids <- colnames(raw_experiment)
cell_ids <- gsub(".*eGFP.","",cell_ids)
cell_info <- data.frame(cell_ids) %>%
  separate(cell_ids, c("Cells", "DiD", NA), "_", remove = FALSE, extra = "merge")
raw_experiment$DiD <- cell_info$DiD
raw_experiment$Cells <- cell_info$Cells

# Subset to single cells and chemokine receptor genes
raw_experiment <-raw_experiment[,raw_experiment$Cells == "1C"]
CCR.sce <- raw_experiment[unique(as.character(CCRs$Mouse_receptor_Ensembl)), ]

# Calculate percentage of CCR positive cells in Proliferative and Dormant tumor cells
CCR.percent.Neg <- rowSums(counts(CCR.sce[,CCR.sce$DiD == "Neg"]) > 0)/ncol(CCR.sce[,CCR.sce$DiD == "Neg"])*100 # Calculates percent of DiD negative cells express each receptor
CCR.percent.Hi <- rowSums(counts(CCR.sce[,CCR.sce$DiD == "Hi"]) > 0)/ncol(CCR.sce[,CCR.sce$DiD == "Neg"])*100 # Calculates percent of DiD Hi cells express each receptor
CCR.percent <- data.frame("MM_Proliferative" = CCR.percent.Neg, "MM_Dormant" = CCR.percent.Hi)

idx <- match(rownames(CCR.percent), CCRs$Mouse_receptor_Ensembl)
CCR.percent$Gene <- CCRs$Mouse_receptor_GeneSymbol [idx]
CCR.percent <- CCR.percent[order(CCR.percent$Gene),]

# Make dotplot relecting percentage of expressing cells per CCR
CCR.percent.melt <- melt(CCR.percent)
colnames(CCR.percent.melt) <- c("Gene", "Cell_type", "Percent")
CCR.percent.melt$Gene <- factor(CCR.percent.melt$Gene, levels = rev(CCR.percent$Gene))
# CCR.percent.melt$Percent[CCR.percent.melt$Percent == 0] <- NA # Removes 0 values from plot if desired

ggplot(CCR.percent.melt, aes(x = Gene, y = Cell_type)) +
  geom_point(aes(size = Percent, color = Percent)) +
  scale_color_viridis_c(option = "B", limits = c(0, 100)) +
  scale_radius(range = c(0, 6), limits = c(0, 100)) +
  coord_flip() +
  theme_minimal() +
  ggsave("MM_CCR_expression.pdf", useDingbats = FALSE)
  


