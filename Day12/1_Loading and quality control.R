library(Seurat)
library(SeuratObject)
library(dplyr)
library(dbplyr)
library(Matrix)
library(patchwork)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(umap)
library(rlang)
library(magrittr)
library(tidyverse)
library(devtools)
library(BiocManager)
library(cli)
library(SingleCellExperiment)
library(AnnotationHub)
library(ensembldb)
library(scales)
library(openxlsx)
library(ggpubr)
library(writexl)
library(limma)
library(sctransform)
library(multtest)
library(metap)

#Load count files :----
counts.CTRLJ121 = Read10X(data.dir = "C:/Users/13700976/Desktop/Analyses single cells Mathilde/DATAS BRUTES/filtered_feature_bc_matrixJ12Cont1/", strip.suffix = TRUE)

counts.CTRLJ122 = Read10X(data.dir = "C:/Users/13700976/Desktop/Analyses single cells Mathilde/DATAS BRUTES/filtered_feature_bc_matrixJ12Cont2/", strip.suffix = TRUE)

counts.JAJ12 = Read10X(data.dir = "C:/Users/13700976/Desktop/Analyses single cells Mathilde/DATAS BRUTES/filtered_feature_bc_matrixJ12JA/", strip.suffix = TRUE)

#Create Seurat object :----
so.CTRLJ121 <- CreateSeuratObject(counts = counts.CTRLJ121, project = "CTRLJ121",min.cells = 3,min.features = 200)
table(so.CTRLJ121$orig.ident) 

so.CTRLJ122 <- CreateSeuratObject(counts = counts.CTRLJ122, project = "CTRLJ122",min.cells = 3,min.features = 200)
table(so.CTRLJ122$orig.ident) 

so.JAJ12 <- CreateSeuratObject(counts = counts.JAJ12, project = "J+AJ12")
table(so.JAJ12$orig.ident) 

CTRLJ12<- merge(so.CTRLJ121, y =so.CTRLJ122, add.cell.ids = c("ctrl1", "ctrl2"), project = "RNA")
head(colnames(CTRLJ12))
JAJ12<- so.JAJ12

#Create new column "condition" 
CTRLJ12$condition <- "CTRL"
JAJ12$condition <- "JNK+ALA"
MergeJ12<- merge(CTRLJ12, y =JAJ12, add.cell.ids = c("CTRL","JNK+ALA"), project = "RNA")

#Add gene number by UMI (Novelty score)
MergeJ12$log10GenesPerUMI <- log10(MergeJ12$nFeature_RNA) / log10(MergeJ12$nCount_RNA)

#Determine mito genes
MergeJ12$mitoRatio <- PercentageFeatureSet(object = MergeJ12, pattern = "^MT-")
MergeJ12$mitoRatio <- MergeJ12@meta.data$mitoRatio / 100

# Visualize the number UMIs/transcripts per cell
MergeJ12@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
MergeJ12@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
MergeJ12@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Visualize the distribution of mitochondrial gene expression detected per cell
MergeJ12@meta.data %>% 
  ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
MergeJ12@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 400) +
  facet_wrap(~orig.ident)

# Filter out low quality cells using selected thresholds - these will change with experiment
Filtered_MergeJ12 <- subset(x = MergeJ12, 
                            subset= (nCount_RNA >= 500) & 
                              (nFeature_RNA >= 400) & 
                              (log10GenesPerUMI > 0.80) & 
                              (mitoRatio < 0.20))

Filtered_MergeJ12@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 400) +
  facet_wrap(~orig.ident)

# Extract counts
counts <- GetAssayData(object = Filtered_MergeJ12, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 5 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 5

# Only keeping those genes expressed in more than 5 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
Filtered_MergeJ12 <- CreateSeuratObject(filtered_counts, meta.data = Filtered_MergeJ12@meta.data)

# ReVisualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
Filtered_MergeJ12@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident)