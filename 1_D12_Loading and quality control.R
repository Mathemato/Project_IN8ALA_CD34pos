library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(sctransform)
library(tidyverse)

#Load count files :----
counts.CTRLD121 = Read10X(data.dir = "C:/Data.path/filtered_feature_bc_matrixJ12Cont1/", strip.suffix = TRUE)
counts.CTRLD122 = Read10X(data.dir = "C:/Data.path/filtered_feature_bc_matrixJ12Cont2/", strip.suffix = TRUE)
counts.IN8ALAD12 = Read10X(data.dir = "C:/Data.path/filtered_feature_bc_matrixJ12JA/", strip.suffix = TRUE)

#Create Seurat object :----
so.CTRLD121 <- CreateSeuratObject(counts = counts.CTRLD121, project = "CTRLD121",min.cells = 3,min.features = 200)
table(so.CTRLD121$orig.ident) 

so.CTRLD122 <- CreateSeuratObject(counts = counts.CTRLD122, project = "CTRLD122",min.cells = 3,min.features = 200)
table(so.CTRLD122$orig.ident) 

so.IN8ALAD12 <- CreateSeuratObject(counts = counts.IN8ALAD12, project = "IN8ALAD12")
table(so.IN8ALAD12$orig.ident) 

CTRLD12<- merge(so.CTRLD121, y =so.CTRLD122, add.cell.ids = c("CTRL1", "CTRL2"), project = "RNA")
head(colnames(CTRLD12))
IN8ALAD12<- so.IN8ALAD12

#Create new column : "condition" 
CTRLD12$condition <- "CTRL"
IN8ALAD12$condition <- "IN8ALA"
MergeD12<- merge(CTRLD12, y =IN8ALAD12, add.cell.ids = c("CTRL","IN8ALA"), project = "RNA")

#Add gene number by UMI (Novelty score)
MergeD12$log10GenesPerUMI <- log10(MergeD12$nFeature_RNA) / log10(MergeD12$nCount_RNA)

#Determine mito genes
MergeD12$mitoRatio <- PercentageFeatureSet(object = MergeD12, pattern = "^MT-")
MergeD12$mitoRatio <- MergeD12@meta.data$mitoRatio / 100

# Visualize the number UMIs/transcripts per cell
MergeD12@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
MergeD12@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
MergeD12@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Visualize the distribution of mitochondrial gene expression detected per cell
MergeD12@meta.data %>% 
  ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
MergeD12@meta.data %>% 
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
Filtered_MergeD12 <- subset(x = MergeD12, 
                            subset= (nCount_RNA >= 500) & 
                              (nFeature_RNA >= 400) & 
                              (log10GenesPerUMI > 0.80) & 
                              (mitoRatio < 0.20))

Filtered_MergeD12@meta.data %>% 
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
counts <- GetAssayData(object = Filtered_MergeD12, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 5 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 5

# Only keeping those genes expressed in more than 5 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
Filtered_MergeD12 <- CreateSeuratObject(filtered_counts, meta.data = Filtered_MergeD12@meta.data)

# ReVisualize the correlation between genes detected and number of UMIs and determine whether there is a strong presence of cells with low numbers of genes/UMIs
Filtered_MergeD12@meta.data %>% 
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

saveRDS(Filtered_MergeD12,"C:/Data.path/Filtered_MergeD12.rds")
