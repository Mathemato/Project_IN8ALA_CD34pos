#Package of interest
library(Seurat)
library(SeuratObject)
library(dplyr)
library(dbplyr)
library(Matrix)
library(patchwork)
library(ggrepel)
library(ggplot2)
library(umap)
library(rlang)
library(magrittr)
library(tidyverse)
library(devtools)
library(BiocManager)
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

#Package to install  
BiocManager::install("multtest") 
#Then install Seurat again, Then install the 3 others : 
BiocManager::install("SingleCellExperiment")
BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")

#Loading of cell ranger files
counts.J0 = Read10X(data.dir = "/Users/13700976/Desktop/Analyses single cells Mathilde/DATAS BRUTES/filtered_feature_bc_matrixJ0/", strip.suffix = TRUE)
counts.CTRL6h = Read10X(data.dir= "/Users/13700976/Desktop/Analyses single cells Mathilde/DATAS BRUTES/filtered_feature_bc_matrix6hCont/", strip.suffix = TRUE)
counts.JNK6h = Read10X(data.dir = "/Users/13700976/Desktop/Analyses single cells Mathilde/DATAS BRUTES/filtered_feature_bc_matrix6hJNK/", strip.suffix = TRUE)
counts.ALA6h = Read10X(data.dir = "/Users/13700976/Desktop/Analyses single cells Mathilde/DATAS BRUTES/filtered_feature_bc_matrix6hALA/", strip.suffix = TRUE)
counts.JA6h = Read10X(data.dir = "/Users/13700976/Desktop/Analyses single cells Mathilde/DATAS BRUTES/filtered_feature_bc_matrix6hJA/", strip.suffix = TRUE)

#Creating SeuratObjects
so.J0 <- CreateSeuratObject(counts = counts.J0, project = "J0")
so.CTRL6h <- CreateSeuratObject(counts = counts.CTRL6h, project = "CTRL6h")
so.ALA6h <- CreateSeuratObject(counts = counts.ALA6h, project = "ALA6h")
so.JA6h <- CreateSeuratObject(counts = counts.JA6h, project = "J+A6h")
so.JNK6h <- CreateSeuratObject(counts = counts.JNK6h, project = "JNK6h")

#Reducing the dataset to be similar to others
so.ALA6h<- so.ALA6h[, sample(colnames(so.ALA6h), size =6000, replace=F)]

#Merging all SeuratObjects
Mergeso<- merge(so.J0, y =c(so.CTRL6h,so.JNK6h,so.ALA6h,so.JA6h), add.cell.ids = c("J0","ctrl", "JNK","ALA","JNK+ALA"), project = "molecules induction")
table(Mergeso$orig.ident)
view(Mergeso@meta.data)

#Check the Quality Control
Mergeso$log10GenesPerUMI <- log10(Mergeso$nFeature_RNA) / log10(Mergeso$nCount_RNA)
Mergeso$mitoRatio <- PercentageFeatureSet(object = Mergeso, pattern = "^MT-")
Mergeso$mitoRatio <- Mergeso@meta.data$mitoRatio / 100
metadatamerge<- Mergeso@meta.data
metadatamergeso$cells<- rownames(metadatamergeso)
metadatamerge$sample <- NA
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^J0_"))] <- "J0"
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^ctrl_"))] <- "CTRL6h"
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^JNK_"))] <- "JNK6h"
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^ALA_"))] <- "ALA6h"
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^JNK+ALA_"))] <- "JNK+ALA6h"
metadatamerge <- metadatamerge %>%
   dplyr::rename(orig.ident = orig.odent,
   nCount_RNA = nUMI,
   nFeature_RNA = nGene)
Mergeso@meta.data <- metadatamerge
Mergeso@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
Mergeso@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
Mergeso@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
Mergeso@meta.data %>% 
  ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
Mergeso@meta.data %>% 
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

#Filter the cells from the bounds determined above
Filtered_Mergeso <- subset(x = Mergeso, subset= (nCount_RNA >= 500) & 
                             (nFeature_RNA >= 400) & 
                             (log10GenesPerUMI > 0.80) & 
                             (mitoRatio < 0.20))
#Extract counts
counts <- GetAssayData(object = Filtered_Mergeso, slot = "counts")

#Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

#Sums all TRUE values and returns TRUE if more than 5 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 5

#Only keeping those genes expressed in more than 5 cells
filtered_counts <- counts[keep_genes, ]

#Reassign to filtered Seurat object
Filtered_Mergeso <- CreateSeuratObject(filtered_counts, meta.data = Filtered_Mergeso@meta.data)

#ReVisualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
Filtered_Mergeso@meta.data %>% 
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