#Package of interest
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(sctransform)
library(tidyverse)

#Loading of cell ranger files
counts.D0 = Read10X(data.dir = "C:/Data.path/filtered_feature_bc_matrixJ0/", strip.suffix = TRUE)
counts.CTRL6h = Read10X(data.dir= "C:/Data.path/filtered_feature_bc_matrix6hCont/", strip.suffix = TRUE)
counts.JNK6h = Read10X(data.dir = "C:/Data.path/filtered_feature_bc_matrix6hJNK/", strip.suffix = TRUE)
counts.ALA6h = Read10X(data.dir = "C:/Data.path/filtered_feature_bc_matrix6hALA/", strip.suffix = TRUE)
counts.IN8ALA6h = Read10X(data.dir = "C:/Data.path/filtered_feature_bc_matrix6hJA/", strip.suffix = TRUE)

#Creating SeuratObjects
so.D0 <- CreateSeuratObject(counts = counts.D0, project = "D0")
so.CTRL6h <- CreateSeuratObject(counts = counts.CTRL6h, project = "CTRL6h")
so.JNK6h <- CreateSeuratObject(counts = counts.JNK6h, project = "JNK6h")
so.ALA6h <- CreateSeuratObject(counts = counts.ALA6h, project = "ALA6h")
so.IN8ALA6h <- CreateSeuratObject(counts = counts.IN8ALA6h, project = "IN8ALA6h")


#Reducing the dataset to be similar to others
so.ALA6h<- so.ALA6h[, sample(colnames(so.ALA6h), size =6000, replace=F)]

#Merging all SeuratObjects
Mergeso<- merge(so.D0, y =c(so.CTRL6h,so.JNK6h,so.ALA6h,so.IN8ALA6h), add.cell.ids = c("D0","CTRL", "JNK","ALA","IN8ALA"), project = "molecules induction")

saveRDS(Mergeso,"C:/Data.path/Mergeso.rds")

#Check the Quality Control
Mergeso$log10GenesPerUMI <- log10(Mergeso$nFeature_RNA) / log10(Mergeso$nCount_RNA)
Mergeso$mitoRatio <- PercentageFeatureSet(object = Mergeso, pattern = "^MT-")
Mergeso$mitoRatio <- Mergeso@meta.data$mitoRatio / 100
metadatamerge<- Mergeso@meta.data
metadatamerge$cells<- rownames(metadatamerge)
metadatamerge$sample <- NA
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^D0_"))] <- "D0"
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^CTRL_"))] <- "CTRL6h"
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^JNK_"))] <- "JNK6h"
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^ALA_"))] <- "ALA6h"
metadatamerge$sample[which(str_detect(metadatamerge$cells, "^IN8ALA_"))] <- "IN8ALA6h"
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

saveRDS(Filtered_Mergeso,"C:/Data.path/Filtered_Mergeso.rds")
