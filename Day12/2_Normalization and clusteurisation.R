# Normalize the counts
seurat_phase <- NormalizeData(Filtered_MergeJ12)

# Load cell cycle markers
load("/Users/13700976/Desktop/Analyses single cells Mathilde/cycle.rda")
# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)  

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 3000, verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "condition",
        split.by = "phase")

# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))
# Plot the PCA 
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")


#Loop pour faire le SCTransform sur tous les échantillons 
MergeJ12_norm<-SCTransform(seurat_phase,vars.to.regress = c("S.Score","G2M.Score"))

#If needed: vars.to.regress = c("S.Score")

# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

DefaultAssay(MergeJ12_norm)<-"SCT"

# Run PCA
MergeJ12_norm <- RunPCA(object = MergeJ12_norm)

# Plot PCA
PCAPlot(MergeJ12_norm,
        split.by = "orig.ident") 

# Run UMAP
MergeJ12_norm <- RunUMAP(MergeJ12_norm, 
                         dims = 1:40,
                         reduction = "pca")
# Plot UMAP                             
DimPlot(MergeJ12_norm)  
MergeJ12_norm$orig.ident<-factor(MergeJ12_norm$orig.ident, levels = c("CTRLJ121","CTRLJ122","J+AJ12"))

DimPlot(MergeJ12_norm,
        split.by = "orig.ident")  
# Explore heatmap of PCs
DimHeatmap(MergeJ12_norm, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
# Printing out the most variable genes driving PCs
print(x = MergeJ12_norm[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = MergeJ12_norm, 
          ndims = 40)

# Determine percent of variation associated with each PC
pct <- MergeJ12_norm[["pca"]]@stdev / sum(MergeJ12_norm[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

# Determine the K-nearest neighbor graph
MergeJ12_norm <- FindNeighbors(object = MergeJ12_norm, 
                               dims = 1:16)
# Determine the clusters for various resolutions                                
MergeJ12_norm <- FindClusters(object = MergeJ12_norm, resolution = 0.6)

MergeJ12_norm <- RunUMAP(MergeJ12_norm, 
                         dims = 1:16,
                         reduction = "pca")
# Plot the UMAP
DimPlot(MergeJ12_norm,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        repel = TRUE)


DimPlot(MergeJ12_norm,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        repel = TRUE,
        split.by = "condition")

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(MergeJ12_norm, 
                     vars = c("ident", "condition")) %>%
  dplyr::count(ident, condition) %>%
  tidyr::spread(ident, n)
# View table
View(n_cells)

write.table(n_cells, file = "C:/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/CTRLJ12 + JA J12/ncells.xls", sep="\t", quote = TRUE, dec = ",", row.names = FALSE, col.names = TRUE)          


# Explore whether clusters segregate by cell cycle phase
DimPlot(MergeJ12_norm,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in MergeJ12_norm@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(MergeJ12_norm, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:20),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(MergeJ12_norm, 
                     vars = columns)

# Extract the UMAP coordinates for the first 10 cells
MergeJ12_norm@reductions$umap@cell.embeddings[1:10, 1:2]

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(MergeJ12_norm, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:13), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)


print(MergeJ12_norm[["pca"]], dims = 1:5, nfeatures = 5)

