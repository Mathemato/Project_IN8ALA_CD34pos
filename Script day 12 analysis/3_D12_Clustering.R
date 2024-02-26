# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Explore heatmap of PCs
DimHeatmap(MergeD12_norm, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
# Printing out the most variable genes driving PCs
print(x = MergeD12_norm[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = MergeD12_norm, 
          ndims = 40)

# Determine percent of variation associated with each PC
pct <- MergeD12_norm[["pca"]]@stdev / sum(MergeD12_norm[["pca"]]@stdev) * 100

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
MergeD12_norm <- FindNeighbors(object = MergeD12_norm, 
                               dims = 1:16)
# Determine the clusters for various resolutions                                
MergeD12_norm <- FindClusters(object = MergeD12_norm, resolution = 0.6)

MergeD12_norm <- RunUMAP(MergeD12_norm, 
                         dims = 1:16,
                         reduction = "pca")
# Plot the UMAP
DimPlot(MergeD12_norm,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        repel = TRUE)


DimPlot(MergeD12_norm,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        repel = TRUE,
        split.by = "condition")

# Dimplot color modification
DimPlot(MergeD12_norm,label = T, label.color = "black", pt.size=0.01, cols = c("tomato1","blue","mediumorchid1","yellow","orchid","gold","chartreuse1","mediumorchid4","chartreuse3","purple2","goldenrod3","red2","goldenrod1","tomato4"))
DimPlot(MergeD12_norm,label = F, label.color = "black", pt.size=0.01, cols = c("tomato1","blue","mediumorchid1","yellow","orchid","gold","chartreuse1","mediumorchid4","chartreuse3","purple2","goldenrod3","red2","goldenrod1","tomato4"))

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
D12n_cells <- FetchData(MergeD12_norm, 
                     vars = c("ident", "condition")) %>%
  dplyr::count(ident, condition) %>%
  tidyr::spread(ident, n)
# View table
View(D12n_cells)

write.table(D12n_cells, file = "C:/Data.path/D12ncells.xls", sep="\t", quote = TRUE, dec = ",", row.names = FALSE, col.names = TRUE)          

# Explore whether clusters segregate by cell cycle phase
DimPlot(MergeD12_norm,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in MergeD12_norm@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(MergeD12_norm, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Define the information in the seurat object of interest
columns <- c(paste0("PC_", 1:20),
             "ident",
             "UMAP_1", "UMAP_2")

# Extract this data from the seurat object
pc_data <- FetchData(MergeD12_norm, 
                     vars = columns)

# Extract the UMAP coordinates for the first 10 cells
MergeD12_norm@reductions$umap@cell.embeddings[1:10, 1:2]

# Add cluster label to center of cluster on UMAP
umap_label <- FetchData(MergeD12_norm, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plot a UMAP plot for each of the PCs
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


print(MergeD12_norm[["pca"]], dims = 1:5, nfeatures = 5)




