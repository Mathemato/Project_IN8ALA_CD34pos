#Package of interest
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(sctransform)
library(tidyverse)
library(cowplot)

#Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)


#Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

#Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)

#Determine percent of variation associated with each PC
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100

#Calculate cumulative percents for each PC
cumu <- cumsum(pct)

#Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

#Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

#last point where change of % of variation is more than 0.1%.
co2

#Minimum of the two calculation
pcs <- min(co1, co2)

pcs

#Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:13)
#Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated, resolution = 0.6)

seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:13,
                             reduction = "pca")
#Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        repel = TRUE)+theme(text = element_text(face="bold"))


DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        repel = TRUE,
        split.by = "orig.ident")

#Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

#View table and save result
View(n_cells)

write.table(n_cells, file = "C:/Data.path/n_cells.xls", sep="\t", quote = TRUE, dec = ",", row.names = FALSE, col.names = TRUE)          

#Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

#Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

#Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:20),
             "ident",
             "UMAP_1", "UMAP_2")

#Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)

#Extract the UMAP coordinates for the first 10 cells
seurat_integrated@reductions$umap@cell.embeddings[1:10, 1:2]

#Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

#Plotting a UMAP plot for each of the PCs
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


print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

saveRDS(seurat_integrated,"C:/data.path/seurat_integrated.rds")
