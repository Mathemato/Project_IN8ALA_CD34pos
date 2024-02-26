# Cluster identification in D12 dataset

# Automated cluster identification through SingleR package and data from Novershtern N et al. (2011)----
library(Seurat)
library(SingleR)
library(celldex)
library(scater)
annotations <- readRDS("C:/Data.path/annotations.rds")

# Loading of reference dataset from celldex package
ref.se <- NovershternHematopoieticData()

# Creation of a matrix with normalized counts using GetAssayData() function for each cluster.
# Example below for the cluster 0 of D12 dataset
Cl0<-subset(MergeD12_norm, subset=seurat_clusters=="0")
cl<-GetAssayData(object=Cl0[['RNA']], slot = "counts")

# SingleR cluster annotation using reference dataset labelling  
pred.cl <- SingleR(test = cl, ref = ref.se, assay.type.test=1,
                   labels = ref.se$label.fine)

# Visualization of cluster identification 
table(pred.cl$labels)
plotScoreHeatmap(pred.cl)
plotDeltaDistribution(pred.cl, ncol=3)

# Quality control of identification 
summary(is.na(pred.cl$pruned.labels))


# Manual checking of cluster identity----
library(Seurat)
library(ggplot2)
library(dplyr)

# Heatmap with genes specifically associated with stemness or lineage committment
DefaultAssay(MergeD12_norm)<-"RNA"
MergeD12_norm$seurat_clusters<-factor(MergeD12_norm$seurat_clusters,levels = c("1","7","9","4","2","5","3","6","8","12","10","11","0","13"))
DoHeatmap(MergeD12_norm,features=c("HOXA9","SELL","SOX4","MPO","PRTN3","AZU1","ELANE","MS4A7","CD14","HLA-DPA1","HLA-DPB1","HLA-DQB1","GATA2","RHEX","CPA3","MS4A2","TPSAB1","APOC1","TFRC","HBD","GATA1","PF4","PPBP","GP9","GP1BB","PLEK"),label = TRUE , group.by = "seurat_clusters",group.bar = TRUE, disp.max = 4,group.colors=c("blue","mediumorchid4","purple2","orchid","mediumorchid1","gold","yellow","chartreuse1","chartreuse3","goldenrod1","goldenrod3","red2","tomato1","tomato4")) + scale_fill_gradientn(colors = c("Green", "Red"))+ theme(text = element_text(size = 30))

DimPlot(MergeD12_norm,label = T, label.color = "black", pt.size=0.01, cols = c("tomato1","blue","mediumorchid1","yellow","orchid","gold","chartreuse1","mediumorchid4","chartreuse3","purple2","goldenrod3","red2","goldenrod1","tomato4"))
DimPlot(MergeD12_norm,label = F, label.color = "black", pt.size=0.01, cols = c("tomato1","blue","mediumorchid1","yellow","orchid","gold","chartreuse1","mediumorchid4","chartreuse3","purple2","goldenrod3","red2","goldenrod1","tomato4"))


# DEG for each cluster relative to all remaining clusters
DefaultAssay(MergeD12_norm)<-"RNA"
NormalizeData(MergeD12_norm)
D12Markers <- FindAllMarkers(object = MergeD12_norm, 
                             only.pos = FALSE,
                             logfc.threshold = 0.3,
                             min.pct = 0.5)
D12Markers<- inner_join(x = D12Markers, 
                        y = annotations[, c("gene_name", "description")],
                        by = c("gene" = "gene_name")) %>% unique()

D12Markers<- D12Markers[ , c(6, 7, 2:4, 1, 5,8)]

D12Markers<- D12Markers %>%
  dplyr::arrange(cluster)

#Add NA column
D12Markers$diffexpressed<- "NO"

#Identify overexpressed genes as "up"
D12Markers$diffexpressed[D12Markers$avg_log2FC > 0.5 & D12Markers$p_val_adj < 0.05] <- "UP"

#Identify downregulated genes as "down"
D12Markers$diffexpressed[D12Markers$avg_log2FC < -0.5 & D12Markers$p_val_adj < 0.05] <- "DOWN"
D12Markers$delabel <- NA
D12Markers$delabel[D12Markers$diffexpressed != "NO"] <- D12Markers$gene[D12Markers$diffexpressed != "NO"]

write.table(D12Markers, file = "C:/Data.path/D12Markers.xls", sep="\t", quote = TRUE, dec = ",", row.names = FALSE, col.names = TRUE)          





