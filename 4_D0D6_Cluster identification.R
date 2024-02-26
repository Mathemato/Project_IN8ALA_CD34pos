# Cluster identification in D0 dataset

# Automated cluster identification through SingleR package and data from Novershtern N et al. (2011)----
library(Seurat)
library(SingleR)
library(celldex)
library(scater)
annotations <- readRDS("C:/Data.path/annotations.rds")

# Loading of reference dataset from celldex package
ref.se <- NovershternHematopoieticData()

# Creation of a matrix with normalized counts using GetAssayData() function for each cluster.
# Example below for the cluster 0 of D0 dataset
Seurat.analyse.D0<-subset(seurat_integrated, subset = orig.ident=="D0")
saveRDS(Seurat.analyse.D0,"C:/Data.path/Seurat_analyse_D0.rds")
cluster0<-subset(Seurat.analyse.D0, subset=seurat_clusters=="0")
cl<-GetAssayData(object=cluster0[['RNA']], slot = "counts")

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
library(ggplot2)
library(dplyr)

# Heatmap with genes specifically associated with stemness or lineage priming in D0 dataset
DefaultAssay(seurat_integrated)<-"RNA"
Seurat.analyse.D0<-subset(seurat_integrated, subset = orig.ident=="D0")
Seurat.analyse.D0 <- ScaleData(Seurat.analyse.D0)

saveRDS(Seurat.analyse.D0,"C:/Data.path/seurat_integrated.rds")

Seurat.analyse.D0$seurat_clusters<-factor(Seurat.analyse.D0$seurat_clusters,levels = c("8","1","12","3","11","4","7","13","14","2","0","9","15","5","6","10","17","16"))
DoHeatmap(Seurat.analyse.D0,features=c("VIM","AVP","PTMS","ID1","CRHBP","APOC1","GATA1","TFRC","TPSAB1","HDC","MS4A3","ALOX5AP","LYZ","FCER1G","IGHM","LTB","HOPX","HOXA9","MME","TNFRSF18","LST1","CD79A","CD3E"),slot = "scale.data",label = TRUE , group.by = "seurat_clusters",group.bar = TRUE, disp.max = 4,group.colors=c("blue3","blueviolet","mediumorchid","mediumorchid2","yellow","gold","orange","red","goldenrod","green4","springgreen3","springgreen4","chartreuse3","green","chartreuse2","springgreen2","aquamarine","springgreen")) + scale_fill_gradientn(colors = c("Green", "Red"))+ theme(text = element_text(size = 30))

# DEG for each cluster relative to all remaining clusters
DefaultAssay(Seurat.analyse.D0)<-"RNA"
D0Markers <- FindAllMarkers(object = Seurat.analyse.D0, 
                            only.pos = FALSE,
                            logfc.threshold = 0.3,
                            min.pct = 0.5)

D0Markers<- inner_join(x = D0Markers, 
                       y = annotations[, c("gene_name", "description")],
                       by = c("gene" = "gene_name")) %>% unique()
D0Markers<- D0Markers %>%
  dplyr::arrange(cluster)
D0Markers$diffexpressed<- "NO"
D0Markers$diffexpressed[D0Markers$avg_log2FC > 0.5 & D0Markers$p_val_adj < 0.05] <- "UP"
D0Markers$diffexpressed[D0Markers$avg_log2FC < -0.5 & D0Markers$p_val_adj < 0.05] <- "DOWN"
D0Markers$delabel <- NA
D0Markers$delabel[D0Markers$diffexpressed != "NO"] <- D0Markers$gene[D0Markers$diffexpressed != "NO"]
D0Markers<- D0Markers[ , c(6, 7, 2:4, 1, 5,8:10)]
write.table(D0Markers, file = "C:/Data.path/D0ClusterMarkers.xls", sep="\t", quote = TRUE, dec = ",", row.names = FALSE, col.names = TRUE)          

# Visualization of D0 dataset
DimPlot(Seurat.analyse.D0,label = T, label.color = "black", pt.size=0.01, cols = c("springgreen3","blueviolet","green4","mediumorchid2","gold","green","chartreuse2","orange","blue3","springgreen4","springgreen2","yellow","mediumorchid4","red","goldenrod","chartreuse3","aquamarine","springgreen"))
DimPlot(Seurat.analyse.D0,label = F, label.color = "black", pt.size=0.01, cols = c("springgreen3","blueviolet","green4","mediumorchid2","gold","green","chartreuse2","orange","blue3","springgreen4","springgreen2","yellow","mediumorchid4","red","goldenrod","chartreuse3","aquamarine","springgreen"))


