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

#View table
View(n_cells)

write.table(n_cells, file = "C:/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/analyse 6h + J0 intégrée/ncells2.xls", sep="\t", quote = TRUE, dec = ",", row.names = FALSE, col.names = TRUE)          


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

#identification of cell types on J0 subset 
Seurat.analyse.J0<-subset(seurat_integrated,subset = orig.ident == 'J0')

#Find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(Seurat.analyse.J0)<-"RNA"
NormalizeData(Seurat.analyse.J0)
MarqueursJ0 <- FindAllMarkers(object = Seurat.analyse.J0, 
                              only.pos = FALSE,
                              logfc.threshold = 0.3,
                              min.pct = 0.5)

MarqueursJ0<- inner_join(x = Seurat.analyse.J0, 
                         y = annotations[, c("gene_name", "description")],
                         by = c("gene" = "gene_name")) %>% unique()
MarqueursJ0<- MarqueursJ0[ , c(6, 7, 2:4, 1, 5,8)]
MarqueursJ0<- MarqueursJ0 %>%
  dplyr::arrange(cluster)
MarqueursJ0$diffexpressed<- "NO"
MarqueursJ0$diffexpressed[MarqueursJ0$avg_log2FC > 0.5 & MarqueursJ0$p_val_adj < 0.05] <- "UP"
MarqueursJ0$diffexpressed[MarqueursJ0$avg_log2FC < -0.5 & MarqueursJ0$p_val_adj < 0.05] <- "DOWN"
MarqueursJ0$delabel <- NA
MarqueursJ0$delabel[MarqueursJ0$diffexpressed != "NO"] <- MarqueursJ0$gene[MarqueursJ0$diffexpressed != "NO"]

#Feature plot thanks to the Nicole Mende (Blood 2022) and Lars Velten (Nature Cell Biology 2017) papers. 
#Module 5 HSC to Neuro prog
FeaturePlot(Seurat.analyse.J0, features = c("CEBPA","CEBPD"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 5 HSC to Neuro prog
FeaturePlot(Seurat.analyse.J0, features = c("CEBPA","CEBPD"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 10 HSC to Pre-B cell prog
FeaturePlot(Seurat.analyse.J0, features = c("EBF1","ID3"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 13 HSC to Pre-B cell/Mono DC  prog
FeaturePlot(Seurat.analyse.J0, features = c("FLT3","SATB1"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 16 HSC to Erythro prog
FeaturePlot(Seurat.analyse.J0, features = c("GATA1","KLF1"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 17 HSC to MK prog
FeaturePlot(Seurat.analyse.J0, features = c("GATA2","NFE2"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 19 HSC to MK prog
FeaturePlot(Seurat.analyse.J0, features = c("GP1BB","PBX1"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 20 HSC to Granulo prog
FeaturePlot(Seurat.analyse.J0, features = c("HDC","LMO4"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 21 HSC primée
FeaturePlot(Seurat.analyse.J0, features = c("HLF","ZFP36L2"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 22 HSC Slow cycling 
FeaturePlot(Seurat.analyse.J0, features = c("HOXA3","HOXB6"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 36 HSC to Neutro 
FeaturePlot(Seurat.analyse.J0, features = c("SPI1","GFI1"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Module 36 HSC to Neutro 
FeaturePlot(Seurat.analyse.J0, features = c("TAL1","HSF1"), keep.scale = c("all"), pt.size = 0.5, slot = "data",label = TRUE,reduction = "umap", blend = TRUE)

#Primitif clusters
FeaturePlot(Seurat.analyse.J0, features = c("ZFP36L2","SPINK2","SELL","AVP","ERG"), pt.size = 0.5, slot = "data",label = FALSE,reduction = "umap", cols = c("lightgrey","#FF0000"), combine = TRUE) & list(labs(
  x="-UMAP2",
  y="UMAP1"),
  theme_bw(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  &
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) + theme(legend.title=element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), 
                                                axis.title = element_blank()))

#Granulo clusters
FeaturePlot(Seurat.analyse.J0, features = c("SERPINB1","PRSS57","ATF3","PLAC8","FCER1A","APOC1","CEBPD"), pt.size = 0.5, slot = "data",label = FALSE,reduction = "umap", cols = c("lightgrey","#FF0000"), combine = TRUE) & list(labs(
  x="-UMAP2",
  y="UMAP1"),
  theme_bw(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  &
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) + theme(legend.title=element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), 
                                                axis.title = element_blank()))

#Lympho clusters
FeaturePlot(Seurat.analyse.J0, features = c("SEMA3C","LTB","CD79A","MS4A1","CD3E","CD79B"), pt.size = 0.5, slot = "data",label = FALSE,reduction = "umap", cols = c("lightgrey","#FF0000"), combine = TRUE) & list(labs(
  x="-UMAP2",
  y="UMAP1"),
  theme_bw(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  &
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) + theme(legend.title=element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), 
                                                axis.title = element_blank()))

#rename clusters 
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "0" = "MPP_LMPP",
                                  "1" = "MPP1",
                                  "2" = "CMP",
                                  "3" = "MEP_proMK",
                                  "4" = "HSC_MPP2",
                                  "5" = "CLP_proB1",
                                  "6" = "LMPP1",
                                  "7" = "CLP_proB2",
                                  "8" = "CLP1",
                                  "9" = "HSC_MPP1",
                                  "10" = "MPP2",
                                  "11" = "proT",
                                  "12" = "proErythro",
                                  "13" = "LMPP2",
                                  "14" = "Neutro",
                                  "15" = "CLP2",
                                  "16" = "GMP")

#Analyse of gene expression in de most primitives clusters (1/4/9/10), between the different culture conditions 
Clusters_prim<-subset(seurat_integrated,idents=c("MPP_LMPP","MPP1","HSC_MPP2","HSC_MPP1"))

DimPlot(Clusters_prim,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        repel = TRUE)

Clusters_prim<-RenameIdents(Clusters_prim,"MPP_LMPP"="Primcells",
                            "MPP1"="Primcells",
                            "HSC_MPP2"="Primcells",
                            "HSC_MPP1"="Primcells" )

Clusters_prim$echantillon <- paste(Idents(Clusters_prim), Clusters_prim$orig.ident, sep = "_")
Clusters_prim$celltype <- Idents(Clusters_prim)

Idents(Clusters_prim) <- "echantillon"

Analyse_prim_cellsCTRL<-FindMarkers(Clusters_prim,ident.1 = "Primcells_CTRL6h",ident.2 = "Primcells_J0",slot = "data",logfc.threshold = 0.45, min.pct = 0.5, only.pos = FALSE)
Analyse_prim_cellsJNK<-FindMarkers(Clusters_prim,ident.1 = "Primcells_JNK6h",ident.2 = "Primcells_J0",slot = "data",logfc.threshold = 0.45, min.pct = 0.5, only.pos = FALSE)
Analyse_prim_cellsALA<-FindMarkers(Clusters_prim,ident.1 = "Primcells_ALA6h",ident.2 = "Primcells_J0",slot = "data",logfc.threshold = 0.45, min.pct = 0.5, only.pos = FALSE)
Analyse_prim_cellsJA<-FindMarkers(Clusters_prim,ident.1 = "Primcells_J+A6h",ident.2 = "Primcells_J0",slot = "data",logfc.threshold = 0.45, min.pct = 0.5, only.pos = FALSE)
