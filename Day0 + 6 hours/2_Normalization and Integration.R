#Normalize the counts
seurat_phase <- NormalizeData(Filtered_Mergeso)

#Load cell cycle markers
load("/Users/13700976/Desktop/Analyses single cells Mathilde/cycle.rda")
# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

#View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)  

#Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 3000, verbose = FALSE)

#Scale the counts
seurat_phase <- ScaleData(seurat_phase)

#Perform PCA
seurat_phase <- RunPCA(seurat_phase)

#Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

#Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

#Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))
# Plot the PCA 
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")


#Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "orig.ident")

split_seurat <- split_seurat[c("J0", "CTRL6h","JNK6h", "ALA6h", "J+A6h")]


#increase the size of the RAM 
options(future.globals.maxSize = 10000 * 1024^2)

#SCTransform on all samples
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score"))}

#Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 2000) 
#Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
#Integrate across conditions
options(future.globals.maxSize = 100000 * 1024^2)

seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

#Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

#Plot PCA
PCAPlot(seurat_integrated,
        split.by = "orig.ident") 

#Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

#Plot UMAP                             
DimPlot(seurat_integrated)  
seurat_integrated$orig.ident<-factor(seurat_integrated$orig.ident, levels = c("J0","CTRL6h","JNK6h","ALA6h","J+A6h"))

DimPlot(seurat_integrated,
        split.by = "orig.ident") 