# Normalization, PCA and UMAP of D12 merge dataset

# Normalize the counts
D12seurat_phase <- NormalizeData(Filtered_MergeD12)

# Analysis of cell cycle and mitoRatio variables to eventually regress

# Analysis of cell cyle variable
# Load cell cycle markers
load("C:/Data.path/cycle.rda")
# Score cells for cell cycle
D12seurat_phase <- CellCycleScoring(D12seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)  

# Identify the most variable genes
D12seurat_phase <- FindVariableFeatures(D12seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 3000, verbose = FALSE)

# Scale the counts
D12seurat_phase <- ScaleData(D12seurat_phase)

# Perform PCA
D12seurat_phase <- RunPCA(D12seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(D12seurat_phase,
        reduction = "pca",
        group.by= "condition",
        split.by = "Phase")

# Check quartile values
summary(D12seurat_phase@meta.data$mitoRatio)

# Analysis of mitoRatio variable
# Turn mitoRatio into categorical factor vector based on quartile values
D12seurat_phase@meta.data$mitoFr <- cut(D12seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.03143, 0.04146, 0.05028, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))
# Plot the PCA 
DimPlot(D12seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")

saveRDS(D12seurat_phase,"C:/Data.path/D12seurat_phase.rds")

# Run SCTransform
MergeD12_norm<-SCTransform(D12seurat_phase,vars.to.regress = c("S.Score","G2M.Score"))

# Run PCA
MergeD12_norm <- RunPCA(object = MergeD12_norm)

# Plot PCA
PCAPlot(MergeD12_norm,
        split.by = "orig.ident") 

# Run UMAP
MergeD12_norm <- RunUMAP(MergeD12_norm, 
                         dims = 1:40,
                         reduction = "pca")
# Plot UMAP                             
DimPlot(MergeD12_norm)  
MergeD12_norm$orig.ident<-factor(MergeD12_norm$orig.ident, levels = c("CTRLD121","CTRLD122","IN8ALAD12"))

DimPlot(MergeD12_norm,
        split.by = "orig.ident") 
saveRDS(MergeD12_norm,"C:/Data.path/MergeD12_norm.rds")

