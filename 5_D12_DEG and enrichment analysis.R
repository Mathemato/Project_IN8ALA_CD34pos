# Analysis of differential gene expression in primitive cells (cluster1) between IN8ALA and CTRL conditions in D12 dataset.

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

annotations <- readRDS("C:/Data.path/annotations.rds")

# DEG between IN8ALA vs CTRL at D12 in primitive cluster (cluster 1)
MergeD12_norm$echantillon <- paste(Idents(MergeD12_norm), MergeD12_norm$condition, sep = "_")
MergeD12_norm$celltype <- Idents(MergeD12_norm)
Idents(MergeD12_norm) <- "echantillon"
MergeD12_norm$echantillon
options(future.globals.maxSize = 10000 * 1024^2)
DefaultAssay(MergeD12_norm)<-"RNA"
PrimD12_IN8ALAvsCTRLdseq2<-FindMarkers(MergeD12_norm, test.use = "DESeq2", ident.1 = c("1_IN8ALA"), ident.2 = c("1_CTRL"), min.pct = 0.5, slot = "data")
PrimD12_IN8ALAvsCTRLdseq2<-tibble::rownames_to_column(PrimD12_IN8ALAvsCTRLdseq2,"gene")
PrimD12_IN8ALAvsCTRLdseq2 <- inner_join(PrimD12_IN8ALAvsCTRLdseq2, annotations, by=c("gene"="gene_name"))
write.table(x = PrimD12_IN8ALAvsCTRLdseq2, file = "/Data.path/PrimD12_IN8ALAvsCTRLdseq2.csv", sep = ",", row.names = T)
saveRDS(PrimD12_IN8ALAvsCTRLdseq2,"C:/Data.path/PrimD12_IN8ALAvsCTRLdseq2.rds")

# GSEA analysis between IN8ALA vs CTRL in primtive cluster (cluster 1) 
PrimD12_IN8ALAvsCTRLdseq2_FC<- PrimD12_IN8ALAvsCTRLdseq2$avg_log2FC
names(PrimD12_IN8ALAvsCTRLdseq2_FC) <- PrimD12_IN8ALAvsCTRLdseq2$gene_id
PrimD12_IN8ALAvsCTRLdseq2_FC<-na.omit(PrimD12_IN8ALAvsCTRLdseq2_FC)
PrimD12_IN8ALAvsCTRLdseq2_FC = sort(PrimD12_IN8ALAvsCTRLdseq2_FC, decreasing = TRUE)
write.csv(PrimD12_IN8ALAvsCTRLdseq2_FC, "C:/Data.path/PrimD12_IN8ALAvsCTRLdseq2_FC.csv", quote=F)

PrimD12_IN8ALAvsCTRLdseq2_GSE <- gseGO(geneList=PrimD12_IN8ALAvsCTRLdseq2_FC, 
                       ont ="ALL", 
                       keyType = "ENSEMBL", 
                       nPerm = 10000, 
                       minGSSize = 3, 
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = org.Hs.eg.db, 
                       pAdjustMethod = "BH")
saveRDS(PrimD12_IN8ALAvsCTRLdseq2_GSE,"C:/Data.path/PrimD12_IN8ALAvsCTRLdseq2_GSE.rds")

PrimD12_IN8ALAvsCTRLdseq2_GSE_result <- PrimD12_IN8ALAvsCTRLdseq2_GSE@result
view(PrimD12_IN8ALAvsCTRLdseq2_GSE_result)
write.csv(PrimD12_IN8ALAvsCTRLdseq2_GSE_result, "C:/Data.path/PrimD12_IN8ALAvsCTRLdseq2_GSE_result.csv", quote=F)

D12categories<-c("spindle assembly", "spindle organization", "actin filament polymerization", "cell junction organization", "actin filament organization", 
              "regulation of cytoskeleton organization", "actin filament-based process", "actin cytoskeleton organization")
dotplot(PrimD12_IN8ALAvsCTRLdseq2_GSE, showCategory=D12categories, split=c(".sign")) + facet_grid(.~.sign)

# Targeted analysis of JNK and WNT pathways in primtive cluster (cluster 1)
D12Primitivepop<-subset(MergeD12_norm,subset = seurat_clusters == '1')
DefaultAssay(D12Primitivepop)<-"RNA"
JNKPATHWAYS <-c("ACTB","ACTG1","CFL1","PFN1","TUBA1A", "VIM","RHOA","ROCK1","CDC42","JUN","FOS","ELK1")
dotJNKPATHWAYS<-DotPlot(D12Primitivepop,
                        features = JNKPATHWAYS, group.by = c("condition"),
                        dot.scale = 10,cols = c("green","red")) + RotatedAxis()
dotJNKPATHWAYS

WNTPATHWAYS <-c("CTNNB1","TCF4","SMAD4","TP53","HSF1","MYC","TLE4","HDAC2")
dotWNTPATHWAYS<-DotPlot(D12Primitivepop,
                        features = WNTPATHWAYS, group.by = c("condition"),
                        dot.scale = 10,cols = c("green","red")) + RotatedAxis()
dotWNTPATHWAYS

