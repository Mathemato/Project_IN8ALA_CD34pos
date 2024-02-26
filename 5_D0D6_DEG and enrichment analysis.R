library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(ggplot2)
annotations <- readRDS("C:/Data.path/annotations.rds")

# 6H test samples versus 6H control sample GSEAnalysis----
# DEG calculation for each 6H TEST versus CTRL culture condition, in the most primitive clusters (1/3/8/12) = Primcells
Clusters_prim<-subset(seurat_integrated,idents=c("1","3","8","12"))

DimPlot(Clusters_prim,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        repel = TRUE)

Clusters_prim<-RenameIdents(Clusters_prim,"1"="Primcells",
                            "3"="Primcells",
                            "8"="Primcells",
                            "12"="Primcells")

Clusters_prim$sample <- paste(Idents(Clusters_prim), Clusters_prim$orig.ident, sep = "_")
Clusters_prim$celltype <- Idents(Clusters_prim)

Idents(Clusters_prim) <- "sample"

saveRDS(Clusters_prim,"C:/Data.path/Clusters_prim.rds")

# JNK_6h vs CTRL_6h:
DefaultAssay(Clusters_prim)<-"RNA"
Analyse_dseq2_prim_cellsJNK<-FindMarkers(Clusters_prim,test.use = "DESeq2",ident.1 = "Primcells_JNK6h",ident.2 = "Primcells_CTRL6h",min.pct = 0.5,slot = "data",only.pos = FALSE)
Analyse_dseq2_prim_cellsJNK<-tibble::rownames_to_column(Analyse_dseq2_prim_cellsJNK,"gene")
Analyse_dseq2_prim_cellsJNK <- inner_join(Analyse_dseq2_prim_cellsJNK, annotations, by=c("gene"="gene_name"))
write.table(x = Analyse_dseq2_prim_cellsJNK, file = "C:/Data.path/Analyse_dseq2_prim_cellsJNK.csv", sep = ",", row.names = T)

JNK_6H_FC<- Analyse_dseq2_prim_cellsJNK$avg_log2FC
names(JNK_6H_FC) <- Analyse_dseq2_prim_cellsJNK$gene_id
JNK_6H_FC<-na.omit(JNK_6H_FC)
JNK_6H_FC = sort(JNK_6H_FC, decreasing = TRUE)

Analyse_prim_cellsJNK_GSE <- gseGO(geneList=JNK_6H_FC, 
                                    ont ="ALL", 
                                    keyType = "ENSEMBL", 
                                    nPerm = 10000, 
                                    minGSSize = 3, 
                                    maxGSSize = 800, 
                                    pvalueCutoff = 0.05, 
                                    verbose = TRUE, 
                                    OrgDb = org.Hs.eg.db, 
                                    pAdjustMethod = "BH")
saveRDS(Analyse_prim_cellsJNK_GSE,"C:/Data.path/Analyse_Prim_cellsJNK_GSE.rds")
Analyse_prim_cellsJNK_GSE_result <- Analyse_prim_cellsJNK_GSE@result
view(Analyse_prim_cellsJNK_GSE_result)
write.csv(Analyse_prim_cellsJNK_GSE_result, "C:/Data.path/Prim_JNK_6h_GSE_result.csv", quote=F)

categoriesJNK<-c("lactate metabolic process","carbohydrate catabolique process","lactate dehydrogenase activity","pyruvate metabolic process", "glycolytic process","carbohydrate catabolic process", "NADH dehydrogenase complex assembly",
                 "mitochondrial respiratory chain complex I assembly","mitochondrial electron transport, NADH to ubiquinone")
dotplot(Analyse_prim_cellsJNK_GSE, showCategory=categoriesJNK, split=c(".sign")) + facet_grid(.~.sign)



# ALA_6h vs CTRL_6h:
DefaultAssay(Clusters_prim)<-"RNA"
Analyse_dseq2_prim_cellsALA<-FindMarkers(Clusters_prim,test.use = "DESeq2",ident.1 = "Primcells_ALA6h",ident.2 = "Primcells_CTRL6h",min.pct = 0.5,slot = "data",only.pos = FALSE)
Analyse_dseq2_prim_cellsALA<-tibble::rownames_to_column(Analyse_dseq2_prim_cellsALA,"gene")
Analyse_dseq2_prim_cellsALA <- inner_join(Analyse_dseq2_prim_cellsALA, annotations, by=c("gene"="gene_name"))
write.table(x = Analyse_dseq2_prim_cellsALA, file = "C:/Data.path/Analyse_dseq2_prim_cellsALA.csv", sep = ",", row.names = T)

ALA_6H_FC<- Analyse_dseq2_prim_cellsALA$avg_log2FC
names(ALA_6H_FC) <- Analyse_dseq2_prim_cellsALA$gene_id
ALA_6H_FC<-na.omit(ALA_6H_FC)
ALA_6H_FC = sort(ALA_6H_FC, decreasing = TRUE)

Analyse_prim_cellsALA_GSE <- gseGO(geneList=ALA_6H_FC, 
                                   ont ="ALL", 
                                   keyType = "ENSEMBL", 
                                   nPerm = 10000, 
                                   minGSSize = 3, 
                                   maxGSSize = 800, 
                                   pvalueCutoff = 0.05, 
                                   verbose = TRUE, 
                                   OrgDb = org.Hs.eg.db, 
                                   pAdjustMethod = "BH")
saveRDS(Analyse_prim_cellsALA_GSE,"C:/Data.path/Analyse_prim_cellsALA_GSE.rds")

Analyse_prim_cellsALA_GSE_result <- Analyse_prim_cellsALA_GSE@result
view(Analyse_prim_cellsALA_GSE_result)
write.csv(Analyse_prim_cellsALA_GSE_result, "C:/Data.path/prim_ALA_6h_GSE_result.csv", quote=F)

categoriesALA<-c("lactate metabolic process","carbohydrate catabolique process","lactate dehydrogenase activity","pyruvate metabolic process", "glycolytic process","carbohydrate catabolic process", "NADH dehydrogenase complex assembly",
                 "mitochondrial respiratory chain complex I assembly","mitochondrial electron transport, NADH to ubiquinone")
dotplot(Analyse_prim_cellsALA_GSE, showCategory=categoriesALA, split=c(".sign")) + facet_grid(.~.sign)

# IN8ALA_6h vs CTRL_6h:
DefaultAssay(Clusters_prim)<-"RNA"
Analyse_dseq2_prim_cellsIN8ALA<-FindMarkers(Clusters_prim,test.use = "DESeq2",ident.1 = "Primcells_IN8ALA6h",ident.2 = "Primcells_CTRL6h",min.pct = 0.5,slot = "data",only.pos = FALSE)
Analyse_dseq2_prim_cellsIN8ALA<-tibble::rownames_to_column(Analyse_dseq2_prim_cellsIN8ALA,"gene")
Analyse_dseq2_prim_cellsIN8ALA <- inner_join(Analyse_dseq2_prim_cellsIN8ALA, annotations, by=c("gene"="gene_name"))
write.table(x = Analyse_dseq2_prim_cellsIN8ALA, file = "C:/Data.path/Analyse_dseq2_prim_cellsIN8ALA_GSE.csv", sep = ",", row.names = T)

IN8ALA_6H_FC<- Analyse_dseq2_prim_cellsIN8ALA$avg_log2FC
names(IN8ALA_6H_FC) <- Analyse_dseq2_prim_cellsIN8ALA$gene_id
IN8ALA_6H_FC<-na.omit(IN8ALA_6H_FC)
IN8ALA_6H_FC = sort(IN8ALA_6H_FC, decreasing = TRUE)

Analyse_prim_cellsIN8ALA_GSE <- gseGO(geneList=IN8ALA_6H_FC, 
                                   ont ="ALL", 
                                   keyType = "ENSEMBL", 
                                   nPerm = 10000, 
                                   minGSSize = 3, 
                                   maxGSSize = 800, 
                                   pvalueCutoff = 0.05, 
                                   verbose = TRUE, 
                                   OrgDb = org.Hs.eg.db, 
                                   pAdjustMethod = "BH")
saveRDS(Analyse_prim_cellsIN8ALA_GSE,"C:/Data.path/Analyse_prim_cellsIN8ALA_GSE.rds")

Analyse_prim_cellsIN8ALA_GSE_result <- Analyse_prim_cellsIN8ALA_GSE@result
view(Analyse_prim_cellsIN8ALA_GSE_result)
write.csv(Analyse_prim_cellsIN8ALA_GSE_result, "C:/Data.path/prim_IN8ALA_6h_GSE_result.csv", quote=F)

categoriesIN8ALA<-c("lactate metabolic process","carbohydrate catabolique process","lactate dehydrogenase activity","pyruvate metabolic process", "glycolytic process","carbohydrate catabolic process", "NADH dehydrogenase complex assembly",
                    "mitochondrial respiratory chain complex I assembly","mitochondrial electron transport, NADH to ubiquinone")
dotplot(Analyse_prim_cellsIN8ALA_GSE, showCategory=categoriesIN8ALA, split=c(".sign")) + facet_grid(.~.sign)

