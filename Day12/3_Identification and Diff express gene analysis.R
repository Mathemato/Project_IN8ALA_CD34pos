# Find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(MergeJ12_norm)<-"RNA"
NormalizeData(MergeJ12_norm)
MarqueursJ12 <- FindAllMarkers(object = MergeJ12_norm, 
                               only.pos = FALSE,
                               logfc.threshold = 0.3,
                               min.pct = 0.5)
MarqueursJ12<- inner_join(x = MarqueursJ12, 
                          y = annotations[, c("gene_name", "description")],
                          by = c("gene" = "gene_name")) %>% unique()

MarqueursJ12<- MarqueursJ12[ , c(6, 7, 2:4, 1, 5,8)]

MarqueursJ12<- MarqueursJ12 %>%
  dplyr::arrange(cluster)

#Volcanoplot----

#Add NA column
MarqueursJ12$diffexpressed<- "NO"

#Identify overexpressed genes as "up"
MarqueursJ12$diffexpressed[MarqueursJ12$avg_log2FC > 0.5 & MarqueursJ12$p_val_adj < 0.05] <- "UP"

#Identify downregulated genes as "down"
MarqueursJ12$diffexpressed[MarqueursJ12$avg_log2FC < -0.5 & MarqueursJ12$p_val_adj < 0.05] <- "DOWN"
MarqueursJ12$delabel <- NA
MarqueursJ12$delabel[MarqueursJ12$diffexpressed != "NO"] <- MarqueursJ12$gene[MarqueursJ12$diffexpressed != "NO"]

write.table(MarqueursJ12, file = "C:/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/CTRLJ12 + JA J12/MarqueursJ12.xls", sep="\t", quote = TRUE, dec = ",", row.names = FALSE, col.names = TRUE)          

#TOP5 genes
topJ12<-MarqueursJ12%>%group_by(cluster)%>%top_n(5,avg_log2FC)
DoHeatmap(MergeJ12_norm,features=topJ12$gene,label = TRUE , group.by = c("seurat_clusters"),group.bar = TRUE,size = 4, disp.min = 0 ) + scale_fill_gradientn(colors = c("black", "Yellow"))+ theme(text = element_text(size = 10))
write.table(topJ12, file = "C:/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/CTRLJ12 + JA J12/TOP5J12.xls", sep="\t", quote = TRUE, dec = ",", row.names = FALSE, col.names = TRUE)          

#Feature plot 

#cluster 0 HSC_MPP
FeaturePlot(MergeJ12_norm, features = c("HOXA9","SPINK2","SELL"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluster 1 GMP
FeaturePlot(MergeJ12_norm, features = c("ELANE","AZU1","PRTN3"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluster 2 PRO_NEUTRO
FeaturePlot(MergeJ12_norm, features = c("S100A8","S100A9","FUCA1"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluster 3 MEP
FeaturePlot(MergeJ12_norm, features = c("HBD","PLEK","GATA1"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluster 4 MPP_3
FeaturePlot(MergeJ12_norm, features = c("C1QTNF4","GINS2","MCM4"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluster 5 MPP_1
FeaturePlot(MergeJ12_norm, features = c("SPINK2","SELL","TOP2A"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluster 7 MPP_2
FeaturePlot(MergeJ12_norm, features = c("SPINK2","SELL","IGLL1"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluster 6 CPA
FeaturePlot(MergeJ12_norm, features = c("HLA-DPA1","HLA-DQB1","HLA-DPB1"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluter 9 PRO_ERY 
FeaturePlot(MergeJ12_norm, features = c("HBB","KLF1","TFRC"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluster 8  PRE_GRANULO
FeaturePlot(MergeJ12_norm, features = c("TPSAB1","MS4A2","CPA3"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
#cluster 10 PRO_GRANULO
FeaturePlot(MergeJ12_norm, features = c("CLC","LMO4","PRSS57"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")
# cluster 11 MEP_MK
FeaturePlot(MergeJ12_norm, features = c("PF4","PLEK","GP9"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")

#cluster prim 
FeaturePlot(MergeJ12_norm, features = c("HOXA9","SPINK2","SELL","TOP2A","IGLL1","C1QTNF4","GINS2","MCM4"), keep.scale = c("all"), pt.size = 0.5,label = TRUE,reduction = "umap")

#Rename clusters 
MergeJ12_norm <- RenameIdents(object = MergeJ12_norm, 
                              "0" = "HSC_MPP",
                              "1" = "GMP",
                              "2" = "PRO_NEUTRO",
                              "3" = "MEP",
                              "4" = "PROG3",
                              "5" = "PROG1",
                              "6" = "CPA",
                              "7" = "PROG2",
                              "8" = "PRE_GRANULO",
                              "9" = "PRO_ERYTHRO",
                              "10" = "PRO_GRANULO",
                              "11" = "MEP_MK")
DimPlot(MergeJ12_norm,label = TRUE, repel = TRUE,label.size = 6)

FeaturePlot(MergeJ12_norm, features = c("HOXA9","SPINK2","SELL","TOP2A","IGLL1","C1QTNF4","GINS2","MCM4"), pt.size = 0.5, slot = "data",label = FALSE,reduction = "umap", cols = c("lightgrey","#FF0000"), combine = TRUE) & list(labs(
  x="-UMAP2",
  y="UMAP1"),
  theme_bw(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  &
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) + theme(legend.title=element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), 
                                                axis.title = element_blank()))

FeaturePlot(MergeJ12_norm, features = c("ELANE","AZU1","PRTN3","CLC","LMO4","PRSS57","TPSAB1","MS4A2","CPA3","S100A8","S100A9","FUCA1"), pt.size = 0.5, slot = "data",label = FALSE,reduction = "umap", cols = c("lightgrey","#FF0000"), combine = TRUE) & list(labs(
  x="-UMAP2",
  y="UMAP1"),
  theme_bw(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  &
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) + theme(legend.title=element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), 
                                                axis.title = element_blank()))

FeaturePlot(MergeJ12_norm, features = c("HBD","PLEK","GATA1","PF4","GP9","HBB","KLF1","TFRC"), pt.size = 0.5, slot = "data",label = FALSE,reduction = "umap", cols = c("lightgrey","#FF0000"), combine = TRUE) & list(labs(
  x="-UMAP2",
  y="UMAP1"),
  theme_bw(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  &
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) + theme(legend.title=element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), 
                                                axis.title = element_blank()))


FeaturePlot(MergeJ12_norm, features = c("HLA-DPA1","HLA-DQB1","HLA-DPB1"), pt.size = 0.5, slot = "data",label = FALSE,reduction = "umap", cols = c("lightgrey","#FF0000"), combine = TRUE) & list(labs(
  x="-UMAP2",
  y="UMAP1"),
  theme_bw(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  &
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) + theme(legend.title=element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), 
                                                axis.title = element_blank()))


MergeJ12_norm$echantillon <- paste(Idents(MergeJ12_norm), MergeJ12_norm$condition, sep = "_")
MergeJ12_norm$celltype <- Idents(MergeJ12_norm)
Idents(MergeJ12_norm) <- "echantillon"
MergeJ12_norm$echantillon
C0_CTRL_JA<-FindMarkers(MergeJ12_norm,ident.1 ="HSC_MPP_CTRL",ident.2 = "HSC_MPP_JNK+ALA", logfc.threshold = 0.2, min.pct = 0.3, only.pos = FALSE )
write.table(C0_CTRL_JA, file = "C:/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/CTRLJ12 + JA J12/C0_CTRL_JA.xls", sep="\t", quote = TRUE, dec = ",", row.names = FALSE, col.names = TRUE)          


#Diff express genes in cluster C0

options(future.globals.maxSize = 10000 * 1024^2)
DefaultAssay(MergeJ12_norm)<-"RNA"
MergeJ12_norm[["RNA"]]@counts <- as.matrix(MergeJ12_norm[["RNA"]]@counts)+1
Prim12h_IN8ALAvsCTRL<-FindMarkers(MergeJ12_norm, test.use = "DESeq2", ident.1 = c("HSC_MPP_JNK+ALA_JNK+ALA"), ident.2 = c("HSC_MPP_CTRL_CTRL"), slot = "data")
write.table(x = Prim12h_IN8ALAvsCTRL, file = "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Prim12h_IN8ALAvsCTRLDESeq2.csv", sep = ",", row.names = T)
write.table(x = dfj0, file = "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Analyse 6h + J0 intégrée/DESeq2J0.xls", sep = ",", row.names = T)

#Add gene name of dataframe annotations 
Prim12h_IN8ALAvsCTRLdseq2<-tibble::rownames_to_column(Prim12h_IN8ALAvsCTRL,"gene")
Prim12h_IN8ALAvsCTRLdseq2 <- inner_join(Prim12h_IN8ALAvsCTRLdseq2, annotations, by=c("gene"="gene_name"))
write.table(x = Prim12h_IN8ALAvsCTRLdseq2, file = "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Analyse 6h + J0 intégrée/Prim12h_IN8ALAvsCTRLdseq2.xls", sep = ",", row.names = T)

#Diff exp gene (DEG)

AllGenes<-as.character(Prim12h_IN8ALAvsCTRLdseq2$gene_id)

Prim12h_IN8ALAvsCTRLdseq2_deg<-dplyr::filter(Prim12h_IN8ALAvsCTRLdseq2, p_val_adj < 0.05)
Prim12h_IN8ALAvsCTRLdseq2_deg_genes<-as.character(Prim12h_IN8ALAvsCTRLdseq2_deg$gene_id)

Prim12h_IN8ALAvsCTRLdseq2_GO<-enrichGO(gene = Prim12h_IN8ALAvsCTRLdseq2_deg_genes, 
                                       universe= AllGenes, 
                                       keyType = "ENSEMBL", 
                                       OrgDb =  org.Hs.eg.db,
                                       ont =  "ALL",
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 0.05,
                                       readable = T)

gsea_go_Prim12h_IN8ALAvsCTRL<-data.frame(Prim12h_IN8ALAvsCTRLdseq2_GO)
write.csv(gsea_go_Prim12h_IN8ALAvsCTRL, "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Résultats_enrichgo_Prim12h_IN8ALAvsCTRL.csv")

# Dotplot
dotplot(Prim12h_IN8ALAvsCTRLdseq2_GO, showCategory=20)

# For KEGG pathway enrichment using the gseKEGG() function, we need to convert
# id types. We can use the bitr function for this (included in clusterProfiler). 
# It is normal for this call to produce some messages / warnings.
# In the bitr function, the param fromType should be the same as keyType from 
# gseGO function above (the annotation source). 
# toType in the bitr function has to be one of the available options from keyTypes
# (org.Dm.eg.db) and must map to one of ‘kegg’, ‘ncbi-geneid’, ‘ncib-proteinid’ or ‘uniprot’ 
# because gseKEGG() only accepts one of these 4 options as it’s keytype parameter. 
# In the case of org.Dm.eg.db, none of those 4 types are available, but ‘ENTREZID’ 
# are the same as ncbi-geneid for org.Dm.eg.db so we use this for toType.
# As our intial input, we use original_gene_list which we created above.

id<-bitr(Prim12h_IN8ALAvsCTRLdseq2$gene_id,fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
Prim12h_IN8ALAvsCTRLdseq2_ID<-inner_join(Prim12h_IN8ALAvsCTRLdseq2,id,by=c("gene_id"="ENSEMBL"), multiple = "all")
Prim12h_IN8ALAvsCTRLdseq2_ID<-Prim12h_IN8ALAvsCTRLdseq2_ID[which(duplicated(Prim12h_IN8ALAvsCTRLdseq2_ID$ENTREZID)== F), ]

GSEA_FC <- Prim12h_IN8ALAvsCTRLdseq2_ID$avg_log2FC

names(GSEA_FC) <- Prim12h_IN8ALAvsCTRLdseq2_ID$ENTREZID

GSEA_FC<-na.omit(GSEA_FC)

GSEA_FC = sort(GSEA_FC, decreasing = TRUE)

Prim12h_IN8ALAvsCTRLdseq2_KEGG<- gseKEGG(geneList = GSEA_FC,
                                         keyType = "kegg",
                                         organism = "hsa",
                                         nPerm = 1000,
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "none")

Prim12h_IN8ALAvsCTRLdseq2_KEGG_results <- Prim12h_IN8ALAvsCTRLdseq2_KEGG@result
View(Prim12h_IN8ALAvsCTRLdseq2_KEGG_results)
write.csv(Prim12h_IN8ALAvsCTRLdseq2_KEGG_results, "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Prim12h_IN8ALAvsCTRLdseq2KEGG.csv", quote=F)

#Visualization KEGG analysis 
gseaplot(Prim12h_IN8ALAvsCTRLdseq2_KEGG, geneSetID = 'hsa00190')
browseKEGG(dfj0dseq2_gseaKEGG, 'hsa00190')

pathview(gene.data = GSEA_FC,
         pathway.id = "hsa04010",
         species = "hsa",
         limit = list(gene = 2, # value gives the max/min limit for foldchanges
                      cpd = 1))


# Dotplot
dotplot(dfj0dseq2_KEGG, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

# GSEA using gene sets associated with BP Gene Ontology terms
GSEA_FC_gsego <- Prim12h_IN8ALAvsCTRLdseq2_ID$avg_log2FC
names(GSEA_FC_gsego) <- Prim12h_IN8ALAvsCTRLdseq2_ID$gene_id
GSEA_FC_gsego<-na.omit(GSEA_FC_gsego)
GSEA_FC_gsego = sort(GSEA_FC_gsego, decreasing = TRUE)

Prim12h_IN8ALAvsCTRLdseq2_GSE <- gseGO(geneList=GSEA_FC_gsego, 
                                       ont ="ALL", 
                                       keyType = "ENSEMBL", 
                                       nPerm = 10000, 
                                       minGSSize = 3, 
                                       maxGSSize = 800, 
                                       pvalueCutoff = 0.05, 
                                       verbose = TRUE, 
                                       OrgDb = org.Hs.eg.db, 
                                       pAdjustMethod = "BH")

Prim12h_IN8ALAvsCTRLdseq2_GSE_result <- Prim12h_IN8ALAvsCTRLdseq2_GSE@result
view(Prim12h_IN8ALAvsCTRLdseq2_GSE_result)
write.csv(Prim12h_IN8ALAvsCTRLdseq2_GSE_result, "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Prim12h_IN8ALAvsCTRLdseq2_GSE_result.csv", quote=F)


gseaplot(dfj0dseq2_GSE, geneSetID = 'GO:0007423')

categories<-c("cytoplasmic translation", "aerobic respiration", "maintenance of location in cell", "antioxidant activity", "cellular oxidant detoxification", "cellular detoxification", 
              "oxidoreductase activity", "oxydative phosphorylation", "hydrogen peroxide metabolic process", "regulation of actin cytoskeleton organization", "lymphocyte differentiation", "mononuclear cell differentiation",
              "hematopoiesis", "myeloid cell homeostasis", "erythrocyte homeostasis", "reactiv oxygen species metabolic process", "actin polymerization or depolymerization", "regulation of actomyosin structure organization",
              "negative regulation of stress fiber assembly", "negative regulation of actin filament bundle assembly", "valine biosynthetic process", "short-chain fatty acid metabolic process")
dotplot(Prim12h_IN8ALAvsCTRLdseq2_GSE, showCategory=categories, split=c(".sign")) + facet_grid(.~.sign)
dotplot(dfj0dseq2_GSE, showCategory=10, split="ONTOLOGY") + facet_grid(.~.sign)
