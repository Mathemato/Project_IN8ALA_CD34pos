#Create Sub-data frames 
Seurat.analyse.J0<-subset(seurat_integrated,subset = orig.ident == 'J0')

Seurat.analyse.CTRL6h<-subset(seurat_integrated,subset = orig.ident == 'CTRL6h')
Seurat.analyse.JNK6h<-subset(seurat_integrated,subset = orig.ident == 'JNK6h')
Seurat.analyse.ALA6h<-subset(seurat_integrated,subset = orig.ident == 'ALA6h')
Seurat.analyse.JA6h<-subset(seurat_integrated,subset = orig.ident == 'J+A6h')

#Genes differentially expressed in clusters HSC_MPP1 & HSC_MPP2 at J0
#Select default dataset, 
#add +1 to the results because of log(0)
#Analyze all cluster marker genes relative to others in the J0 condition
options(future.globals.maxSize = 10000 * 1024^2)
DefaultAssay(Seurat.analyse.J0)<-"RNA"
Seurat.analyse.J0[["RNA"]]@counts <- as.matrix(Seurat.analyse.J0[["RNA"]]@counts)+1
dfj0<-FindMarkers(Seurat.analyse.J0, test.use = "DESeq2", ident.1 = c("9","4"), ident.2 = c("0","1","2","3","5","6","7","8","10","11","12","13","14","15","16"), slot = "data")
write.table(x = dfj0, file = "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Analyse 6h + J0 intégrée/DESeq2J0.csv", sep = ",", row.names = T)
write.table(x = dfj0, file = "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Analyse 6h + J0 intégrée/DESeq2J0.xls", sep = ",", row.names = T)

#add name gene column 
dfj0dseq2<-tibble::rownames_to_column(dfj0,"gene")
dfj0dseq2 <- inner_join(dfj0dseq2, annotations, by=c("gene"="gene_name"))
write.table(x = dfj0dseq2, file = "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Analyse 6h + J0 intégrée/dfj0dseq2.xls", sep = ",", row.names = T)

# List gene creation necessary for the analysis----

AllGenes<-as.character(dfj0dseq2$gene_id)

# Create DEG (diff exp genes) list
dfj0dseq2_deg<-dplyr::filter(dfj0dseq2, p_val_adj < 0.05)
dfj0dseq2_deg_genes<-as.character(dfj0dseq2_deg$gene_id)

# GSEA analysis with Gene ontology (GO)----
#Use the gseGO function for all the variations (+ and -) and enrichGO to get only the enriched pathway on the samples
# ont one of “BP” (biological process), “MF” (molecular function), “CC” (cellular component) or “ALL”
# pAdjustMethod one of “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none"

dfj0dseq2_GO<-enrichGO(gene = dfj0dseq2_deg_genes, 
                       universe= AllGenes, 
                       keyType = "ENSEMBL", 
                       OrgDb =  org.Hs.eg.db,
                       ont =  "ALL",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = T)

gsea_go_J0<-data.frame(dfj0dseq2_GO)
write.csv(gsea_go_J0, "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Analyse 6h + J0 intégrée/Résultats_enrichgo_J0.csv")

#GSE GO vizualisation----
# Dotplot
dotplot(dfj0dseq2_GO, showCategory=20)

#Create a netplot
dfj0dseq2_FC <- dfj0dseq2_deg$avg_log2FC
names(dfj0dseq2_FC) <- dfj0dseq2_deg$gene
dfj0dseq2_FC_list<-na.omit(dfj0dseq2_FC)
dfj0dseq2_FC_list = sort(dfj0dseq2_FC_list, decreasing = TRUE)
cnetplot(dfj0dseq2_GO, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=dfj0dseq2_FC_list, 
         vertex.label.font=6)

#To perform GSEA analysis of KEGG gene sets, 
#clusterProfiler requires the genes to be identified 
#using Entrez IDs for all genes in our results dataset. 
#We also need to remove the NA values and duplicates (due to gene ID conversion) prior to the analysis:

# Realisation GSEA à partir des gene sets des voies KEGG----

# Preparation des donnees
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

# Conversion of gene IDs for the gseKEGG function, in ENTREZID format
# Some genes will be lost because not all of them will be converted
# First convert all genes in ENTREZID format, into a new df
id<-bitr(dfj0dseq2$gene_id,fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)

# Then associate ENTREZID IDs with dfj0dseq2, via ENSEMBL IDs, 
# in a new df using the inner_join function
dfj0dseq2_ID<-inner_join(dfj0dseq2,id,by=c("gene_id"="ENSEMBL"), multiple = "all")
dfj0dseq2_ID<-dfj0dseq2_ID[which(duplicated(dfj0dseq2_ID$ENTREZID)== F), ]

# Vector creation with log2FoldChange data
GSEA_FC <- dfj0dseq2_ID$avg_log2FC

# Name the vector with "ENTER ids".
names(GSEA_FC) <- dfj0dseq2_ID$ENTREZID
GSEA_FC<-na.omit(GSEA_FC)

# Sort list in descending order (required for ClusterProfiler)
GSEA_FC = sort(GSEA_FC, decreasing = TRUE)

# Perform analysis and create gseaKEGG object
dfj0dseq2_KEGG<- gseKEGG(geneList = GSEA_FC,
                         keyType = "kegg",
                         organism = "hsa",
                         nPerm = 1000,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "none")

# Extract the result (modulated channels) and save it
dfj0dseq2_KEGG_results <- dfj0dseq2_KEGG@result
View(dfj0dseq2_KEGG_results)
write.csv(dfj0dseq2_KEGG_results, "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Analyse 6h + J0 intégrée/dfj0dseq2KEGG.csv", quote=F)

#Vizualisation ----

# View GSEA KEGG analysis on 1 KEGG gene set
gseaplot(dfj0dseq2_KEGG, geneSetID = 'hsa00190')

# View a KEGG track (send via web browser)
browseKEGG(dfj0dseq2_gseaKEGG, 'hsa00190')

# Creation of a loop to display modified tracks on a map
# first unload dplyr to avoid conflicts
detach("package:dplyr", unload=TRUE) 
get_kegg_plots <- function(x) {
  pathview(gene.data = kegg_gene_list, pathway.id = gseaKEGG_results$ID[x], species = "hsa", 
           limit = list(gene = 2, cpd = 1))
}

purrr::map(1:length(gseaKEGG_results$ID), get_kegg_plots)

# Pathview visualization of GSEA KEGG analysis within a KEGG map
pathview(gene.data = GSEA_FC,
         pathway.id = "hsa04010",
         species = "hsa",
         limit = list(gene = 2, # value gives the max/min limit for foldchanges
                      cpd = 1))


# Dotplot
dotplot(dfj0dseq2_KEGG, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

# GSEA using gene sets associated with BP Gene Ontology terms
# Creation d'un vecteur avec les donnees de log2FoldChange
GSEA_FC_gsego <- dfj0dseq2_ID$avg_log2FC
names(GSEA_FC_gsego) <- dfj0dseq2_ID$gene_id
GSEA_FC_gsego<-na.omit(GSEA_FC_gsego)
GSEA_FC_gsego = sort(GSEA_FC_gsego, decreasing = TRUE)

dfj0dseq2_GSE <- gseGO(geneList=GSEA_FC_gsego, 
                       ont ="ALL", 
                       keyType = "ENSEMBL", 
                       nPerm = 10000, 
                       minGSSize = 3, 
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = org.Hs.eg.db, 
                       pAdjustMethod = "BH")

dfj0dseq2_GSE_result <- dfj0dseq2_GSE@result
view(dfj0dseq2_GSE_result)
write.csv(dfj0dseq2_GSE_result, "/Users/13700976/Desktop/Analyses single cells Mathilde/Analyses multiples/Analyse 6h + J0 intégrée/dfj0dseq2_GSE_result.csv", quote=F)


gseaplot(dfj0dseq2_GSE, geneSetID = 'GO:0007423')

categories<-c("cytoplasmic translation", "aerobic respiration", "aerobic electron transport chain", "oxidative phosphorylation",
              "oxidoreductase complex", "oxidoreductase activity", "mitochondrial respirasome", "ATP metabolic process")
dotplot(dfj0dseq2_GSE, showCategory=categories, split=c(".sign")) + facet_grid(.~.sign)
dotplot(dfj0dseq2_GSE, showCategory=10, split="ONTOLOGY") + facet_grid(.~.sign)