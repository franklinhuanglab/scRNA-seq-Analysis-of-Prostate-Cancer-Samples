
setwd("../organoid/")
#### V1: 5251
dge=list()
dge_temp$type<-dge_temp$ID
dge_temp$ID<-paste0(as.vector(dge_temp@active.ident),"_",dge_temp$ID)
dge[[1]]=SubsetData(dge_temp,cells = colnames(dge_temp)[dge_temp$type=="BE"])
#dge[[1]]=SubsetData(dge_temp,cells = colnames(dge_temp)[dge_temp$type=="BE"])
#dge[[2]]=SubsetData(dge_temp,cells = colnames(dge_temp)[dge_temp$type=="BE_Henry"])
dge[[2]]=SubsetData(dge_organoid,cells = colnames(dge_organoid)[dge_organoid@active.ident=="BE_organoid"])
for (i in 1:length(dge)) {
dge[[i]] <- NormalizeData(dge[[i]], verbose = FALSE)
dge[[i]] <- FindVariableFeatures(dge[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
dge.anchors <- FindIntegrationAnchors(object.list = dge, dims = 1:50)
dge <- IntegrateData(anchorset = dge.anchors, dims = 1:50)

DefaultAssay(dge) <- 'RNA'
genes.use = Seurat:::CheckFeatures(object.name = dge,data  = dge[['RNA']]@data, features = rownames(dge[['RNA']]@data),verbose = F)
dge =  ScaleData(dge, features = genes.use,vars.to.regress = c('nCount_RNA'))


DefaultAssay(dge) <- "integrated"
# dge$ID<-paste0(dge@active.ident,dge$ID)
# dge$ID[dge$ID=="Club_organoid_PCa"]="Club_organoid"
# dge$ID[dge$ID=="BE_organoidBE_organoid"]="BE_organoid"
dge$type[!(dge$type %in% c("BE","BE_Henry"))]="Organoid"

# Run the standard workflow for visualization and clustering
dge <- ScaleData(dge, verbose = FALSE)
dge <- RunPCA(dge, npcs = 100, verbose = FALSE)
# t-SNE and Clustering
dge <- RunUMAP(dge, reduction = "pca", dims = 1:100)
dge <- FindNeighbors(dge, reduction = "pca", dims = 1:100)
dge <- FindClusters(dge, resolution = 0.4)
# Visualization
dge$type="Tumor_Tissue"
dge$type[dge$ID %in% c("BE_organoid","Tumor_organoid","Club_organoid","hMSC_organoid","Proliferation_organoid","Hillock_organoid")]="Organoid"
p1 <- DimPlot(dge, reduction = "umap", group.by = "type")
p2 <- DimPlot(dge, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

ggsave(file="UMAP_combined_tumor_organoid.pdf",width = 40,height = 20,units = "cm")
DimPlot(dge, reduction = "umap", split.by = "type")
ggsave(file="UMAP_separated_bytumor_organoid.pdf",width = 40,height = 20,units = "cm")



Tab<-table(dge$type,Idents(object = dge))
write.table(Tab,"table_initialclustering.txt",sep="\t",row.names = T,col.names = T)

Tab<-table(dge$ID,Idents(object = dge))
write.table(Tab,"table_ID.txt",sep="\t",row.names = T,col.names = T)

DefaultAssay(dge) <- "RNA"
dge$ID=as.vector(dge$ERGtype)



library(ggpubr)
Featurename<-c("AR_independence","Proliferation","KEGG_cell_cycle","KEGG_DNA_replication","GO_epithelial_proliferation","GO_prostate_gland_development","GO_prostate_gland_growth",
               "Biocarta_cellcycle","Biocarta_MCM","HALLMARK_Epithelial_Mesenchymal_transition","HWANG_PCA","hMSC","LIU_VAV3_PCA_Carcinogenesis","KEGG_cell_cycle","BE","Club","LE","Hillock","ERGpos_Tumor")
Features<-list()
Features[[1]] <- as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/AR_independence.txt",header=FALSE,sep="\n")$V1)
Features[[2]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Proliferation.txt",header=FALSE,sep="\n")$V1)
Features[[3]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/KEGG_CELL_CYCLE.txt",header=FALSE,sep="\n")$V1)
Features[[4]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/KEGG_DNA_replication.txt",header=FALSE,sep="\n")$V1)
Features[[5]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/GO_Epithelial_Cell_Proliferation.txt",header=FALSE,sep="\n")$V1)
Features[[6]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/GO_Prostate_Gland_development.txt",header=FALSE,sep="\n")$V1)
Features[[7]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/GO_Prostate_Gland_Growth.txt",header=FALSE,sep="\n")$V1)
Features[[8]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Biocarta_Cellcycle.txt",header=FALSE,sep="\n")$V1)
Features[[9]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Biocarta_MCM.txt",header=FALSE,sep="\n")$V1)
Features[[10]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt",header=FALSE,sep="\n")$V1)
Features[[11]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/HWANT_Prostate_Cancer.txt",header=FALSE,sep="\n")$V1)
Features[[12]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/hMSC.txt",header=FALSE,sep="\n")$V1)
Features[[13]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/LIU_VAV3_Prostate_Carcinogenesis.txt",header=FALSE,sep="\n")$V1)
Features[[14]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/KEGG_CELL_CYCLE.txt",header=FALSE,sep="\n")$V1)
Features[[15]] <- as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/BE.txt",header=FALSE,sep="\n")$V1)
Features[[16]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Club.txt",header=FALSE,sep="\n")$V1)
Features[[17]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/LE.txt",header=FALSE,sep="\n")$V1)
Features[[18]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Hillock.txt",header=FALSE,sep="\n")$V1)
Features[[19]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_combined/TCGA/ERGpos_Tumor_Markers.txt",header=FALSE,sep="\n")$V1)


#dge$ID<-as.vector(dge@active.ident)
#dge<-SetIdent(dge,value = as.vector(dge$seurat_clusters))
for (i in 1:length(Features)) {
  gene.set <- Features[[i]][1:length(Features[[i]])]
  dge <- AddModuleScore(object = dge,assay = "RNA", features = list(gene.set), ctrl = 5, name = Featurename[i])
  V1<-as.vector(FetchData(object = dge,vars = paste0(Featurename[i],"1")))
  FeaturePlot(dge,features=paste0(c(Featurename[i]),"1"))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))
  ggsave(file=paste0("../Signature_normalcells/Modulescore_cluster_",Featurename[i],"_feature.eps"),width = 20,height = 20,units = "cm")
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() +geom_boxplot(outlier.size = 0.1)+ #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
    stat_compare_means(method = "anova",label.x = 3,label.y = max(V1)+0.05)+
    ggtitle(Featurename[i])
  ggsave(file=paste0("../Signature_normalcells/Modulescore_cluster_",Featurename[i],"_boxplot.eps"),width = 20,height = 20,units = "cm")
  
}


for (i in 1:length(unique(dge@active.ident)))  {
#markers_temp <- FindConservedMarkers(dge, ident.1 = unique(dge@active.ident)[i], grouping.var = "ID")
markers_temp <- FindConservedMarkers(dge, ident.1 = unique(dge@active.ident)[i], grouping.var = "type")
markers_temp$gene<-as.vector(rownames(markers_temp))
top10_1 <- markers_temp %>% top_n(n = 10, wt = PCa_avg_logFC)
top10_2 <- markers_temp %>% top_n(n = 10, wt = Organoid_avg_logFC)
#top10_1 <- markers_temp %>% top_n(n = 10, wt = Normal_avg_logFC)
#top10_2 <- markers_temp %>% top_n(n = 10, wt = PCA_avg_logFC)
genelist1<-unique(c(top10_1$gene,top10_2$gene))

cluster_temp <- subset(dge, idents = unique(dge@active.ident)[i])
cluster_temp<-SetIdent(cluster_temp,value = as.vector(cluster_temp$type))
#cluster_temp<-SetIdent(cluster_temp,value = as.vector(cluster_temp$type))
markers_temp <- FindMarkers(cluster_temp, ident.1 = "PCa", ident.2 = NULL, only.Tumor = TRUE)
markers_temp$gene<-as.vector(rownames(markers_temp))
list_tmp <- markers_temp  %>% top_n(n = 10, wt = abs(avg_logFC))
genelist2<-list_tmp$gene

markers_temp <- FindMarkers(cluster_temp, ident.1 = "Organoid", ident.2 = NULL, only.Tumor = TRUE)
markers_temp$gene<-as.vector(rownames(markers_temp))
list_tmp <- markers_temp  %>% top_n(n = 10, wt = abs(avg_logFC))
genelist3<-list_tmp$gene
  
genelist_temp<-unique(c(genelist1,genelist2,genelist3))
avg.cluster_temp <- log1p(AverageExpression(cluster_temp)$RNA)
avg.cluster_temp$gene <- rownames(avg.cluster_temp)
p1 <- ggplot(avg.cluster_temp, aes(PCa,Organoid)) + geom_point() + ggtitle(paste0("Cluster_",unique(dge@active.ident)[i]))+stat_cor(label.x = 0, label.y = 6.5) +
  stat_regline_equation(label.x = 1, label.y = 6)
LabelPoints(plot = p1, points = genelist_temp, repel = TRUE)

ggsave(file=paste0("Scatter_expression_cluster_",unique(dge@active.ident)[i],".eps"),width = 20,height = 20,units = "cm")
DimPlot(cluster_temp,group.by = "type")
ggsave(file=paste0("UMAP_cluster_",unique(dge@active.ident)[i],".eps"),width = 20,height = 20,units = "cm")

all_markers <- FindAllMarkers(cluster_temp,assay = "RNA", only.Tumor = TRUE, min.pct = 0.001, logfc.threshold = 0.001,test.use = "wilcox")
mito.genes <- grep(pattern = "^MT-", x = rownames(cluster_temp@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(cluster_temp@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)
all_markers<-all_markers[!(all_markers$gene %in% remove_gene),]
top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
dge_temp<-SubsetData(cluster_temp,max.cells.per.ident = 50)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
DoHeatmap(dge_temp,assay = "RNA", features = top10$gene,raster=F) + theme(axis.text.y = element_text(size = 5))
ggsave(file=paste0("Heatmap_bytype_cluster_",unique(dge@active.ident)[i],".eps"),width = 20,height = 20,units = "cm")
}

#####generate signature
dge<-dge_pureERGstatus
dge$ID<-paste0(dge$ERGstatus,"_",dge$ID)
dge<-SetIdent(dge,value = as.vector(dge$ID))
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
for (i in 1:length(unique(dge@active.ident))) {
  markers_temp <- FindMarkers(dge, ident.1 = unique(dge@active.ident)[i], ident.2 = NULL, only.pos = TRUE)
  markers_temp$gene<-as.vector(rownames(markers_temp))
  markers_temp<-markers_temp[order(markers_temp$avg_logFC,decreasing = T),]
  sig_cluster<-markers_temp$gene[markers_temp$p_val_adj<0.05]
  sig_cluster<-sig_cluster[1:min(length(sig_cluster),50)]
  write.table(sig_cluster,paste0("",unique(dge@active.ident)[i],"_Signature.txt"),sep = "\t",quote = F,col.names = F,row.names = F)
}


pbmc.markers <- FindAllMarkers(dge,assay = "RNA",only.pos = TRUE, min.pct = 0.001, logfc.threshold = 0.001,test.use = "wilcox")
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
top100 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_diff)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
dge_temp_2<-SubsetData(dge,max.cells.per.ident = 50)
dge_temp_2 <- ScaleData(object=dge_temp_2,features=rownames(dge_temp_2))
DoHeatmap(dge_temp_2,assay = "RNA",features = top10$gene,raster=F) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_10.pdf",width = 20,height = 10,units = "cm",limitsize = FALSE)
DotPlot(dge, features = unique(top5$gene),assay = "RNA") + RotatedAxis()
ggsave(file="dot_plot_5.eps",width = 20,height = 15,units = "cm",limitsize = FALSE)

write.table(pbmc.markers,"Markers_bycellstate.txt",sep="\t",col.names = T,row.names = T)
dge<-readRDS("./stromal_pca.rds")
dge$ID_2<-as.vector(dge$seurat_clusters)
dge$ID_2="Others"
dge$ID_2[dge$seurat_clusters=="1"]="Cluster_1"
dge<-SetIdent(dge,value = as.vector(dge$ID_2))
dge<-SubsetData(dge,max.cells.per.ident = min(table(dge@active.ident)))
#########Signature score
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/prostate.gmt",header=FALSE,sep="\t")
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_combined/c6_geneset.gmt",header=FALSE,sep="\t")
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/cytokine.gmt",header=FALSE,sep="\t")
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/AR.txt",header=FALSE,sep="\t")
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Inflammation.gmt",header=FALSE,sep="\t")
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/h.all.v7.0.symbols.gmt",header=FALSE,sep="\t")
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/c2.cp.kegg.v7.0.symbols.gmt",header=FALSE,sep="\t")

for (i in 1:nrow(tmp)) {
  name_feature=as.character(tmp[i,1])
  feature_tmp<-na.omit(tmp[i,])
  Features<-as.vector(feature_tmp[3:ncol(feature_tmp)])
  
  dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature,assay = "RNA")
  #names(x = dge[[]])

  V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
  FeaturePlot(dge,features=c(paste0(name_feature,"1")))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(min(V1), max(V1)))
  ggsave(file=paste0("Featureplot_",name_feature,".eps"),width = 15,height = 15,units = "cm")
  
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    #geom_signif(comparisons = list(c("Tumor_AFR","ERGpos_Tumor")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
    geom_signif(comparisons = list(c("BE","BE_Henry")),map_signif_level=TRUE,y_position = max(V1)+0.05,test = "wilcox.test") +
    ggtitle(name_feature)+ RotatedAxis()
  
  ggsave(file=paste0("Violinplot_",name_feature,".eps"),width = 15,height = 15,units = "cm")
}



cluster.averages <- AverageExpression(dge,assays = "integrated",slot = "scale.data",return.seurat = T)
mtx<-cluster.averages$RNA

Feature<-NULL
Feature$names=c("KLK3","KLK2","ACPP","MSMB","AR","LSAMP","KRT5","KRT15","KRT17","DST","TP63","ERG","GOLM1","EEF2","STEAP4","PCA3","AMACR","SPON2","PCAT4","PIGR","MMP7","LTF","CP","IGFBP3","GLUL","KLF5","SLC14A1","DDIT4")

Type<-Feature$type[order(Feature$type)]
Names<-Feature$names[order(Feature$type)]
Feature$names=as.vector(Names)
Feature$type=as.vector(Type)
Feature<-as.data.frame(Feature)

mtx<-cluster.averages$RNA
mtx<-mtx[rownames(mtx) %in% Feature$names,]
mtx<-log2(mtx+1)
mtx<-mtx[Feature$names,]


mtx<-na.omit(mtx[order(match(rownames(mtx),Feature$names)),])
mtx<-mtx[,c(1,11,12,13,14,15,10,3,5,7,9,2,4,6,8)]
#for club, use this:
mtx<-mtx[,c(10,11,12,13,14,9,2,4,6,8,1,3,5,7)]
mtx_norm<-na.omit(t(apply(mtx, 2, function(x)(x-min(x))/(max(x)-min(x)))))
metadata<-data.frame(as.character(Feature$type),row.names = Feature$names)
colnames(metadata)<-"Marker_Type"
require(pheatmap)
pheatmap(mtx_norm,cluster_rows = F,cluster_cols = F,filename = "E_supervisedheatmap.pdf",width = 10,height = 10)

######## Heatmap of cell type signature score

# pheatmap(mtx_norm,cluster_rows = F,cluster_cols = F)
# ,filename = "E_supervisedheatmap.pdf",width = 10,height = 10)
FeaturePlot(Club_PCA,features=c("AR"))
ggsave(file="E_AR.eps",width = 20,height = 20,units = "cm")
VlnPlot(dge,features=c("AR"),assay = "RNA",group.by = "ID_bytype")
ggsave(file="E_AR_vlnplot.eps",width = 40,height = 20,units = "cm")

EGFR_pathway<-c("PIK3R1","PIK3CA","PRKCA","SOS1","MAPK3","ELK1","SHC1","HRAS","EGF","CSNK2A1",
                "GRB2","JAK1","JUN","PLCG1","PRKCB","MAPK8","MAP2K1","RAF1","RASA1","MAP2K4",
                "SRF","STAT3","STAT5A","EGFR","FOS","STAT1")
EGFR_pathway<-as.vector(unlist(read.csv("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/geneset/BIOCARTA_EGF_PATHWAY.txt",header=FALSE,sep="\t")))
AR_pathway<-as.vector(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/AR.txt",header=FALSE,sep="\t")))
EMT_pathway<-as.vector(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt",header=FALSE,sep="\t")))
KEGG_hormone_pathway<-as.vector(unlist(read.csv("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/geneset/KEGG_STEROID_HORMONE_BIOSYNTHESIS.txt")))
BIOCARTA_TGFB_PATHWAY<-as.vector(unlist(read.csv("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/geneset/BIOCARTA_TGFB_PATHWAY.txt")))
WNT_pathway<-as.vector(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/KEGG_WNT_SIGNALING_PATHWAY.txt",header=FALSE,sep="\t")))
Biocarta_RAS_pathway<-as.vector(unlist(read.csv("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/geneset/BIOCARTA_RAS_PATHWAY.txt",header=FALSE,sep="\t")))
Hallmark_RAS<-as.vector(unlist(read.csv("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/geneset/HALLMARK_KRAS_SIGNALING_UP.txt",header=FALSE,sep="\t")))
Hallmark_AKT<-as.vector(unlist(read.csv("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/geneset/HALLMARK_PI3K_AKT_MTOR_SIGNALING.txt",header=FALSE,sep="\t")))

ASC_pathway<-Feature_ASC

dge_temp<-readRDS("../E_integrated.rds")
DefaultAssay(dge_temp)="RNA"
Features=c("BIOCARTA_TGFB_PATHWAY1","EMT_signature1","TGFB1","TGFB2","TGFBR1","TGFBR2","TGFBR3","TGFB3","HSD17B3")
#dge_temp <- AddModuleScore(object = dge_temp, features = list(AR_pathway), name ="AR_signature",assay = "RNA")
dge_temp<-AddModuleScore(object = dge_temp,features=list(EMT_pathway),name="EMT_pathway",assay = "RNA")
dge_temp<-AddModuleScore(object = dge_temp,features=list(EGFR_pathway),name="EGFR_pathway",assay = "RNA")
dge_temp<-AddModuleScore(object = dge_temp,features=list(WNT_pathway),name="WNT_pathway",assay = "RNA")
dge_temp<-AddModuleScore(object = dge_temp,features=list(Biocarta_RAS_pathway),name="Biocarta_RAS",assay = "RNA")
dge_temp<-AddModuleScore(object = dge_temp,features=list(Hallmark_RAS),name="Hallmark_RAS",assay = "RNA")
dge_temp<-AddModuleScore(object = dge_temp,features=list(CDC24_feature),name="CDC24_pathway",assay = "RNA")

VlnPlot(dge_temp,features=c("WNT_pathway1"),group.by = "ID")
VlnPlot(dge_temp,features=c("EMT_pathway1"),group.by = "ID")
VlnPlot(dge_temp,features=c("Biocarta_RAS1"),group.by = "ID")
VlnPlot(dge_temp,features=c("Hallmark_RAS1"),group.by = "ID")
VlnPlot(dge_temp,features=c("EGFR_pathway1"),group.by = "ID")
VlnPlot(dge_temp,features=c("KRAS.PROSTATE_UP.V1_UP1"),group.by = "ID")

targeted_cluster<-c("6")
dge_temp$target<-"Others"
dge_temp$target[dge_temp@active.ident %in% targeted_cluster]<-"Target"

Features=c("WNT_pathway1","EMT_pathway1","Biocarta_RAS1","Hallmark_RAS1","EGFR_pathway1","CDC24_pathway1")
Features=c("ZBTB16","STEAP4","FKBP5")
for (i in 1:length(Features)){
V1<-as.vector(FetchData(object = dge_temp,vars = Features[i]))
#V1<-as.vector(FetchData(object = dge_temp,vars = "KEGG_hormone_pathway1"))
#FeaturePlot(dge_temp,features="AR")#+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, 0.6))
#FeaturePlot(dge_temp,features="KEGG_hormone_pathway1")#+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, 0.6))

#ggsave(file="AR_expression_E_integrated.eps",width = 20,height = 20,units = "cm")
#ggsave(file="KEGG_hormone_pathway_E_integratde.eps",width = 20,height = 20,units = "cm")
V1<-as.vector(FetchData(object = dge_temp,vars = Features[i]))
BOX_df<-NULL
BOX_df$id<-dge_temp$target
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
  #geom_signif(comparisons = list(c("BE","BE_Normal")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  #geom_signif(comparisons = list(c("Club","Club_Normal")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  #geom_signif(comparisons = list(c("LE","LE_Normal")),map_signif_level=TRUE,y_position = max(V1)+0.25)
  geom_signif(comparisons = list(c("Target","Others")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  ggtitle(Features[i])+ RotatedAxis()

#ggsave(file="KEGG_hormone_pathway_vlnplot.pdf",width = 15,height = 10,units = "cm")
ggsave(file=paste0(Features[i],"_vlnplot_PCa.pdf"),width = 30,height = 20,units = "cm")
}

dge<-club_PCA
dge$target<-"Others"
dge$target[dge@active.ident=="0"]="Target"
dge <- AddModuleScore(object = dge, features = list(KEGG_hormone_pathway), name ="KEGG_hormone_pathway",assay = "RNA")
V1<-as.vector(FetchData(object = dge,vars = "AR_pathway1"))
FeaturePlot(dge,features="KEGG_hormone_pathway1")
setwd("../../../Seqwell_combined/combined_analysis/")
ggsave(file="PCa_E_KEGG_hormone_featureplot.eps",width = 20,height = 20,units = "cm")
BOX_df<-NULL
BOX_df$id<-as.vector(dge@active.ident)
BOX_df$id<-as.vector(dge$target)
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() +
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
  geom_signif(comparisons = list(c("Target","Others")),map_signif_level=TRUE,y_position = max(V1)+0.25,test = "wilcox.test",test.args = "greater") +
  ggtitle("AR expression")+ RotatedAxis()
ggsave(file="BE_6_AR_expression_integrated.eps",width = 20,height = 20,units = "cm")
mean(BOX_df$value[BOX_df$id=="Target"])
mean(BOX_df$value[BOX_df$id=="Others"])

#########Heatmap of signature score for each cell state
tmp<-list()
for (i in 1:6) {
  #tmp[[i]]<-read.csv(paste0("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/BE_signature/BE_",i-1,"_Signature.txt"),header=FALSE,sep="\n")
  #tmp[[i]]<-read.csv(paste0("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/LE_signature/LE_",i-1,"_Signature.txt"),header=FALSE,sep="\n")
  tmp[[i]]<-read.csv(paste0("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/Club_signature/Club_",i-1,"_Signature.txt"),header=FALSE,sep="\n")
  #tmp[[i]]<-read.csv(paste0("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/Tumor_signature/",i-1,"_Signature.txt"),header=FALSE,sep="\n")
}


for (i in 1:length(tmp)) {
  name_feature=paste0("Cluster_",i-1,"Club")
  feature_tmp<-unlist(tmp[[i]])
  Features<-as.vector(feature_tmp)
  dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature,assay = "RNA")
  
  #dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature,assay = "RNA")
  #names(x = dge[[]])
  
  V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
  FeaturePlot(dge,features=c(paste0(name_feature,"1")))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, 1.5))
  ggsave(file=paste0("Featureplot_club_",name_feature,"_CPRC_club.eps"),width = 20,height = 20,units = "cm")
} 
###Club
Featurelist<-c("KLK3","KLK2","LTF","LCN2","SCGB3A1","PSCA","MYC","KRT23","KRT7","CCND1","LAMC2","KRT5","KRT13","SNCG","GPX2","CRYAB","HSPA8","HSPD1","HSPH1","SERPINH1","FKBP4","FHL2","CXCL2","ATF3","FOS","JUN","PIGR","MMP7","CP","AR","NKX3-1")
CellstateID<-c(rep(0,3),rep(1,3),rep(2,6),rep(3,3),rep(4,7),rep(5,4),rep(6,5))
###BE
Featurelist<-c("CCL2","MT2A","KRT14","TPM2","DKK1","ATF3","EGFR","JUN","AZGP1","MMP7","MT1X","WFDC2","HSPH1","HSPA8","HSPD1","PLAU","FOSL1","CCL20","G0S2","KRT13","SCGB1A1","CXCL17","SERPINB1","ACPP","KLK3","KLK2","FOS","KIT","KRT5","SPINK2","IL1RL1","IL32")
CellstateID<-c(rep(0,3),rep(1,3),rep(2,6),rep(3,3),rep(4,7),rep(5,4),rep(6,5))
###LE
Featurelist<-c("ORM1","ORM2","FOLH1","VEGFA","MME","GOLM1","TGFB3","IGFBP5","NKX3-1","NIPAL3","TMEFF2","NCOA4","NEFH",
               "SPON2","S100A6","CXCL8","KRT19","KRT17","TFF1","JUN","ATF3","FOS","EGR1","FOSB",
               "RHOB","KLF6","LTF","TGM2","CD74","HLA-DRA","PIGR","LCN2","SFTPA2","GDF15","PDLIM5",
               "CST1","MT1X","MT2A","NEAT1","KRT8","AR")
Featurelist<-CDC24_feature
Featurelist<-c("IL17D","IL17RA","IL17RB","IL17RC","IL17RD","IL17RE",
               "IL17A","IL17B","IL17C","IL17F")
Featurelist<-c("CD247","ENTPD1","GZMA","GZMB","GZMK")
Featurelist<-c("TGFB1","TGFB2","TGFBR1","TGFBR2")
for (i in 1:length(Featurelist)) {
  name_feature=Featurelist[i]
  #FeaturePlot(LE_Henry,features=Featurelist[i])+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, 5))
  #ggsave(file=paste0(name_feature,"_Henry_LE_Cluster.eps"),width = 20,height = 20,units = "cm")
  V1<-as.vector(FetchData(object = dge_temp,vars = name_feature))
  BOX_df<-NULL
  #BOX_df$id<-dge_temp$target
  BOX_df$id<-as.vector(dge_temp@active.ident)
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    geom_signif(comparisons = list(c("ERGneg_Stromal","ERGpos_Stromal")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
    geom_signif(comparisons = list(c("ERGneg_CD4","ERGpos_CD4")),map_signif_level=TRUE,y_position = max(V1)+0.5) +
    geom_signif(comparisons = list(c("ERGneg_Tumor","ERGpos_Tumor")),map_signif_level=TRUE,y_position = max(V1)+0.75)
  #geom_signif(comparisons = list(c("Target","Others")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  ggtitle(name_feature)+ RotatedAxis()
  ggsave(file=paste0(name_feature,".pdf"),width = 30,height = 20,units = "cm")
} 


dge_Tcell$ERGtype<-"ERGneg"
dge_Tcell$ERGtype[dge_Tcell$orig.ident %in% c("AUG_PB1","MAY_PB1","PR5249_T","PR5249_N","PR5251_T","PR5251_N","PR5199","PR5269")]="ERGpos"

dge_temp<-SubsetData(dge_Tcell,ident.use = c("CD8","CD8_cluster3","CD8_cluster2"))
Featurelist<-c("CD247","ENTPD1","GZMA","GZMB","GZMK")
Featurelist<-c("TIGIT","HAVCR2","CTLA4","LAG3","ITGA4")
Featurelist<-c("TGFB1","TGFB2","CXCR4","IFNG","IFNGR1","IFNGR2","TGFBR1","TGFBR2")
Featurelist<-c("PDCD1","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1","HLA-DRB5","STAT1","MT2A","HLA-A","HLA-B","HLA-C",
               "HLA-E","HLA-F","HLA-G","ICAM1","CD44","B2M","GBP2","GBP4","GBP5")




cells_pure<-colnames(dge_ERGstatus)[dge_ERGstatus$orig.ident %in% c("MAY_PB2","PR5186","PR5196","PR5199","PR5269","PR5254_T","PR5254_N","PR5261_T","PR5261_N")]
dge_pureERGstatus<-SubsetData(dge_ERGstatus,cells = cells_pure)

dge_temp<-SubsetData(dge_pureERGstatus,ident.use = c("ERGneg_Stromal","ERGpos_Stromal"))
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp@active.ident))
source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
dge_temp<-SubsetData(dge_tumor,cells = colnames(dge_tumor)[dge_tumor$ID=="ERGneg_Tumor"])
Data<-dge_temp
name1<-"Ident"
name2<-"ERGnegTumor_ARstatus"
Prepare.for.GSEA(Data,name1,name2)


FC_df<-NULL
for (i in 1:length(Featurelist)) {
  name_feature=Featurelist[i]
  #FeaturePlot(LE_Henry,features=Featurelist[i])+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, 5))
  #ggsave(file=paste0(name_feature,"_Henry_LE_Cluster.eps"),width = 20,height = 20,units = "cm")
  V1<-as.vector(FetchData(object = dge_temp,vars = name_feature))
  BOX_df<-NULL
  #BOX_df$id<-dge_temp$target
  BOX_df$id<-as.vector(dge_temp$ERGtype)
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  Pos_mean=mean(BOX_df$value[BOX_df$id=="Pos"])
  Neg_mean=mean(BOX_df$value[BOX_df$id=="Neg"])
  FC=log(Neg_mean/Pos_mean)/log(2)
  FC_temp<-c(Pos_mean,Neg_mean,FC)
  FC_df<-rbind(FC_df,FC_temp)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    geom_signif(comparisons = list(c("Neg","Pos")),map_signif_level=TRUE,y_position = max(V1)+0.25,test = "t.test")
  #geom_signif(comparisons = list(c("Target","Others")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  ggtitle(name_feature)+ RotatedAxis()
  ggsave(file=paste0(name_feature,"_in_CD4.pdf"),width = 20,height = 15,units = "cm")
} 

FC_df<-as.data.frame(FC_df)
colnames(FC_df)=c("Pos_mean","Neg_mean","log2FC")
rownames(FC_df)=as.vector(Featurelist)
write.table(FC_df,"CD4_DEG_FC.txt",sep="\t",quote = F,col.names = T,row.names = T)

FeaturePlot(LE_PCA,features="AR")+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, 3))
ggsave(file="AR_PCA_LE_Cluster.eps",width = 20,height = 20,units = "cm")

write.table(Club_PCA@reductions$umap@cell.embeddings,"../Club_integrated/Club_PCA_UMAP.txt",sep="\t",col.names = T,row.names = T)




#####volcano plot
dge_temp<-club_PCA
dge_temp<-SubsetData(dge_E,ident.use = c("ERGpos_Tumor","LE"))
targeted_cluster<-c("0","5")
dge_temp$target<-"Others"
dge_temp$target[dge_temp$seurat_clusters %in% targeted_cluster]<-"Target"
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$target))
all_markers <-FindAllMarkers(dge_temp,min.diff.pct=0.001,min.pct =  0.001,test.use = "wilcox",assay = "RNA")
mito.genes <- grep(pattern = "^MT-", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
all_markers<-all_markers[!(all_markers$gene %in% remove_gene),]
write.table(all_markers,"Club_0_Others_markers.txt",sep="\t",row.names=F, col.names = T)
all_markers<-all_markers[1:(nrow(all_markers)/2),]
DE_select<-all_markers$p_val_adj<0.05
length(which(DE_select))
all_markers$threshold <- DE_select
ggplot(all_markers) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj), colour=threshold)) +
  ggtitle(paste0(targeted_cluster,"_vs_rest overexpression")) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
## Sort by ordered padj
all_markers_reorder <- all_markers[order(abs(all_markers$avg_logFC),decreasing = T), ] 

## Create a column to indicate which genes to label
all_markers_reorder$genelabels <- ""
all_markers_reorder$genelabels[1:20] <- rownames(all_markers_reorder)[1:20]
all_markers_reorder$genelabels[all_markers_reorder$gene =="KLK2"]="KLK2"
all_markers_reorder$genelabels[all_markers_reorder$gene =="KLK3"]="KLK3"
all_markers_reorder$genelabels[all_markers_reorder$gene =="ACPP"]="ACPP"
all_markers_reorder$genelabels[all_markers_reorder$gene =="MSMB"]="MSMB"
all_markers_reorder$genelabels[all_markers_reorder$gene =="NKX3-1"]="NKX3-1"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="SCGB1A1"]="SCGB1A1"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="CD74"]="CD74"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="LTF"]="LTF"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="KRT8"]="KRT8"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="KRT18"]="KRT18"
ggplot(all_markers_reorder) +
  geom_point(aes(x = avg_logFC, y = -log10(p_val_adj), colour = threshold)) +
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = ifelse(genelabels !="", genelabels,""))) +
  ggtitle(paste0(targeted_cluster,"_vs_rest overexpression")) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave(file="LE_12_vs_rest_volcano.eps",width = 20,height = 10,units = "cm")


VlnPlot(dge_temp,features=c("SCGB.genes1"))

dge_infer_mtx <- as.matrix(GetAssayData(dge_organoid,assay = "RNA", slot = "counts"))
write.table(dge_infer_mtx,"LE_counts.txt",sep="\t",row.names = T,col.names = T,quote=F)
write.table(dge@meta.data,"LE_metadata.txt",sep="\t",row.names = T,col.names = T)
write.table(dge@reductions$umap@cell.embeddings,"LE_UMAP.txt",sep="\t",row.names = T,col.names = T)


IL6_pathway<-c("FOS","JUN","CEBPB","IL6ST","IL6R","MAP2K1","HRAS","MAPK3","RAF1","SHC1","ELK1","PTPN11","SRF","GRB2","SOS1","JAK2","IL6")
BCR_pathway<-c("FOS","JUN","CD79A","MAP2K1","PPP3CB","HRAS","MAPK3","RAF1","SHC1","ELK1","VAV1","MAPK8","MAP3K1","PPP3CC","NFATC1","GRB2","PRKCB","PRKCA","SOS1","PLCG1","CD79B","SYK","CALM3","LYN")

club_0_sig<-read.table("./Club_signature/Club_0_Signature.txt",sep="\n")

dge<-readRDS("./Tumor_integrated.rds")
setwd("../Tumor/")
mkdir("Immune")
setwd("./Immune/")
targeted_pathway<-BCR_pathway
DimPlot(dge,group.by = "ID")
DefaultAssay(dge)="integrated"
dge <- ScaleData(dge, verbose = FALSE)
dge <- RunPCA(dge,features =c(t(targeted_pathway)),  verbose = T)
DimPlot(dge,reduction = "pca")
print(dge[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(dge,reduction = "pca")

DefaultAssay(dge)<-"RNA"
FeatureScatter(dge, "FOS", "JUN")
pcs.plot=paste("PC",1:2,sep="")
FeaturePlot(dge,pcs.plot, pt.size = 2)
FeaturePlot(dge, features = c("FOS", "JUN"), blend = TRUE)
FeaturePlot(dge, features = c("KLK3", "AR"), blend = TRUE)
FeaturePlot(dge, features = c("FOS", "IL6ST"), blend = TRUE)
FeaturePlot(dge, features = c("LTF", "BIRC3"), blend = TRUE)
ggsave(file="./LTF_BIRC3_coexpression.eps",width = 20,height = 20,units = "cm")
RidgePlot(object = dge, feature = c("MSMB","MGP","TIMP1"))
ggsave(file="./Ridgeplot_markers_0.eps",width = 20,height = 20,units = "cm")

targeted_pathway<-read.table("./Club_signature/Club_3_Signature.txt",sep="\n")
Filename1=c("BIOCARTA_CXCR4_PATHWAY","BIOCARTA_EGF_PATHWAY","BIOCARTA_IGF1_PATHWAY","BIOCARTA_IL2_PATHWAY",
           "BIOCARTA_IL6_PATHWAY","BIOCARTA_INSULIN_PATHWAY","KEGG_ANTIGEN","KEGG_BCELL_SIGNALING","PID_IL6_PATHWAY",
           "PID_TCR_CALCIUM_PATHWAY","REACTOME_ANTIGEN","REACTOME_INTERFEREN_GAMMA","REACTOME_NOTCH2",
           "REACTOME_PD1","REACTOME_TCR","REACTOME_WNT_CANCER")
Filename2=c("REACTOME_RRNA_PROCESSING","REACTOME_METABOLISM_RNA","REACTOME_SIGNALING_BY_ROBO",
           "BIOCARTA_P53_PATHWAY","REACTOME_TRANSLATION","KEGG_RIBOSOME","REACTOME_INTERFERON",
           "REACTOME_IMMUNE_CYTOKINE","REACTOME_CELLULAR_HEAT_STRESS","PID_IL1_PATHWAY",
           "REACTOME_ATTENUATION_PHASE","BIOCARTA_DEATH_PATHWAY","KEGG_MAPK_PATHWAY",
           "REACTOME_EGFR_SIGNALING","BIOCARTA_IL1R_PATHWAY","KEGG_CYTOKINE_INTERACTION",
           "REACTOME_CHEMOKINE_RECEPTORS","PID_IL4_2PATHWAY","KEGG_NKCELL_CYTOTOXICITY",
           "BIOCARTA_CTLA4_PATHWAY","PID_KIT_PATHWAY","PID_TCR_PATHWAY","BIOCARTA_NKT_PATHWAY")
Filename3=c("PID_BMP_PATHWAY","EGF_SIGNALING_UP","TOMLINS_PROSTATE_CANCER_UP","SETLUR_TMPRSS2_ERG_FUSION","LIU_PROSTATE_CANCER")

Filename4<-c("SECRETORY_PATHWAY","GO_INTERFERON_ALPHA_SECRETION","GO_INTERFERON_BETA_SECRETION",
             "GO_INTERFERON_GAMMA_SECRETION","GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY",
             "GO_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY","HALLMARK_PROTEIN_SECRETION")

Filename_common<-set_immune_stromal_tumor

Filename<-c(Filename1,Filename2,Filename3)

#dge_temp<-SubsetData(dge,cells = colnames(dge)[dge$ID=="ERGneg_Tumor"])
dge_temp<-dge
DefaultAssay(dge_temp)="RNA"
dge_temp<-readRDS("../LE_integrated.rds")
dge_temp<-dge_ERGneg
dge_temp$ID="Others"
dge_temp$ID[dge_temp@active.ident %in% c("2","1")]="Cluster 1&2"
Filename<-as.vector(set_immune_stromal_tumor)
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$ID))
for (i in 1:length(Filename)) {
targeted_pathway<-read.table(paste0("~/Desktop/Seqwell_combined/combined_analysis/geneset/",Filename[i],".txt"),sep="\n")
targeted_pathway<-read.table("~/Desktop/Seqwell_combined/combined_analysis/LE_Signature.txt",sep="\n")
pathway_name<-as.vector(unlist(targeted_pathway)[1])
targeted_pathway<-as.vector(unlist(targeted_pathway)[3:length(unlist(targeted_pathway))])
#dge_temp<-SubsetData(dge,ident.use = "3")


dge_temp <- AddModuleScore(object = dge_temp,assay = "RNA", features = list(targeted_pathway), ctrl = 5, name = pathway_name)

V1<-as.vector(FetchData(object = dge_temp,vars = paste0(pathway_name,"1")))
FeaturePlot(dge_temp,features=paste0(pathway_name,"1"))#+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(min(V1), max(V1)))
ggsave(file=paste0("./",pathway_name,"_inCBE_feature.eps"),width = 20,height = 20,units = "cm")
BOX_df<-NULL
BOX_df$id<-dge_temp@active.ident
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
  geom_signif(comparisons = list(c("4","6")),map_signif_level=TRUE,y_position = max(V1)+0.15)+
  ggtitle(pathway_name)
ggsave(file=paste0(pathway_name,"inBE_boxplot.eps"),width = 20,height = 20,units = "cm")


# #DefaultAssay(dge_temp)="RNA"
# matrix<-dge_temp@assays$RNA@data
# matrix_mod<-as.matrix(matrix)
# gene<-matrix_mod[rownames(matrix_mod) %in% unlist(targeted_pathway),]
# #rcorr(t(gene), type = c("pearson"))
# cormat <- round(cor(t(gene)),2)
# cormat[is.na(cormat)] <- 0
# 
# reorder_cormat <- function(cormat){
#   # Use correlation between variables as distance
#   dd <- as.dist((1-cormat)/2)
#   hc <- hclust(dd)
#   cormat <-cormat[hc$order, hc$order]
}
cormat <- reorder_cormat(cormat)
library(reshape2)
melted_cormat <- melt(cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()
ggsave(file=paste0(pathway_name,"_correlation_inERGneg.eps"),width = 20,height = 20,units = "cm")


#####Visualize GSEA results
install.packages("readxl")
library(readxl)
library(forcats)

dot_df_1<-read.table("./GSEA_Club_0_Hallmark_UP.txt",sep = "\t",header = T)
dot_df_2<-read.table("./GSEA_Club_0_Hallmark_UP.txt",sep = "\t",header = T)

dot_df_1<-read.table("./ERGneg_ERGpos_Tumor_up.txt",sep = "\t",header = T)
dot_df_2<-read.table("./ERGneg_ERGpos_Tumor_down.txt",sep = "\t",header = T)

dot_df<-rbind(dot_df_1[1:20,],dot_df_2[1:20,])

dot_df<-rbind(dot_df_1,dot_df_2)
dot_df<-dot_df[dot_df$FDR.q.val<0.1,]
dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"
dot_df_ERGneg_ERGpos_Tumor<-dot_df
set_tumor<-as.vector(dot_df$NAME)



dot_df_1<-read.table("./Mix_ERGpos_CD4_up.txt",sep = "\t",header = T)
dot_df_2<-read.table("./Mix_ERGpos_CD4_down.txt",sep = "\t",header = T)
dot_df_3<-read.table("./Mix_ERGpos_CD8_up.txt",sep = "\t",header = T)
dot_df_4<-read.table("./Mix_ERGpos_CD8_down.txt",sep = "\t",header = T)
dot_df_5<-read.table("./Mix_ERGpos_Stromal_up.txt",sep = "\t",header = T)
dot_df_6<-read.table("./Mix_ERGpos_Stromal_down.txt",sep = "\t",header = T)

dot_df<-rbind(dot_df_1,dot_df_2)
dot_df<-dot_df[dot_df$FDR.q.val<0.1,]
dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"
dot_df_Mix_ERGpos_CD4<-dot_df[dot_df$FDR.q.val<0.1,]
dot_df<-rbind(dot_df_3,dot_df_4)
dot_df<-dot_df[dot_df$FDR.q.val<0.1,]
dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"
dot_df_Mix_ERGpos_CD8<-dot_df[dot_df$FDR.q.val<0.1,]
dot_df<-rbind(dot_df_5,dot_df_6)
dot_df<-dot_df[dot_df$FDR.q.val<0.1,]
dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"
dot_df_Mix_ERGpos_Stromal<-dot_df[dot_df$FDR.q.val<0.1,]


# dot_df<-dot_df_tumor[dot_df_tumor$NAME %in% set_immune_stromal_tumor,]
# dot_df<-dot_df_Stromal[dot_df_Stromal$NAME %in% set_immune_stromal_tumor,]


#set_Tumor_up<-as.vector(dot_df_Mix_ERGpos_Tumor$NAME[dot_df_Mix_ERGpos_Tumor$type=="upregulated"])
#set_Tumor_down<-as.vector(dot_df_Mix_ERGpos_Tumor$NAME[dot_df_Mix_ERGpos_Tumor$type=="downregulated"])
set_CD4_up<-as.vector(dot_df_Mix_ERGpos_CD4$NAME[dot_df_Mix_ERGpos_CD4$type=="upregulated"])
set_CD4_down<-as.vector(dot_df_Mix_ERGpos_CD4$NAME[dot_df_Mix_ERGpos_CD4$type=="downregulated"])
set_CD8_up<-dot_df_Mix_ERGpos_CD8$NAME[dot_df_Mix_ERGpos_CD8$type=="upregulated"]
set_CD8_up<-as.vector(set_CD8_up)
set_CD8_down<-dot_df_Mix_ERGpos_CD8$NAME[dot_df_Mix_ERGpos_CD8$type=="downregulated"]
set_CD8_down<-as.vector(set_CD8_down)
set_Stromal_up<-as.vector(dot_df_Mix_ERGpos_Stromal$NAME[dot_df_Mix_ERGpos_Stromal$type=="upregulated"])
set_Stromal_down<-as.vector(dot_df_Mix_ERGpos_Stromal$NAME[dot_df_Mix_ERGpos_Stromal$type=="downregulated"])

x=list(Tumor=set_Tumor_down, CD4=set_CD4_down, CD8=set_CD8_down,Stromal=set_Stromal_down)
x=list(Tumor=set_Tumor_up, CD4=set_CD4_up, CD8=set_CD8_up,Stromal=set_Stromal_up)
#MSET(x, 2232, FALSE)
res=supertest(x, n=2232)
pdf("Intersection_Mix_ERGpos_Tumor_Tcell_Stromal_DOWN.pdf")
plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.15))
dev.off()
pdf("Intersection_Mix_ERGpos_Tumor_Tcell_Stromal_DOWN_barchart.pdf")
plot(res, Layout="landscape", degree=1:4, sort.by="size", margin=c(0.5,5,1,2))
dev.off()

###Normal vs PCA: p<1e-10
###ERGpos vs ERGneg: Tumor, BE, LE, Club: p = 0.0005

set_immune_stromal_tumor<-intersect(intersect(set_Tumor_up,set_CD8_up),set_Stromal_up)

set_immune_stromal_tumor<-intersect(intersect(set_Tumor_up,set_CD4_up),set_Stromal_up)

dot_df<-dot_df_CD4_ERGneg[dot_df_CD4_ERGneg$NAME %in% set_immune_stromal_tumor,]

dot_df<-dot_df_ERGneg_Mix_CD4[dot_df_ERGneg_Mix_CD4$NAME %in% set_immune_stromal_tumor,]

dot_df<-dot_df_ERGneg_ERGpos_CD8[dot_df_ERGneg_ERGpos_CD8$NAME %in% set_immune_stromal_tumor,]

info_column<-as.vector(dot_df$LEADING.EDGE)
DF_nes<-strsplit(info_column, ", ")

dot_df$tags=0
for (i in 1:length(DF_nes)) {
  dot_df$tags[i]<-as.numeric(as.vector(gsub("%", "", as.vector(gsub("tags=", "", DF_nes[[i]][1])))))
}
dot_df$GeneRatio<-dot_df$SIZE*dot_df$tags/100
dot_df$Counts<-dot_df$SIZE*dot_df$tags/100
dot_df$GeneRatio<-dot_df$tags/100

p <- ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(NAME, GeneRatio))) + 
  geom_point(aes(size = Counts, color = FDR.q.val)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.25), low="red") +
  ylab(NULL) +
  ggtitle("C2CP GSEA")

p + facet_grid(.~type)
ggsave(file=paste0("./Upregulated_ERGneg_Mix_Tumor_CD4_CD8_Stromal.pdf"),width = 35,height = 15,units = "cm")


set_stromal_tumor<-intersect(dot_df_Stromal$NAME,dot_df_tumor$NAME)
set_CD4_tumor<-intersect(dot_df_CD4$NAME,dot_df_tumor$NAME)
set_CD8_tumor<-intersect(dot_df_CD8$NAME,dot_df_tumor$NAME)
set_Tcell_tumor<-intersect(dot_df_alltcell$NAME,dot_df_tumor$NAME)
set_immune_stromal_tumor<-intersect(set_stromal_tumor,set_CD4_tumor)


#AR_pathway,EGFR_pathway
SCGB.genes <- grep(pattern = "^SCGB", x = rownames(club_PCA@assays$RNA@data), value = TRUE)
pathway_name="EGFR_pathway"
#[1] "SCGB1D2"  "SCGB2A1"  "SCGB1A1"  "SCGB3A1"  "SCGB2A2"  "SCGB2B2"  "SCGB3A2"  "SCGB1B2P"
club_integrated <- AddModuleScore(object = club_integrated,assay = "RNA", features = list(EGFR_pathway), ctrl = 5, name = pathway_name)

dge_temp<-dge_E

DefaultAssay(dge_temp)="RNA"
targeted_cluster<-c("6")
dge_temp$target<-"Others"
dge_temp$target[dge_temp$seurat_clusters %in% targeted_cluster]<-"Target"
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$target))

V1<-as.vector(FetchData(object = dge_temp,vars = paste0(pathway_name,"1")))
V1<-as.vector(FetchData(object = dge_temp,vars = "KRT8"))
FeaturePlot(dge_temp,features="AR")+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))
BOX_df<-NULL
BOX_df$id<-BE_PCA@active.ident
BOX_df$id<-dge_temp$target
BOX_df$id<-dge_temp$ID
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
  #stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
  geom_signif(comparisons = list(c("LE","Club")),map_signif_level=TRUE,y_position = max(V1)+0.15)+ 
  #geom_signif(comparisons = list(c("ERGneg_Tumor","Club")),map_signif_level=TRUE,y_position = max(V1)+0.35)+ 
  #geom_signif(comparisons = list(c("ERGpos_Tumor","Club")),map_signif_level=TRUE,y_position = max(V1)+0.55)+ 
  #geom_signif(comparisons = list(c("BE","Club")),map_signif_level=TRUE,y_position = max(V1)+0.15)+ 
  ggtitle("KRT8_expression")
ggsave(file=paste0("../KRT8","_boxplot_E.pdf"),width = 10,height = 10,units = "cm")

dge<-readRDS("./ERGneg/ERG_LE.rds")
 DimPlot(dge)
 unique(dge$ID)
 unique(dge$orig.ident)
 dge$orig.ident[dge$orig.ident %in% c("PR5249_T","PR5249_N")]="PR5249"
 dge$orig.ident[dge$orig.ident %in% c("PR5251_T","PR5251_N")]="PR5251"
 dge$orig.ident[dge$orig.ident %in% c("PR5254_T","PR5254_N")]="PR5254"
 dge$orig.ident[dge$orig.ident %in% c("PR5261_T","PR5261_N")]="PR5261"
 dge<-SetIdent(dge,value = as.vector(dge$orig.ident))
 DimPlot(dge)
 dge$ID_new<-as.vector(dge$orig.ident)
 dge$ID_new[dge$ID=="LE"]="LE"
 DimPlot(dge,group.by = "ID_new")
 dge<-SetIdent(dge,value = as.vector(dge$ID_new))
 DimPlot(dge)
ggsave(file="./ERGneg_LE_patient.eps",width = 20,height = 20,units = "cm")
FeaturePlot(dge,features="AR")+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))


VlnPlot(dge_temp,features=c("AHNAK","SCUBE2","EGR1","CP","IGFBP3"),group.by = "ID")
ggsave("DEG_group2.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("KLF6","ATF3","PIGR","MMP7","OLFM4"),group.by = "ID")
ggsave("DEG_group3.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("RNA28S5","PABPC1","GNAS","DSP","CTSB","HLA-B"),group.by = "ID")
ggsave("DEG_group1.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("FTH1","S100A6","EEF1A1","RPL10","RPLP1"),group.by = "ID")
ggsave("DEG_group4.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("RPL13","RPL13A","TMSB4X","RPL15","RPL18A"),group.by = "ID")
ggsave("DEG_group5.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("FTL","MALAT1","NEAT1","SAT1","LCN2"),group.by = "ID")
ggsave("DEG_group6.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("HSPB1","SCGB3A1","B2M"),group.by = "ID")
ggsave("DEG_group7.eps",width = 30,height = 30,units = "cm")


VlnPlot(dge_temp,features=c("FOS","JUN","TXNIP","ACPP","PABPC1"),group.by = "ID")
ggsave("DEG_group1.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("CTSB","DSP","GNAS","KLK3","AHNAK"),group.by = "ID")
ggsave("DEG_group2.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("PSAP","SLC14A1","EGR1","ATF3","FOSB"),group.by = "ID")
ggsave("DEG_group3.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("BTG2","MSMB","KRT15","OLFM4","RNA28S5"),group.by = "ID")
ggsave("DEG_group4.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("MALAT1","RPLP1","RPL13","RPL10","KRT17"),group.by = "ID")
ggsave("DEG_group5.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("CXCL8","EEF1A1","RPL13A","RPL7","RPL34"),group.by = "ID")
ggsave("DEG_group6.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("FTH1","S100A6","ACTB","S100A2","RPL32","RPL3"),group.by = "ID")
ggsave("DEG_group7.eps",width = 30,height = 30,units = "cm")
VlnPlot(dge_temp,features=c("FTL","RPL18A","RPL15","HSPB1"),group.by = "ID")
ggsave("DEG_group8.eps",width = 30,height = 30,units = "cm")



dge_temp<-BE_PCA
DefaultAssay(dge_temp)="RNA"
targeted_cluster<-c("6")
dge_temp$target<-"Others"
dge_temp$target[dge_temp$seurat_clusters %in% targeted_cluster]<-"Target"
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$target))

Featurelist<-intersect(Nelson_pathway,as.vector(rownames(dge_temp)))
Featurelist<-c(intersect(Hallmark_AKT,as.vector(rownames(dge_temp))),"TNK2")
Featurelist<-AR_pathway
FC_df<-NULL
for (i in 1:length(Featurelist)) {
  name_feature=Featurelist[i]
  V1<-as.vector(FetchData(object = dge_temp,vars = name_feature))
  BOX_df<-NULL
  BOX_df$id<-dge_temp$target
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)

  Target_mean=mean(BOX_df$value[BOX_df$id=="Target"])
  Others_mean=mean(BOX_df$value[BOX_df$id=="Others"])
  FC=log(Target_mean/Others_mean)/log(2)
  FC_temp<-c(Target_mean,Others_mean,FC)
  FC_df<-rbind(FC_df,FC_temp)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    #geom_signif(comparisons = list(c("ERGneg","ERGpos")),map_signif_level=TRUE,y_position = max(V1)+0.25)
    geom_signif(comparisons = list(c("Target","Others")),map_signif_level=TRUE,y_position = max(V1)+0.25)
  #geom_signif(comparisons = list(c("Target","Others")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  ggtitle(name_feature)+ RotatedAxis()
  ggsave(file=paste0(name_feature,"_in_BE_PCA.pdf"),width = 20,height = 15,units = "cm")
} 
FC_df<-as.data.frame(FC_df)
#colnames(FC_df)=c("Pos_mean","Neg_mean","log2FC")
colnames(FC_df)=c("Target_mean","Others_mean","log2FC")
rownames(FC_df)=as.vector(Featurelist)
write.table(FC_df,"Club_PCa_Nelson_AR_markers.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(FC_df,"Club_PCa_AR_markers.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(FC_df,"BE_PCa_AKT_markers.txt",sep="\t",quote = F,col.names = T,row.names = T)

dge<-club_PCA
min=-0.2754
max=0.3598
Nelson_pathway<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/geneset/NELSON_RESPONSE_TO_ANDROGEN_UP.txt",header=FALSE,sep="\n")$V1)
dge_temp <- AddModuleScore(object = dge_temp,assay = "RNA", features = list(Nelson_pathway), ctrl = 5, name = "Nelson_AR")
dge_temp$AR_pathway1<-(dge_temp$Hallmark_AR-mean(dge_temp$Hallmark_AR))/sd(dge_temp$Hallmark_AR)
dge_temp$AR_pathway1<-(dge_temp$AR_pathway1-min(dge_temp$AR_pathway1))/(max(dge_temp$AR_pathway1)-min(dge_temp$AR_pathway1))*(0.5+0.2754)-0.2754
ggsave("Club_0_AR_Nelson_signature_inPCa.eps")
VlnPlot(dge_temp,features=c("AR_pathway1"),group.by = "target")
VlnPlot(dge_temp,features=c("Nelson_AR1"),group.by = "target")
ggsave("BE_6_AR_Nelson_signature_inPCa.eps")


-0.2899403 
V1<-as.vector(FetchData(object = TCGA,vars = "",slot = "data"))
V1<-V1[,1]
cor.test(V1,as.numeric(TCGA$Gleasonscore))
cor.test(V1,as.numeric(TCGA$OS_month))



cells_ARpos_tumor<-colnames(SubsetData(dge_tumor,subset.name = "AR",low.threshold = 0))
dge_tumor$ARstatus<-"Neg"
dge_tumor$ARstatus[colnames(dge_tumor) %in% cells_ARpos_tumor]="Pos"
DimPlot(dge_tumor,group.by = "ARstatus")
dge_tumor<-SetIdent(dge_tumor,value = as.vector(dge_tumor$ARstatus))

all_markers <-FindAllMarkers(dge_temp,min.diff.pct=0.001,min.pct =  0.001,test.use = "wilcox",assay = "RNA")
mito.genes <- grep(pattern = "^MT-", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
all_markers<-all_markers[!(all_markers$gene %in% remove_gene),]
all_markers_AR_ERGneg<-all_markers

VlnPlot(dge_temp,features=c("Luminal1","BE1","NE1","Club1"),pt.size = 0)
V1<-as.vector(FetchData(object = dge_temp,vars = "ERBB2"))
BOX_df<-NULL
BOX_df$id<-as.vector(dge_temp@active.ident)
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)

ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
  #geom_signif(comparisons = list(c("ERGneg","ERGpos")),map_signif_level=TRUE,y_position = max(V1)+0.25)
  geom_signif(comparisons = list(c("Pos","Neg")),map_signif_level=TRUE,y_position = max(V1)+0.25)+
  ggtitle("BE_signature")+ RotatedAxis()


Filename=c("REACTOME_SIGNALING_BY_FGFR","REACTOME_SIGNALING_BY_FGFR_IN_DISEASE",
           "BIOCARTA_MAPK_PATHWAY","KEGG_MAPK_SIGNALING_PATHWAY",
           "REACTOME_ONCOGENIC_MAPK_SIGNALING","BIOCARTA_AKT_PATHWAY",
           "BIOCARTA_MTOR_PATHWAY","PID_PI3KCI_AKT_PATHWAY")


for (i in 1:length(Filename)) {
  targeted_pathway<-read.table(paste0("~/Desktop/Seqwell_combined/combined_analysis/geneset/",Filename[i],".txt"),sep="\n")
  pathway_name<-Filename[i]
  targeted_pathway<-as.vector(unlist(targeted_pathway)[2:length(unlist(targeted_pathway))])
  dge_temp <- AddModuleScore(object = dge_temp,assay = "RNA", features = list(targeted_pathway), ctrl = 5, name = pathway_name)
  #dge_temp<-SubsetData(dge,ident.use = "3")
  #df=rbind(intersect(targeted_pathway,rownames(dge_temp)))
  V1<-as.vector(FetchData(object = dge_temp,vars = paste0(pathway_name,"1")))
  BOX_df<-NULL
  BOX_df$id<-as.vector(dge_temp@active.ident)
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    #geom_signif(comparisons = list(c("ERGneg","ERGpos")),map_signif_level=TRUE,y_position = max(V1)+0.25)
    geom_signif(comparisons = list(c("Pos","Neg")),map_signif_level=TRUE,y_position = max(V1)+0.25)+
    ggtitle(pathway_name)+ RotatedAxis()
  ggsave(paste0(pathway_name,"_ERGneg.pdf"))
}








