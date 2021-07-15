
setwd("../organoid/")
#### V1: 5251
dge=list()
#dge_temp$type<-dge_temp$ID
#dge_temp$ID<-paste0(as.vector(dge_temp@active.ident),"_",dge_temp$ID)
dge[[1]]=SubsetData(dge_E,cells = colnames(dge_E)[dge_E$orig.ident %in% 
                                                    c("PR5249_N", "PR5249_T", "PR5251_N", "PR5251_T", "PR5254_N", "PR5254_T", "PR5261_N","PR5261_T","PR5269")])
#dge[[1]]=SubsetData(dge_temp,cells = colnames(dge_temp)[dge_temp$type=="BE"])
#dge[[2]]=SubsetData(dge_temp,cells = colnames(dge_temp)[dge_temp$type=="BE_Henry"])
dge[[2]]=SubsetData(dge_organoid,cells = colnames(dge_organoid)[dge_organoid$orig.ident %in% 
                                                                  c("HNW_PR5249_N_org_pool1_pool2_S42_L003","HNW_PR5251_N_org_pool1_S35_L002",
                                                                    "HNW_PR5251_N_org_pool2_S41_L003","HNW_PR5251_T_org_pool1_pool2_S40_L003",
                                                                    "HNW_PR5254_N_org_pool1_S39_L003","HNW_PR5254_T_org_S38_L003",
                                                                    "PR5249_N_org_pool3_S57_L004",
                                                                    "PR5251_T_org_pool3_S62_L004","PR5254_N_org_pool2_pool3_S65_L004",
                                                                    "PR5261_N_V1_org_S69_L004","PR5261_T_V1_org_S55_L004",
                                                                    "PR5269_T_org_S60_L004","HNW_PR5261_N_V2_org_S37_L003",
                                                                    "PR5261_T_V2_org_S58_L004")])

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
dge <- FindVariableFeatures(object = dge)
# dge$ID<-paste0(dge@active.ident,dge$ID)
# dge$ID[dge$ID=="Club_organoid_PCa"]="Club_organoid"
# dge$ID[dge$ID=="BE_organoidBE_organoid"]="BE_organoid"
dge$type[!(dge$type %in% c("BE","BE_Henry"))]="Organoid"

# Run the standard workflow for visualization and clustering
dge <- ScaleData(dge, verbose = FALSE)
dge <- RunPCA(dge, npcs = 100, verbose = FALSE)
# t-SNE and Clustering
dge <- RunUMAP(dge, reduction = "pca", dims = 1:100)
dge <- FindNeighbors(dge, reduction = "pca", dims = 1:50)
dge <- FindClusters(dge, resolution = 0.4)
# Visualization
dge$type[dge$ID=="BE_organoid"]="Organoid"
p1 <- DimPlot(dge, reduction = "umap", group.by = "type")
p2 <- DimPlot(dge, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DefaultAssay(dge) <- "integrated"
dge <- ScaleData(dge, verbose = FALSE)
dge <- RunPCA(dge, npcs = 100, verbose = FALSE)
dge <- RunUMAP(dge, reduction = "pca", dims = 1:100)
dge <- FindNeighbors(dge, reduction = "pca", dims = 1:100)
dge <- FindClusters(dge, resolution =c(0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,2,2.5))

dge <- FindClusters(dge, resolution =1.0)

pdf("PR5316_stability.pdf",width = 10,height = 15)
clustree(dge, prefix = "integrated_snn_res.")
clustree(dge, prefix = "RNA_snn_res.")
dev.off()


ggsave(file="UMAP_combined.eps",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "umap", split.by = "type")
ggsave(file="UMAP_separated_bytype.eps",width = 40,height = 20,units = "cm")



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

BE_markers<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_combined/TCGA/BE_Markers.txt",header=FALSE,sep="\n")$V1)
Club_markers<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Club.txt",header=FALSE,sep="\n")$V1)
dge<-SubsetData(dge_overlap,ident.use = c("BE","BE_organoid","Club","Club_organoid"))
dge<-SetIdent(dge,value = as.vector(dge@active.ident))
dge <- AddModuleScore(object = dge,assay = "RNA", features = list(BE_markers), ctrl = 5, name = "BE_signature")
dge <- AddModuleScore(object = dge,assay = "RNA", features = list(Club_markers), ctrl = 5, name = "Club_signature")

#dge$ID<-as.vector(dge@active.ident)
#dge<-SetIdent(dge,value = as.vector(dge$seurat_clusters))
for (i in 1:length(Features)) {
  gene.set <- Features[[i]][1:length(Features[[i]])]
  dge <- AddModuleScore(object = dge,assay = "RNA", features = list(gene.set), ctrl = 5, name = Featurename[i])
  V1<-as.vector(FetchData(object = dge,vars = paste0(Featurename[i],"1")))
  FeaturePlot(dge,features=paste0(c(Featurename[i]),"1"))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))
  ggsave(file=paste0("../Signature_normalcells/Modulescore_cluster_",Featurename[i],"_feature.eps"),width = 20,height = 20,units = "cm")
  V1<-as.vector(FetchData(object = dge,vars ="Club_signature1"))
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
    #stat_compare_means(method = "anova",label.x = 3,label.y = max(V1)+0.05)+
    geom_signif(comparisons = list(c("Club","Club_organoid")),map_signif_level=TRUE,y_position = max(V1)+0.45,test.args = "less") +
    ggtitle(Featurename[i])
  ggsave(file="Club_signaturescore_overlap.eps",width = 20,height = 20,units = "cm")
  #ggsave(file=paste0("../Signature_normalcells/Modulescore_cluster_",Featurename[i],"_boxplot.eps"),width = 20,height = 20,units = "cm")
  
}


for (i in 6:length(unique(dge@active.ident)))  {
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

all_markers <- FindAllMarkers(cluster_temp,assay = "RNA", only.Tumor = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
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
for (i in 1:length(unique(dge@active.ident))) {
  markers_temp <- FindMarkers(dge, ident.1 = unique(dge@active.ident)[i], ident.2 = NULL, only.pos = TRUE)
  markers_temp$gene<-as.vector(rownames(markers_temp))
  markers_temp<-markers_temp[order(markers_temp$avg_logFC,decreasing = T),]
  sig_cluster<-markers_temp$gene[markers_temp$p_val_adj<0.05]
  sig_cluster<-sig_cluster[1:min(length(sig_cluster),100)]
  write.table(sig_cluster,paste0("",unique(dge@active.ident)[i],"_Signature.txt"),sep = "\t",quote = F,col.names = F,row.names = F)
}

DefaultAssay(dge)="RNA"
pbmc.markers <- FindAllMarkers(dge,assay = "RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
top100 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_diff)
top100 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top15 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
##0,1,2,4,5
dge_temp<-SubsetData(dge,ident.use = "2")
#dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
marker_temp<-top10$gene[top10$cluster=="2"]
marker_temp_back<-top15$gene[top15$cluster=="2"]
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$type))
#marker_temp[7]="SLC14A1"
#DotPlot(dge_temp, features = marker_temp,assay = "RNA",scale = F) + RotatedAxis()+coord_flip()
DotPlot(dge_temp, features = marker_temp[1:10], split.by = "type",scale = F,cols = c("blue", "red")) + RotatedAxis()+coord_flip()
ggsave(file="BE_cluster2_dotplot.eps",width = 10,height = 15,units = "cm",limitsize = FALSE)


dge_temp<-SubsetData(dge,max.cells.per.ident = 100)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))

levels(x = dge_temp) <- c("BE","BE_organoid","Club","Club_organoid","ERGpos_Tumor","ERGneg_Tumor","Tumor_organoid")

DoHeatmap(dge_temp,assay = "RNA",features = top10$gene,raster=F) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_10_PR5269.pdf",width = 30,height = 20,units = "cm",limitsize = FALSE)
DotPlot(dge, features = unique(top5$gene),assay = "RNA") + RotatedAxis()
DotPlot(dge, features = c("KRT5","KRT15","KRT17","TP63","DST","PIGR","LTF","MMP7","LCN2","SCGB1A1"),assay = "RNA") + RotatedAxis()
ggsave(file="dot_plot_5_overlap.pdf",width = 40,height = 15,units = "cm",limitsize = FALSE)

write.table(top100,"Markers_bycelltype_overlap.txt",sep="\t",col.names = T,row.names = T)
dge<-readRDS("./stromal_pca.rds")
dge$ID_2<-as.vector(dge$seurat_clusters)
dge$ID_2="Others"
dge$ID_2[dge$seurat_clusters=="1"]="Cluster_1"
dge<-SetIdent(dge,value = as.vector(dge$ID_2))
dge<-SubsetData(dge,max.cells.per.ident = min(table(dge@active.ident)))

dge$sample<-as.vector(dge$orig.ident)
dge$sample[dge$sample %in% c("HNW_PR5249_N_org_pool1_pool2_S42_L003","PR5249_N_org_pool3_S57_L004")]="PR5249_N_organoid"
dge$sample[dge$sample %in% c("HNW_PR5251_N_org_pool1_S35_L002","HNW_PR5251_N_org_pool2_S41_L003")]="PR5251_N_organoid"
dge$sample[dge$sample %in% c("HNW_PR5251_T_org_pool1_pool2_S40_L003","PR5251_T_org_pool3_S62_L004")]="PR5251_T_organoid"
dge$sample[dge$sample %in% c("HNW_PR5254_N_org_pool1_S39_L003","PR5254_N_org_pool2_pool3_S65_L004")]="PR5254_N_organoid"
dge$sample[dge$sample %in% c("HNW_PR5254_T_org_S38_L003")]="PR5254_T_organoid"
dge$sample[dge$sample %in% c("PR5261_N_V1_org_S69_L004","HNW_PR5261_N_V2_org_S37_L003")]="PR5261_N_organoid"
dge$sample[dge$sample %in% c("PR5261_T_V1_org_S55_L004","PR5261_T_V2_org_S58_L004")]="PR5261_T_organoid"
dge$sample[dge$sample %in% c("PR5269_T_org_S60_L004")]="PR5269_organoid"

dge<-SetIdent(dge,value = as.vector(dge$sample))
levels(x = dge) <- c("PR5249_T","PR5249_N","PR5249_N_organoid",
                     "PR5251_T","PR5251_T_organoid","PR5251_N","PR5251_N_organoid",
                     "PR5254_T","PR5254_T_organoid","PR5254_N","PR5254_N_organoid",
                     "PR5261_T","PR5261_T_organoid","PR5261_N","PR5261_N_organoid",
                     "PR5269","PR5269_organoid")

dge_PR5261<-SubsetData(dge,ident.use = c("PR5261_T","PR5261_N","PR5261_N_organoid","PR5261_T_organoid"))
dge_PR5261<-SetIdent(dge_PR5261,value = as.vector(dge_PR5261$ID))
dge_PR5261<-SubsetData(dge_PR5261,ident.remove = c("ERGpos_Tumor","ERGneg_Tumor","Tumor_organoid"))
dge_PR5261<-SetIdent(dge_PR5261,value = as.vector(dge_PR5261$ID))

dge<-SubsetData(dge_overlap,ident.use = c("BE","BE_organoid","Club","Club_organoid"))
dge<-SetIdent(dge,value = as.vector(dge@active.ident))
pbmc.markers <- FindAllMarkers(dge,assay = "RNA",only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.01,test.use = "wilcox")
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
top10_dge <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top20_dge <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top5_dge <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
dge_temp<-SubsetData(dge,max.cells.per.ident = 100)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))

av.exp <- AverageExpression(dge)$RNA
av.exp<-dge@assays$RNA@scale.data
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, c(as.vector(colnames(dge))))
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0,limit = c(-1,1))
ggsave(file="./PCa_organoid_correlation_PR5261.pdf",width = 20,height = 20,units = "cm")


#levels(x = dge_temp) <- c("BE","BE_organoid","Club","Club_organoid","ERGpos_Tumor","ERGneg_Tumor","Tumor_organoid")
levels(x = dge_temp) <- c("BE","BE_organoid","Club","Club_organoid")

DoHeatmap(dge_temp,assay = "RNA",features = top20_dge$gene,raster=F) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_10_dge_BE_Club.pdf",width = 30,height = 20,units = "cm",limitsize = FALSE)
DotPlot(dge, features = unique(top5_dge$gene),assay = "RNA") + RotatedAxis()
ggsave(file="dot_plot_5_dge.pdf",width = 30,height = 15,units = "cm",limitsize = FALSE)
DotPlot(dge, features = c("CAV1","KRT5","KRT15","KRT17","TP63","DST","PIGR","LTF","MMP7","LCN2","SCGB1A1","CP","PSCA"),assay = "RNA") + RotatedAxis()
ggsave(file="dot_plot_5_knownmarkers_dge.pdf",width = 20,height = 15,units = "cm",limitsize = FALSE)

#########Signature score
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/prostate.gmt",header=FALSE,sep="\t")
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_combined/c6_geneset.gmt",header=FALSE,sep="\t")
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/cytokine.gmt",header=FALSE,sep="\t")

for (i in 1:nrow(tmp)) {
  name_feature=as.character(tmp[i,1])
  feature_tmp<-na.omit(tmp[i,])
  Features<-as.vector(feature_tmp[3:ncol(feature_tmp)])
  
  dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature,assay = "RNA")
  #names(x = dge[[]])

  V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
  FeaturePlot(dge,features=c(paste0(name_feature,"1")))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))
  ggsave(file=paste0("Featureplot_",name_feature,".png"),width = 20,height = 20,units = "cm")
  
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    #geom_signif(comparisons = list(c("Tumor_AFR","ERGpos_Tumor")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
    #geom_signif(comparisons = list(c("ERGpos","ERGneg")),map_signif_level=TRUE,y_position = max(V1)+0.45) +
    ggtitle(name_feature)+ RotatedAxis()
  
  ggsave(file=paste0("Violinplot_",name_feature,".png"),width = 20,height = 20,units = "cm")
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
tmpname<-c("Proliferation","Cellcycle","DNA_replication","EMT","hMSC","BE","Club","LE","Hillock")
tmp<-list()
tmp[[1]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Proliferation.txt",header=FALSE,sep="\n")$V1)
tmp[[2]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/KEGG_CELL_CYCLE.txt",header=FALSE,sep="\n")$V1)
tmp[[3]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/KEGG_DNA_replication.txt",header=FALSE,sep="\n")$V1)
tmp[[4]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt",header=FALSE,sep="\n")$V1)
tmp[[5]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/hMSC.txt",header=FALSE,sep="\n")$V1)
tmp[[6]] <- as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/BE.txt",header=FALSE,sep="\n")$V1)
tmp[[7]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Club.txt",header=FALSE,sep="\n")$V1)
tmp[[8]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/LE.txt",header=FALSE,sep="\n")$V1)
tmp[[9]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Hillock.txt",header=FALSE,sep="\n")$V1)
Mean_score_celltype<-NULL
for (i in 1:length(tmp)) {
  name_feature=tmpname[i]
  feature_tmp<-unlist(tmp[[i]])
  Features<-as.vector(feature_tmp)
  
  dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature,assay = "RNA")
  #names(x = dge[[]])
  
  V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
  BOX_df<-NULL
  BOX_df$id<-as.vector(dge@active.ident)
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  df_temp<-aggregate(BOX_df$value, list(BOX_df$id), mean)
  Mean_score_celltype<-rbind(Mean_score_celltype,as.numeric(t(df_temp)[2,]))
}
colnames(Mean_score_celltype)<-as.vector(t(df_temp)[1,])
#colnames(Mean_score)<-dge$ERGtype
rownames(Mean_score_celltype)<-tmpname
pheatmap(Mean_score_celltype,cluster_rows = T,cluster_cols = T,filename = "Organoid_bycelltype_heatmap.pdf",width = 10,height = 10)

#########Heatmap of signature score for each cell state
tmp<-list()
for (i in 1:6) {
  #tmp[[i]]<-read.csv(paste0("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/BE_signature/BE_",i-1,"_Signature.txt"),header=FALSE,sep="\n")
  #tmp[[i]]<-read.csv(paste0("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/LE_signature/LE_",i-1,"_Signature.txt"),header=FALSE,sep="\n")
  tmp[[i]]<-read.csv(paste0("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/Club_signature/Club_",i-1,"_Signature.txt"),header=FALSE,sep="\n")
  #tmp[[i]]<-read.csv(paste0("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/Tumor_signature/",i-1,"_Signature.txt"),header=FALSE,sep="\n")
}

Mean_score_BE<-NULL
for (i in 1:length(tmp)) {
  name_feature=paste0("Cluster_",i-1,"Tumor")
  feature_tmp<-unlist(tmp[[i]])
  Features<-as.vector(feature_tmp)
  TCGA <- AddModuleScore(object = TCGA, features = list(Features), name = name_feature,assay = "RNA")
  
  #dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature,assay = "RNA")
  #names(x = dge[[]])

  V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
  FeaturePlot(dge,features=c(paste0(name_feature,"1")))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))
  ggsave(file=paste0("Featureplot_",name_feature,"_organoid_BE.eps"),width = 20,height = 20,units = "cm")
  
  BOX_df<-NULL
  #BOX_df$id<-dge$orig.ident
  #BOX_df$id<-colnames(dge)
  #BOX_df$id<-dge$Gleason
  BOX_df$id<-as.vector(dge@active.ident)
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  df_temp<-aggregate(BOX_df$value, list(BOX_df$id), mean)
  Mean_score_BE<-rbind(Mean_score_BE,as.numeric(t(df_temp)[2,]))

  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    #geom_signif(comparisons = list(c("Tumor_AFR","ERGpos_Tumor")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
    #geom_signif(comparisons = list(c("ERGpos","ERGneg")),map_signif_level=TRUE,y_position = max(V1)+0.45) +
    ggtitle(name_feature)+ RotatedAxis()
  ggsave(file=paste0("Violinplot_",name_feature,"_BE.eps"),width = 30,height = 20,units = "cm")
}
colnames(Mean_score_BE)<-paste0(as.vector(t(df_temp)[1,]),"_organoid")
#colnames(Mean_score)<-dge$ERGtype
rownames(Mean_score_BE)<-paste0(0:(nrow(Mean_score_BE)-1),"_PCa")
pheatmap(Mean_score_BE,cluster_rows = T,cluster_cols = T,filename = "Organoid_heatmap_Clubcluster.pdf",width = 10,height = 10)


###for club
Mean_score<-Mean_score[,c(1,2,3,4,5,14,7,9,11,13,6,8,10,12)]
require(pheatmap)
pheatmap(Mean_score_BE,cluster_rows = T,cluster_cols = T,filename = "Organoid_heatmap_Clubcluster.pdf",width = 10,height = 10)
mtx_norm<-Mean_score
mtx_norm<-na.omit(apply(Meanscore, 2, function(x)(x-min(x))/(max(x)-min(x))))
pheatmap(mtx_norm,cluster_rows = F,cluster_cols = F,filename = "Tumor_cluster_heatmap_TCGA_bygleason.pdf",width = 10,height = 10)

# pheatmap(mtx_norm,cluster_rows = F,cluster_cols = F)
# ,filename = "E_supervisedheatmap.pdf",width = 10,height = 10)
FeaturePlot(dge,features=c("AR"))
ggsave(file="E_AR.eps",width = 20,height = 20,units = "cm")
VlnPlot(dge,features=c("AR"),assay = "RNA",group.by = "ID_bytype")
ggsave(file="E_AR_vlnplot.eps",width = 40,height = 20,units = "cm")

EGFR_pathway<-c("PIK3R1","PIK3CA","PRKCA","SOS1","MAPK3","ELK1","SHC1","HRAS","EGF","CSNK2A1",
                "GRB2","JAK1","JUN","PLCG1","PRKCB","MAPK8","MAP2K1","RAF1","RASA1","MAP2K4",
                "SRF","STAT3","STAT5A","EGFR","FOS","STAT1")
dge <- AddModuleScore(object = dge, features = list(EGFR_pathway), name ="EGFR_pathway",assay = "RNA")
V1<-as.vector(FetchData(object = dge,vars = "EGFR_pathway1"))
FeaturePlot(dge,features="EGFR_pathway1")+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))
ggsave(file="E_EGFR_Featureplot.eps",width = 20,height = 20,units = "cm")
BOX_df<-NULL
BOX_df$id<-dge$ID_bytype
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
  #geom_signif(comparisons = list(c("Tumor_AFR","ERGpos_Tumor")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  #geom_signif(comparisons = list(c("6","4")),map_signif_level=TRUE,y_position = max(V1)+0.45) +
  ggtitle("EGFR Signature Score")+ RotatedAxis()
ggsave(file="EGFR_Violinplot.eps",width = 20,height = 20,units = "cm")

V1<-as.vector(FetchData(object = dge,vars = "ERG"))
FeaturePlot(dge_E,features="ERG")+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(min(V1), max(V1)))
ggsave(file="ERG_Featureplot.eps",width = 20,height = 20,units = "cm")

V1<-as.vector(FetchData(object = dge,vars = "ETV1"))
FeaturePlot(dge_E,features="ETV1")+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(min(V1), max(V1)))
ggsave(file="ETV1_Featureplot.eps",width = 20,height = 20,units = "cm")

V1<-as.vector(FetchData(object = dge,vars = "ETV5"))
FeaturePlot(dge_E,features="ETV5")+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(min(V1), max(V1)))
ggsave(file="ETV5_Featureplot.eps",width = 20,height = 20,units = "cm")

V1<-as.vector(FetchData(object = dge,vars = "ERG"))
FeaturePlot(dge_E,features="ERG")+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(min(V1), max(V1)))
ggsave(file="ERG_Featureplot.eps",width = 20,height = 20,units = "cm")


ERGtype<-read.table("../../TCGA/ERG_fusion_status.tsv",sep="\t",header = T)
dge<-SubsetData(TCGA,cells = colnames(dge)[colnames(dge) %in% ERGtype$Sample.ID])
dge$ERGtype<-"None"
for (i in 1:ncol(dge)) {
  dge$ERGtype[i]=as.character(ERGtype$ERG..FUSION[ERGtype$Sample.ID==colnames(dge)[i]])
}

paste(as.character(text1[1,]),sep = "\t")

rownames(Mean_score_BE)<-paste0("BE_",rownames(Mean_score_BE))
rownames(Mean_score_LE)<-paste0("LE_",rownames(Mean_score_LE))
rownames(Mean_score_Club)<-paste0("Club_",rownames(Mean_score_Club))
Meanscore<-rbind(Mean_score_BE,Mean_score_LE,Mean_score_Club)
mtx_norm<-na.omit(apply(Meanscore, 2, function(x)(x-min(x))/(max(x)-min(x))))

pheatmap(Meanscore,cluster_rows = T,cluster_cols = T,filename = "E_cluster_heatmap_TCGA_raw.pdf",width = 10,height = 10)

Mean_score_BE_new<-(Mean_score_BE-min(Mean_score_BE))/(max(Mean_score_BE)-min(Mean_score_BE))
pheatmap(Mean_score_BE_new,cluster_rows = F,cluster_cols = F,filename = "BE_cluster_heatmap_TCGA_zscore.pdf",width = 10,height = 10)

Meanscore<-(Meanscore-mean(Meanscore))/sd(Meanscore)

my_sample_col <- data.frame(sample = as.vector(colnames(mtx_norm)))
row.names(my_sample_col) <- make.unique(as.vector(colnames(mtx_norm)))
colnames(my_sample_col)="Gleason_score"
pheatmap(mtx_norm,annotation_col = my_sample_col,cluster_rows = T,cluster_cols = T,filename = "E_cluster_heatmap_TCGA_normalized.pdf",width = 10,height = 10)

write.table(Meanscore,"TCGA_Signaturescore.txt",sep="\t",col.names = T,row.names = T)

dge_temp<-dge
dge_temp$Cluster_6BE1[!(dge_temp@active.ident=="BE_organoid")]=-0.2
V1<-as.vector(FetchData(object = dge_temp,vars = "Cluster_6BE1"))
FeaturePlot(dge_temp,features=c("Cluster_6BE1"))+scale_color_gradientn( colours = c("grey","green","yellow","red"),  limits = c(min(V1), max(V1)))
ggsave(file="BE_cluster6.eps",width = 20,height = 20,units = "cm")

dge_temp<-dge
dge_temp$Cluster_4BE1[!(dge_temp@active.ident=="BE_organoid")]=-0.2
V1<-as.vector(FetchData(object = dge_temp,vars = "Cluster_4BE1"))
FeaturePlot(dge_temp,features=c("Cluster_4BE1"))+scale_color_gradientn( colours = c("grey","green","yellow","red"),  limits = c(min(V1), max(V1)))
ggsave(file="BE_cluster4.eps",width = 20,height = 20,units = "cm")

dge_temp<-dge
dge_temp$Cluster_0Club1[!(dge_temp@active.ident=="Club_organoid")]=min(dge_temp$Cluster_0Club1)
V1<-as.vector(FetchData(object = dge_temp,vars = "Cluster_0Club1"))
FeaturePlot(dge_temp,features=c("Cluster_0Club1"))+scale_color_gradientn( colours = c("grey","green","yellow","red"),  limits = c(min(V1), max(V1)))
ggsave(file="Club_cluster0.eps",width = 20,height = 20,units = "cm")

dge_temp<-dge
dge_temp$Cluster_3Club1[!(dge_temp@active.ident=="Club_organoid")]=min(dge_temp$Cluster_3Club1)
V1<-as.vector(FetchData(object = dge_temp,vars = "Cluster_3Club1"))
FeaturePlot(dge_temp,features=c("Cluster_3Club1"))+scale_color_gradientn( colours = c("grey","green","yellow","red"),  limits = c(min(V1), max(V1)))
ggsave(file="Club_cluster3.eps",width = 20,height = 20,units = "cm")

dge<-readRDS("./Club_integrated.rds")
dge<-SetIdent(dge,value = as.vector(dge$type))
pbmc.markers <- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.001, logfc.threshold = 0.001)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge_1_pca@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)
pbmc.markers<-pbmc.markers[!(rownames(pbmc.markers) %in% remove_gene),]
#dge_subset is 12 idents, 50 cells each.
top50_Club_Organoid <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top50_Club_Tumor <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top20_Tumor <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

genelist<-top20$gene
genelist[1]="SGK1"
dge_temp<-SubsetData(dge,max.cells.per.ident = 100)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
heatmap_gene<-top50_Club_Pair$gene[top50_Club_Pair$p_val_adj<=0.1]
DoHeatmap(dge_temp, features = c(genelist),raster = F) + theme(axis.text.y = element_text(size = 8))#+ scale_fill_gradientn(colors = c("white", "red"))
ggsave(file="./Heatmap_Club_PCA_Organoid.pdf",width = 20,height = 20,units = "cm")
write.table(top50_Club_Organoid,"Club_organoid_top50_DEG.txt",col.names = T,row.names = T)

levels(x = dge_temp) <- c("ERGneg_Tumor","ERGpos_Tumor","ERGneg_Stromal","ERGpos_Stromal","ERGneg_CD4","ERGpos_CD4")
levels(x = dge_temp)<-c("BE","BE_Normal","Club","Club_Normal","LE","LE_Normal","Hillock_Normal","ERGpos_Tumor","ERGneg_Tumor")
DoHeatmap(dge_temp, features = c(top20_CD4$gene,top20_Stromal$gene,top20_Tumor$gene),raster = F) + theme(axis.text.y = element_text(size = 5))#+ scale_fill_gradientn(colors = c("white", "red"))
DoHeatmap(dge_temp, features = c(top50_Club_Tumor$gene),raster = F) + theme(axis.text.y = element_text(size = 6))#+ scale_fill_gradientn(colors = c("white", "red"))
DoHeatmap(dge_temp, features = unique(c(top50_BE$gene,top50_Club$gene,top50_LE$gene)),raster = F) + theme(axis.text.y = element_text(size = 6))#+ scale_fill_gradientn(colors = c("white", "red"))

ggsave(file="./Heatmap_Club_Tumor_purple.pdf",width = 20,height = 30,units = "cm")
ggsave(file="./Heatmap_ERGstatus_purple_byERGtype_top20.eps",width = 20,height = 25,units = "cm")
write.table(top50_Club_Tumor,"Club_Tumor_top50_DEG.txt",col.names = T,row.names = T)

dge <- ScaleData(object=dge,features=rownames(dge))
DotPlot(dge, features = c("SCGB1A1","PIGR","MMP7","CP","CD74","KRT5","KRT15","KRT17","DST","TP63"),assay = "RNA") + RotatedAxis()
ggsave(file="dot_plot_5.eps",width = 20,height = 15,units = "cm",limitsize = FALSE)

cluster.averages_organoid <- AverageExpression(dge,assays = "RNA",return.seurat = T)
genelist_6<-genelist
genelist_5<-genelist_6[1:50]
genelist_temp<-Response_features
DoHeatmap(cluster.averages_organoid,assay = "RNA",slot = "data", features = unlist(genelist_temp),raster = F, size = 3,draw.lines = FALSE,lines.width = 0.5)+ 
  scale_fill_gradientn(colors = c("blue", "red"))+ theme(axis.text.y = element_text(size = 4))
DoHeatmap(cluster.averages_organoid,assay = "RNA",features = unlist(c("IGFBP3","NLRP1","GSTP1","GLUL","KRT7","MUC4","LAMB3","CHGB","DDIT4","LY6D","DLPI","SDC1")),raster = F, size = 3,draw.lines = FALSE,lines.width = 0.5)+ 
  scale_fill_gradientn(colors = c("white", "red"))+ theme(axis.text.y = element_text(size = 4))

dge<-SubsetData(dge_temp,max.cells.per.ident = 100)

DoHeatmap(dge,assay = "RNA",slot = "scale.data", features = c(Response_features),raster = F, size = 3,draw.lines = FALSE)#+ scale_fill_gradientn(colors = c("blue", "red"))
FeaturePlot(dge_temp,features=c("AR_signature1"))
ggsave(file="./heatmap_NE.pdf",width = 10,height = 15,units = "cm")


VlnPlot(dge,group.by = "type",features=c("KRT6A"),assay = "RNA")
ggsave(file="./BE_integrated/KRT6A_vlnplot.eps",width = 10,height = 10,units = "cm")
FeaturePlot(dge,features=c("KRT6A"))
ggsave(file="./KRT6A_Featureplot.eps",width = 10,height = 10,units = "cm")
FeaturePlot(dge,features=c("KRT14"))
ggsave(file="../BE_integrated/KRT14_Featureplot.eps",width = 10,height = 10,units = "cm")
FeaturePlot(dge,features=c("PSCA"))
ggsave(file="./PSCA_Featureplot.eps",width = 10,height = 10,units = "cm")
FeaturePlot(dge,features=c("TFF3"))
ggsave(file="./TFF3_Featureplot.eps",width = 10,height = 10,units = "cm")
FeaturePlot(dge,features=c("KRT13"))
ggsave(file="./KRT13_Featureplot.eps",width = 10,height = 10,units = "cm")
FeaturePlot(dge,features=c("LTF"))
ggsave(file="./LTF_Featureplot.eps",width = 10,height = 10,units = "cm")



av.exp <- AverageExpression(dge)$RNA
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, c("BE","BE_organoid","Hillock_organoid","Club","Club_organoid","hMSC_organoid","Proliferation_organoid","LE","ERGpos_Tumor","ERGneg_Tumor","Tumor_organoid","NE"))
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0.75,limit = c(0.5,1))
ggsave(file="./PCa_organoid_correlation.pdf",width = 20,height = 20,units = "cm")
PCA_organoid<-dge
dge<-readRDS("../organoid/BE_integrated.rds")

V1<-as.vector(FetchData(object = dge,vars = "BE1"))
BOX_df<-NULL
BOX_df$id<-as.vector(dge@active.ident)
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
  geom_signif(comparisons = list(c("Club","Club_organoid")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  geom_signif(comparisons = list(c("BE","BE_organoid")),test = "wilcox.test",map_signif_level=TRUE,y_position = max(V1)+0.25,test.args = "greater") +
  ggtitle("BE1")+ RotatedAxis()
ggsave(file="BE_signature_vlnplot.eps",width = 20,height = 20,units = "cm")


markers<-read.table("~/Desktop/Seqwell_0913/human_output/prostate/geneset/Club.txt",sep="\n")
markers<-read.table("~/Desktop/Seqwell_combined/combined_analysis/PCA_signature/Club_Signature.txt",sep="\n")
markers<-as.vector(markers$V1)

for (i in 2:length(markers)){
  dge_temp<-SubsetData(club_organoid,subset.name = markers[i],low.threshold = 0)
  print(ncol(dge_temp)/ncol(club_organoid))
}





