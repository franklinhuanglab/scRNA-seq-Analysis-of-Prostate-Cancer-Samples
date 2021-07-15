dge_EUR_LE<-SubsetData(dge_LE,cells = colnames(dge_LE)[!(dge_LE$orig.ident %in% c("AUG_PB1","MAY_PB1","MAY_PB2"))])
dge_EUR_LE$orig.ident[dge_EUR_LE$orig.ident %in% c("PR5249_T","PR5249_N")]="PR5249"
dge_EUR_LE$orig.ident[dge_EUR_LE$orig.ident %in% c("PR5251_T","PR5251_N")]="PR5251"
dge_EUR_LE$orig.ident[dge_EUR_LE$orig.ident %in% c("PR5254_T","PR5254_N")]="PR5254"
dge_EUR_LE$orig.ident[dge_EUR_LE$orig.ident %in% c("PR5261_T","PR5261_N")]="PR5261"
table(dge_EUR_LE$orig.ident)
#PR5186 PR5196 PR5199 PR5249 PR5251 PR5254 PR5261 PR5269 
#9    135     19    828   1585    132   1073    329 

#### V1: 5251
dge=list()
dge[[1]]=SubsetData(dge_E,ident.use = "ERGneg_Tumor")
dge[[2]]=SubsetData(dge_E,ident.use = "ERGpos_Tumor")
for (i in 1:length(dge)) {
dge[[i]] <- NormalizeData(dge[[i]], verbose = FALSE)
dge[[i]] <- FindVariableFeatures(dge[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
dge.anchors <- FindIntegrationAnchors(object.list = dge, dims = 1:50)
dge <- IntegrateData(anchorset = dge.anchors, dims = 1:50)

DefaultAssay(dge) <- "integrated"

# Run the standard workflow for visualization and clustering
dge <- ScaleData(dge, verbose = FALSE)
dge <- RunPCA(dge, npcs = 100, verbose = FALSE)
# t-SNE and Clustering
dge <- RunUMAP(dge, reduction = "pca", dims = 1:100)
dge <- FindNeighbors(dge, reduction = "pca", dims = 1:100)
dge <- FindClusters(dge, resolution = 1)
# Visualization

p1 <- DimPlot(dge, reduction = "umap", group.by = "ID")
p2 <- DimPlot(dge, reduction = "umap", label = TRUE)

plot_grid(p1, p2)
ggsave(file="UMAP_combined.eps",width = 40,height = 20,units = "cm")
DimPlot(dge, reduction = "umap", split.by = "ID")
ggsave(file="UMAP_separated.eps",width = 40,height = 20,units = "cm")

dge$ID[dge$ID=="Epithelial"]="LE"

dge$orig.ident[dge$orig.ident %in% c("PR5249_T","PR5249_N")]="PR5249"
dge$orig.ident[dge$orig.ident %in% c("PR5251_T","PR5251_N")]="PR5251"
dge$orig.ident[dge$orig.ident %in% c("PR5254_T","PR5254_N")]="PR5254"
dge$orig.ident[dge$orig.ident %in% c("PR5261_T","PR5261_N")]="PR5261"


dge$type<-as.vector(dge$ID)
dge$type[dge$type %in% c("BE","Club","LE","ERGneg_Tumor","ERGpos_Tumor")]="PCA"
dge$type[dge$type %in% c("BE_normal","Club_Normal","LE_Normal","Hillock_Normal")]="Normal"

dge$ID[dge$ID=="ERGpos_Tumor"]="Tumor"
dge$ID[dge$ID=="ERGneg_Tumor"]="Tumor"

Tab<-table(dge$orig.ident,Idents(object = dge))
write.table(Tab,"table_initialclustering.txt",sep="\t",row.names = T,col.names = T)

Tab<-table(dge$ID,Idents(object = dge))
write.table(Tab,"table_ID.txt",sep="\t",row.names = T,col.names = T)

DefaultAssay(dge) <- "RNA"
dge$ID[dge$ID=="Stroma_Henry"]="Normal"
dge$ID[dge$ID %in% c("Endothelial","Smooth_muscle","Fibroblast")]="PCA"


dge$ID[dge$ID=="BE"]="PCA"
dge$ID[dge$ID=="BE_Henry"]="Normal"

for (i in 1:length(unique(dge@active.ident)))  {
markers_temp <- FindConservedMarkers(dge, ident.1 = unique(dge@active.ident)[i], grouping.var = "ID")
#markers_temp <- FindConservedMarkers(dge, ident.1 = unique(dge@active.ident)[i], grouping.var = "type")
markers_temp$gene<-as.vector(rownames(markers_temp))
top10_1 <- markers_temp %>% top_n(n = 10, wt = Normal_avg_logFC)
top10_2 <- markers_temp %>% top_n(n = 10, wt = PCA_avg_logFC)
#top10_1 <- markers_temp %>% top_n(n = 10, wt = Normal_avg_logFC)
#top10_2 <- markers_temp %>% top_n(n = 10, wt = PCA_avg_logFC)
genelist1<-unique(c(top10_1$gene,top10_2$gene))

cluster_temp <- subset(dge, idents = unique(dge@active.ident)[i])
cluster_temp<-SetIdent(cluster_temp,value = as.vector(cluster_temp$ID))
#cluster_temp<-SetIdent(cluster_temp,value = as.vector(cluster_temp$type))
markers_temp <- FindMarkers(cluster_temp, ident.1 = "Normal", ident.2 = NULL, only.Tumor = TRUE)
markers_temp$gene<-as.vector(rownames(markers_temp))
list_tmp <- markers_temp  %>% top_n(n = 10, wt = abs(avg_logFC))
genelist2<-list_tmp$gene

markers_temp <- FindMarkers(cluster_temp, ident.1 = "PCA", ident.2 = NULL, only.Tumor = TRUE)
markers_temp$gene<-as.vector(rownames(markers_temp))
list_tmp <- markers_temp  %>% top_n(n = 10, wt = abs(avg_logFC))
genelist3<-list_tmp$gene
  
genelist_temp<-unique(c(genelist1,genelist2,genelist3))
avg.cluster_temp <- log1p(AverageExpression(cluster_temp)$RNA)
avg.cluster_temp$gene <- rownames(avg.cluster_temp)
p1 <- ggplot(avg.cluster_temp, aes(Normal,PCA)) + geom_point() + ggtitle(paste0("Cluster_",unique(dge@active.ident)[i]))+stat_cor(label.x = 0, label.y = 6.5) +
  stat_regline_equation(label.x = 1, label.y = 6)
LabelPoints(plot = p1, points = genelist_temp, repel = TRUE)

ggsave(file=paste0("Scatter_expression_cluster_",unique(dge@active.ident)[i],".eps"),width = 20,height = 20,units = "cm")
DimPlot(cluster_temp,group.by = "ID")
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
ggsave(file=paste0("Heatmap_byID_cluster_",unique(dge@active.ident)[i],".eps"),width = 20,height = 20,units = "cm")
}
DefaultAssay(dge) <- "RNA"
pbmc.markers <- FindAllMarkers(dge,only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1,test.use = "wilcox")
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
dge_temp<-SubsetData(dge,max.cells.per.ident = 100)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
DoHeatmap(dge_temp,assay = "RNA", features = top5$gene,raster=F) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_10.pdf",width = 20,height = 20,units = "cm",limitsize = FALSE)
ggsave(file="Heatmap_10.pdf",width = 15,height = 10,units = "cm",limitsize = FALSE)
DotPlot(dge, features = unique(top10$gene),assay = "RNA") + RotatedAxis()
ggsave(file="dot_plot_bypatient.eps",width = 55,height = 15,units = "cm",limitsize = FALSE)
write.table(pbmc.markers,"All_markers.txt",sep="\t",col.names = T,row.names = T)
for (i in 1:length(unique(dge@active.ident))) {
  markers_temp <- FindMarkers(dge, ident.1 = unique(dge@active.ident)[i], ident.2 = NULL, only.pos = TRUE)
  markers_temp$gene<-as.vector(rownames(markers_temp))
  markers_temp<-markers_temp[order(markers_temp$avg_logFC,decreasing = T),]
  sig_cluster<-markers_temp$gene[markers_temp$p_val_adj<0.05]
  sig_cluster<-sig_cluster[1:min(length(sig_cluster),100)]
  write.table(sig_cluster,paste0("",unique(dge@active.ident)[i],"_Signature.txt"),sep = "\t",quote = F,col.names = F,row.names = F)
}

BuildClusterTree(dge)
PlotClusterTree(dge)

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

ggsave(file="Club_stability.eps",width = 20,height = 40,units = "cm",limitsize = FALSE)
clustree(dge, prefix = "integrated_snn_res.",
         node_colour = "LTF", node_colour_aggr = "median")

clustree_overlay(dge, prefix = "integrated_snn_res.", x_value = "umap1", y_value = "umap2")


ERG_pos<-SubsetData(dge,cells = colnames(dge)[dge$ID=="ERGpos_Tumor"])
ERG_neg<-SubsetData(dge,cells=colnames(dge)[dge$ID=="ERGneg_Tumor"])

DefaultAssay(ERG_pos) <- "integrated"
ERG_pos <- ScaleData(ERG_pos, verbose = FALSE)
ERG_pos <- RunPCA(ERG_pos, npcs = 100, verbose = FALSE)
ERG_pos <- RunUMAP(ERG_pos, reduction = "pca", dims = 1:100)

DefaultAssay(ERG_neg) <- "integrated"
ERG_neg <- ScaleData(ERG_neg, verbose = FALSE)
ERG_neg <- RunPCA(ERG_neg, npcs = 100, verbose = FALSE)
ERG_neg <- RunUMAP(ERG_neg, reduction = "pca", dims = 1:100)

DimPlot(ERG_neg)
ggsave(file="UMAP_ERG_neg.eps",width = 20,height = 20,units = "cm")

DimPlot(ERG_pos)
ggsave(file="UMAP_ERG_pos.eps",width = 20,height = 20,units = "cm")

Fisher_table<-table(dge$ID,dge@active.ident)
n_ERG_neg<-sum(Fisher_table[1,])
n_ERG_pos<-sum(Fisher_table[2,])
Fisher_great<-NULL
Fisher_less<-NULL
for (i in 1:ncol(Fisher_table)) {
  dat = rbind(c(Fisher_table[1,i],(n_ERG_neg-Fisher_table[1,i])),
                 + c(Fisher_table[2,i],(n_ERG_pos-Fisher_table[2,i])) )
test <- fisher.test(dat,alternative = "greater")
test2 <- fisher.test(dat,alternative = "less")
Fisher_great<-cbind(Fisher_great,test$p.value)
Fisher_less<-cbind(Fisher_less,test2$p.value)
}
p.adjust(Fisher_great,n = 8,method = "fdr")
p.adjust(Fisher_less,n = 8,method = "fdr")







dge<-SetIdent(dge,value = as.vector(dge$ID_2))
dge_temp<-AverageExpression(dge)
av.exp <- AverageExpression(dge)$integrated
allmarkers<-read.table("../BE_integrated/all_BE_markers.txt",sep="\t",header = T,row.names = 1)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
genelist<-as.vector(top10$gene)
av.exp<-av.exp[rownames(av.exp) %in% genelist,]
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, c('0','1','2','3','4','5','6','7','8'))
cor.df <- tidyr::gather(data = cor.exp, y, correlation, c('0_BE', '0_BE_Henry',
                                                          '1_BE', '1_BE_Henry',
                                                          '2_BE', '2_BE_Henry',
                                                          '3_BE', '3_BE_Henry',
                                                          '4_BE', '4_BE_Henry',
                                                          '5_BE', '5_BE_Henry',
                                                          '6_BE', '6_BE_Henry',
                                                          '7_BE', '7_BE_Henry',
                                                          '8_BE', '8_BE_Henry',))
ggplot(cor.df, aes(x, y, fill = correlation)) +geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")





