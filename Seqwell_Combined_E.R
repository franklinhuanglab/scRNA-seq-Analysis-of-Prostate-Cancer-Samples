
library(DESeq2)
library(edgeR)
library(DEFormats)
library(Seurat)
library(cowplot)
library(dplyr)
library(GEOquery)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggpubr)
library(ggsignif)
theme_set(theme_cowplot())
#samples <- as.vector(read.csv("../Filename.txt",header=FALSE,sep="\n")$V1)
#l=length(samples)

setwd("/Users/hsong/Desktop/Seqwell_combined/")
mkdir("E_analysis")
setwd("./E_analysis/")

mkdir("./E_analysis/BE_analysis")
setwd("./E_analysis/BE_analysis/")

mkdir("./E_analysis/LE_analysis")
setwd("./E_analysis/LE_analysis/")

mkdir("./LE_cluster")
setwd("./LE_cluster")
#Kidney samples

#merge 2 objects
dge_temp$patient=as.vector(dge_temp$orig.ident)
dge_temp$patient[dge_temp$orig.ident %in% c("PR5249_T","PR5249_N")]="PR5249"
dge_temp$patient[dge_temp$orig.ident %in% c("PR5251_T","PR5251_N")]="PR5251"
dge_temp$patient[dge_temp$orig.ident %in% c("PR5254_T","PR5254_N")]="PR5254"
dge_temp$patient[dge_temp$orig.ident %in% c("PR5261_T","PR5261_N")]="PR5261"
VlnPlot(dge_E,features=c("KLK3"),group.by = "patient")

dge<-SubsetData(pca_epithelial,ident.use = c("ERGpos_Tumor","ERGneg_Tumor"))

dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)

dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = dge), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = dge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge = 0,ynudge = 0)
CombinePlots(plots = list(plot1, plot2))
ggsave(file="VariableFeature_2000.pdf",width = 40,height = 20,units = "cm")

all.genes <- rownames(x = dge)
dge <- ScaleData(object=dge,features=VariableFeatures(object = dge))
dge <- RunPCA(object = dge, features = VariableFeatures(object = dge),npcs = 100)

# Create a vector of gene names for downstream analysis
#gene.names.vector <- as.character("gene names")
# Run PCA on selected genes
#dge <- RunPCA(object = dge, features = gene.names.vector, do.print = TRUE, pcs.print = 1:5, 
#                 genes.print = 5)

slot(dge[["pca"]], "misc")
print(x = dge[["pca"]], dims = 1:5, nfeatures = 5)
mat <- Seurat::GetAssayData(dge,assay="RNA",slot="scale.data")
pca <-dge[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)

DimPlot(object = dge, reduction = "pca",group.by = "orig.ident")
ggsave(file="PCA.pdf",width = 30,height = 30,units = "cm")

png("PCA_heatmap.png")
DimHeatmap(object = dge, dims = 1:10, cells = length(dge$nCount_RNA), balanced = TRUE)
dev.off()
png("Elbowplot_dge.png")
ElbowPlot(dge,ndims = 100)
dev.off()

dge <- JackStraw(dge, num.replicate = 100,dims = 75)
dge <- ScoreJackStraw(dge, dims = 1:75) 
JackStrawPlot(dge, dims = 1:75)
ggsave(file="Strawplot.pdf",width = 20,height = 20,units = "cm")
# This is for Epithelial cells
n_pc=100
resolut=2

##test1:
n_pc=100
resolut=1
# This is for BE
n_pc=26
resolut=1
# This is for LE
n_pc=34
resolut=0.8
# This is for Tumor
n_pc=51
resolut=1

dge <- FindNeighbors(dge, dims = 1:n_pc,k.param = 30)
dge <- FindClusters(dge, resolution =c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,2,2.5))



pdf("Stability.pdf",width = 10,height = 15)
clustree(dge, prefix = "RNA_snn_res.")
dev.off()

dge <- FindClusters(dge, resolution = resolut)
#head(Idents(dge), 10)

dge=RunTSNE(dge,dims.use = 1:n_pc,perplexity=40,seed.use = 10)
TSNEPlot(dge,pt.size = 1,label=T)
ggsave(file=paste0("TSNE_raw_",resolut,".pdf"),width = 20,height = 20,units = "cm")


# DimPlot(dge_ERGstatus)
# dge_temp<-dge_ERGstatus
# cell_tumor<-colnames(dge_ERGstatus)[dge_ERGstatus@active.ident %in% c("ERGneg_Tumor","ERGpos_Tumor")]
# dge_tumor<-SubsetData(dge_ERGstatus,cells = cell_tumor)
# DimPlot(dge_tumor)
# dge_tumor_UMAP<-dge_tumor@reductions$umap@cell.embeddings
# dge_tumor_UMAP<-as.data.frame(dge_tumor_UMAP)
# cells_outlier<-NULL
# for (i in 1:nrow(dge_tumor_UMAP)) {
#   if(dge_tumor_UMAP$UMAP_1[i]<(0)) {
#     cells_outlier<-c(cells_outlier,rownames(dge_tumor_UMAP)[i])
#   }
# }
# dge_temp<-SubsetData(dge_ERGstatus,cells = colnames(dge_ERGstatus)[!(colnames(dge_ERGstatus) %in% cells_outlier)])
# dge_ERGstatus <- AddModuleScore(object = dge_ERGstatus, features = list(PD1_pathway), ctrl = 5, name = "PD1_pathway")
# dge_ERGstatus <- AddModuleScore(object = dge_ERGstatus, features = list(Interferon_gamma_pathway), ctrl = 5, name = "Interferon_gamma_pathway")
# 
# cell_tumor<-colnames(dge_ERGstatus)[dge_ERGstatus@active.ident %in% c("ERGneg_Tumor","ERGpos_Tumor")]
# dge_tumor<-SubsetData(dge_ERGstatus,cells = cell_tumor)
# FeaturePlot(dge_ERGstatus,features=c("PD1_pathway1"))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(-0.5, 0.5))
# ggsave(file="PD1_ERGstatus.eps",width = 20,height = 20,units = "cm")
# FeaturePlot(dge_ERGstatus,features=c("Interferon_gamma_pathway1"))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(-0.5, 0.5))
# ggsave(file="Interferon_gamma_group.eps",width = 20,height = 20,units = "cm")

dge <- RunUMAP(dge, dims = 1:n_pc)
DimPlot(dge, reduction = "umap",label=T)
ggsave(file=paste0("UMAP_raw_",resolut,".pdf"),width = 20,height = 20,units = "cm")

DimPlot(dge, reduction = "umap",group.by="orig.ident",label=F)
ggsave(file="Umap_group.pdf",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "tsne",group.by="orig.ident",label=F)
ggsave(file="TSNE_group.pdf",width = 20,height = 20,units = "cm")

# Print number of cells per cluster
Tab<-table(dge$orig.ident,Idents(object = dge))
write.table(tab,"table_initialclustering.txt",sep="\t",row.names = T,col.names = T)
#clustering dendrogram
dge<- BuildClusterTree(dge)
png("Cluster_tree.png")
PlotClusterTree(dge)
dev.off()

saveRDS(dge, file = "LE_aftercluster.rds")

dge<-SetIdent(dge,value = as.vector(dge$orig.ident))
pbmc.markers <- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.001, logfc.threshold = 0.001)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# write.table(pbmc.markers,"all_BEcell_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")
# write.table(pbmc.markers,"all_BEcell_id.txt",col.names = T,row.names = T,quote = F,sep="\t")
# write.table(pbmc.markers,"all_Ecell_markers_bypatient.txt",col.names = T,row.names = T,quote = F,sep="\t")
# 
write.table(pbmc.markers,"all_LEcell_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")
# write.table(pbmc.markers,"all_LEcell_id.txt",col.names = T,row.names = T,quote = F,sep="\t")
# write.table(pbmc.markers,"all_Lcell_markers_bypatient.txt",col.names = T,row.names = T,quote = F,sep="\t")
write.table(pbmc.markers,paste0("all_Ecell_markers_",resolut,".txt"),col.names = T,row.names = T,quote = F,sep="\t")


top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
dge_temp<-SubsetData(dge,max.cells.per.ident = 100)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
DoHeatmap(dge_temp, features = top10$gene,raster = F) + theme(axis.text.y = element_text(size = 5))#+ scale_fill_gradientn(colors = c("white", "red"))
ggsave(file=paste0("Heatmap_10_cluster_",resolut,".pdf"),width = 30,height = 30,units = "cm")
ggsave(file="Heatmap_10_id.pdf",width = 100,height = 100,units = "cm")
ggsave(file="Heatmap_10_patient.pdf",width = 100,height = 100,units = "cm")

top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
features_top10 <- as.vector(unique(top10$gene))

DotPlot(dge, features = c(features_top10),assay = "RNA") + RotatedAxis()

ggsave(file="dot_plot_all_PCa.eps",width = 50,height = 25,units = "cm")

DotPlot(dge, features = c(features_top5,"ERG"),assay = "RNA") + RotatedAxis()

dge<-SetIdent(dge,value = dge$seurat_clusters)
unique(dge$seurat_clusters)
#For dge_E_id
new.cluster.ids <- rep("LE",25)
BE_clusters<-c(0,3,24)
Club_clusters<-c(4)
###15 might be benign? test in infercnv
###10 might be TMPRSS2-ERG+
###9 16 2,7,21,not clear
LE_clusters<-c(1,2,8,11,6,7,10,21)
ERGpos_clusters<-c(5,12,15,17,20)
ERGneg_clusters<-c(14,23,16,9,22,18,13,19)
new.cluster.ids[BE_clusters+1]="BE"
new.cluster.ids[Club_clusters+1]="Club"
new.cluster.ids[ERGpos_clusters+1]="ERGpos_Tumor"
new.cluster.ids[ERGneg_clusters+1]="ERGneg_Tumor"

names(new.cluster.ids) <- levels(dge)
dge <- RenameIdents(dge, new.cluster.ids)
DimPlot(dge,label = T)
ggsave(file="E_marked_umap.eps",width = 20,height = 20,units = "cm")
DimPlot(dge,label = T,reduction = "tsne")
ggsave(file="E_marked_tsne.eps",width = 20,height = 20,units = "cm")
#For BE
new.cluster.ids <- c("BE","BE","KRT13+","BE","KLK3","KRT13+","BE_cluster6","BE_cluster7","BE_cluster8","BE_cluster9")
names(new.cluster.ids) <- levels(dge)
dge <- RenameIdents(dge, new.cluster.ids)
saveRDS(dge, file = "After_BE.rds")
dge_KLK3BE<-SubsetData(dge,ident.use = "KLK3")
dge_KRT13BE<-SubsetData(dge,ident.use = "KRT13+")
VlnPlot(dge,features=c("KRT5","KRT15","KRT17","TP63","KLK3","KRT8","KRT18","KRT13","APOBEC3A","CSTB","LYPD3","SERPINB1"))
ggsave(file="Marker.eps",width = 50,height = 50,units = "cm")


#For LE
new.cluster.ids <- c("LE","LE","LE","LE","LE","LE","LE","LE","LE","LE","KRT5+LE","LE")
names(new.cluster.ids) <- levels(dge)
dge <- RenameIdents(dge, new.cluster.ids)
saveRDS(dge, file = "After_LE.rds")

Tab<-table(dge$orig.ident,Idents(object = dge))
write.table(Tab,"../E_analysis/table_IDclustering.txt",sep="\t",row.names = T,col.names = T)

DimPlot(dge, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(file="TSNE_marked.eps",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "umap",label = TRUE, pt.size = 0.5)
ggsave(file="Umap_marked.eps",width = 20,height = 20,units = "cm")


#######Check some markers among ID clusters
marker_list=c("KRT5","KRT15","KRT17","TP63","KLK3","KRT8","KRT18","KRT13","APOBEC3A","CSTB","LYPD3","SERPINB1")
for (i in 1:length(marker_list)) {
  V1<-as.vector(FetchData(object = dge,vars = marker_list[i],slot = "data"))
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(notch = F) + geom_jitter(shape=16,position = position_jitter(0.1))+
    #geom_signif(comparisons = list(c("Normal","Tumor")),map_signif_level=TRUE,y_position = max(V1)+0.5) +stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    ggtitle(marker_list[i])
  #ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(notch = F) + geom_jitter(shape=16,position = position_jitter(0.1))+geom_signif(comparisons = list(c("Club_PCA","Club_Normal")),map_signif_level=TRUE,y_position = max(V1)+0.5)+
  #  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
  #  ggtitle(marker_list[i])
  ggsave(file=paste0(marker_list[i],"_marker_boxplot.eps"),width = 20,height = 20,units = "cm")
  
}


######### Finished
dge$ID<-as.vector(dge@active.ident)
dge$ID[dge$ID %in% c("LE","BE","Club")]="Epithelial"
dge$refineID<-as.vector(dge@active.ident)
dge<-SetIdent(dge,value = dge$ID)
DimPlot(dge, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(file="TSNE_coarse.eps",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "umap",label = TRUE, pt.size = 0.5)
ggsave(file="Umap_coarse.eps",width = 20,height = 20,units = "cm")

dge_temp<-dge_Club_0
dge_temp<-SubsetData(dge,cells = colnames(dge)[dge$orig.ident %in% c("PR5249_T","PR5249_N","PR5251_T",
                                                                    "PR5251_N","PR5254_T","PR5254_N",
                                                                    "PR5261_T","PR5261_N")])
dge_temp$pairtype<-"PairN"
dge_temp$pairtype[dge_temp$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T")]="PairT"
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$pairtype))

dge_temp$type<-"PCa"
dge_temp$type[dge_temp$orig.ident =="Henry_lung"]="Normal"
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$type))

dge_subset<-SubsetData(dge,max.cells.per.ident = 100)
source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
Data<-SubsetData(dge_E,ident.use = c("KRT13pos_BE","KRT13pos_Club"))
Data<-SubsetData(dge_tumor,cells = colnames(dge_tumor)[dge_tumor$orig.ident %in% c("PR5251_T","PR5251_N")])

Data<-dge_temp
Data<-SetIdent(Data,value = as.vector(Data$ERGtype))
dge<-SubsetData(dge_E,ident.use = c("Club"))
dge$ERGtype<-"ERGneg"
dge$ERGtype[dge$orig.ident %in% c("AUG_PB1","MAY_PB1","PR5249","PR5251","PR5199","PR5269")]="ERGpos"
dge$ID<-as.vector(dge$ERGtype)
dge<-SetIdent(dge,value = as.vector(dge$ID))
Data<-SubsetData(dge_temp,max.cells.per.ident = 2667)
Data<-dge
name1<-"Ident"
name2<-"PR5254_Tcell"
Prepare.for.GSEA(Data,name1,name2)

dge<-SubsetData(dge,ident.use = c("1","8","11","6","10","13","19","18","7","21","22","2","9","16","23","14","5","15","12","17","20"))
dge<-SetIdent(dge,value = as.vector(dge$seurat_clusters))
dge_temp<-SubsetData(dge_E,ident.use = c("LE","ERGpos_Tumor","ERGneg_Tumor"))

dge_temp<-SubsetData(dge_E,ident.use = c("LE","ERGpos_Tumor","ERGneg_Tumor"))
dge_temp$ID=dge_temp$patient
dge_temp$ID[!(dge_temp@active.ident %in% c("ERGpos_Tumor","ERGneg_Tumor"))]="Others"

dge_temp<-SubsetData(dge_E,ident.remove = c("LE"))
dge_infer_mtx <- as.matrix(GetAssayData(dge_temp, slot = "counts"))
write.table(dge_infer_mtx,"Enew_counts_byid.txt",sep="\t",row.names = T,col.names = T,quote=F)
phenotype<-NULL
phenotype$cellname<-colnames(dge_temp)
phenotype$ID<-as.vector(dge_temp$ID)
phenotype$ID<-as.vector(dge_organoid@active.ident)
phenotype<-as.data.frame(phenotype)
#View(phenotype)
write.table(phenotype,"Enew_phenotype_byid.txt",sep="\t",row.names = F,col.names = F,quote=F)

dge_BE<-SubsetData(dge,ident.use = BE_clusters)
dge_Club<-SubsetData(dge,ident.use = Club_clusters)
dge_ERGpos<-SubsetData(dge,ident.use = ERGpos_clusters)
dge_ERGnegandbenignLE<-SubsetData(dge,ident.remove =c("0","3","6","7","11","13") )
dge_LE<-SubsetData(dge,ident.use = c("6","7","11","13","0","3","6","7","11","13"))

saveRDS(dge_Club,"dge_Club.rds")
saveRDS(dge_BE,"dge_BE.rds")
saveRDS(dge_ERGpos,"dge_ERGpos.rds")
saveRDS(dge_LE,"dge_LE.rds")


##### Volcano plot of tumor vs Benign
dge_LE<-SubsetData(dge,cells = cell_LE)
if((length(cell_tumor_ERGp)/length(cell_tumor))>=0.9) {
  print("Majority: ERG+ tumor cells",quote=F)
}
if((length(cell_tumor_ERGn)/length(cell_tumor))>=0.9) {
  print("Majority: ERG- tumor cells",quote=F)
}
normal_LE<-setdiff(colnames(dge_LE),cell_tumor)
dge_LE<-SetIdent(object = dge_LE,cells=cell_tumor,value = "Tumor")
dge_LE<-SetIdent(object = dge_LE,cells=normal_LE,value = "Benign")
dge_LE_backup<-dge_LE
dge_LE <- NormalizeData(object = dge_LE, normalization.method = "LogNormalize", scale.factor = 10000)
dge_LE <- FindVariableFeatures(object = dge_LE, selection.method = "vst", nfeatures = 10000)
dge_LE <- ScaleData(object=dge_LE,features=VariableFeatures(object = dge))


dge_temp<-SubsetData(seurat_all,ident.use = c("LE","ERGneg_Tumor"))
all_markers <-FindAllMarkers(dge_temp,slot="data",min.diff.pct=0.1,min.pct =  0.1,test.use = "wilcox",assay = "RNA")
mito.genes <- grep(pattern = "^MT-", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
all_markers<-all_markers[!(all_markers$gene %in% remove_gene),]
write.table(all_markers,"E_markers.txt",sep="\t",row.names=F, col.names = T)
all_markers<-all_markers[1:(nrow(all_markers)/2),]
DE_select<-all_markers$p_val_adj<0.05
length(which(DE_select))
all_markers$threshold <- DE_select
ggplot(all_markers) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj), colour=threshold)) +
  ggtitle("PCA vs Normal overexpression") +
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
all_markers_reorder$genelabels[1:30] <- rownames(all_markers_reorder)[1:30]
# all_markers_reorder$genelabels[all_markers_reorder$gene =="KLK2"]="KLK2"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="KLK3"]="KLK3"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="PIGR"]="PIGR"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="MMP7"]="MMP7"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="CP"]="CP"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="SCGB1A1"]="SCGB1A1"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="CD74"]="CD74"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="LTF"]="LTF"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="KRT8"]="KRT8"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="KRT18"]="KRT18"
ggplot(all_markers_reorder) +
  geom_point(aes(x = avg_logFC, y = -log10(p_val_adj), colour = threshold)) +
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = ifelse(genelabels !="", genelabels,""))) +
  ggtitle("ERGneg_vs_ERGpos overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave(file="ERGneg_vs_ERGpos_volcano.eps",width = 20,height = 20,units = "cm")

################ Prepare for GSEA analysis##########
source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
dge_E_Henry_subset
Data<-dge_temp
name1<-"Ident"
name2<-"ERG_compare"
Prepare.for.GSEA(Data,name1,name2)

dge_E_Henry_subset<-SubsetData(dge_E_Henry,max.cells.per.ident = 1000)
dge_E_Henry_subset$ID<-as.vector(dge_E_Henry_subset@active.ident)
dge_E_Henry_subset$ID[dge_E_Henry_subset$ID=="LE"]="LE_normal"
dge_E_Henry_subset$ID[dge_E_Henry_subset$ID=="BE"]="BE_normal"
dge_E_Henry_subset$ID[dge_E_Henry_subset$ID=="ClubOE"]="Club_normal"
dge_E_Henry_subset<-SetIdent(dge_E_Henry_subset,value = as.vector(dge_E_Henry_subset$ID))

dge_E_id<-readRDS("Initial_E_aftercluster.rds")
dge_E_id<-SetIdent(dge_E_id,value = as.vector(dge_E_id$refineID))
dge_E_id_subset<-SubsetData(dge_E_id,max.cells.per.ident = 1000)
dge_E_id_subset$ID<-as.vector(dge_E_id_subset@active.ident)
dge_E_id_subset$ID[dge_E_id_subset$ID=="LE"]="LE_PCA"
dge_E_id_subset$ID[dge_E_id_subset$ID=="BE"]="BE_PCA"
dge_E_id_subset$ID[dge_E_id_subset$ID=="Club"]="Club_PCA"
dge_E_id_subset<-SetIdent(dge_E_id_subset,value = as.vector(dge_E_id_subset$ID))

dge=list()
dge[[1]]=dge_E_Henry_subset
dge[[2]]=dge_E_id_subset
for (i in 1:length(dge)) {
  dge[[i]] <- NormalizeData(dge[[i]], verbose = FALSE)
  dge[[i]] <- FindVariableFeatures(dge[[i]], selection.method = "vst", nfeatures = 10000, verbose = FALSE)
}
E.anchors <- FindIntegrationAnchors(object.list = dge, dims = 1:30,anchor.features = 2000)
E_integrated <- IntegrateData(anchorset =E.anchors, dims = 1:30)

DefaultAssay(E_integrated) <- "integrated"
E_integrated <- ScaleData(E_integrated, verbose = FALSE)
E_integrated <- RunPCA(E_integrated, npcs = 50, verbose = FALSE)
E_integrated <- RunUMAP(E_integrated, reduction = "pca", dims = 1:30)
DimPlot(E_integrated, reduction = "umap",label = T)
E_integrated$source<-as.vector(E_integrated@active.ident)
ggsave(file="E_downsampled_Integrated_UMAP.eps",width = 20,height = 20,units = "cm")


BE_integrated<-SubsetData(E_integrated,ident.use = c("BE_normal","BE_PCA"))
Data<-BE_integrated
name1<-"Ident"
name2<-"BE_integrated"
Prepare.for.GSEA(Data,name1,name2)

Club_integrated<-SubsetData(E_integrated,ident.use = c("Club_normal","Club_PCA"))
Data<-Club_integrated
name1<-"Ident"
name2<-"Club_integrated"
Prepare.for.GSEA(Data,name1,name2)

LE_integrated<-SubsetData(E_integrated,ident.use = c("LE_normal","LE_PCA"))
Data<-LE_integrated
name1<-"Ident"
name2<-"LE_integrated"
Prepare.for.GSEA(Data,name1,name2)



library(dplyr)
library(tidyr)
library(broom)
raw_E_score<-read.csv("dge_LE_cluster_ssGSEA.PROJ.gct",sep="\t",header = F)
cellname<-na.omit(raw_E_score[1,3:ncol(raw_E_score)])
raw_E_score<-raw_E_score[-c(1),]
raw_E_score<-raw_E_score[,-c(2)]

raw_E_score<-t(raw_E_score)
raw_E_cluster<-read.csv("dge_LE_cluster_phenotypes.cls",sep=" ",header=F)
raw_E_score_df<-as.data.frame(raw_E_score[2:5648,1:ncol(raw_E_score)])

rownames(raw_E_score_df)<-as.vector(t(cellname))
colnames(raw_E_score_df)<-raw_E_score[1,1:ncol(raw_E_score)]
raw_E_score_df$CLUSTER<-t(raw_E_cluster)
raw_E_score_df <- as.data.frame(sapply(raw_E_score_df, as.numeric))

varNames <- names(raw_E_score_df)[-length(raw_E_score_df)]
mean <- colMeans(raw_E_score_df[,1:length(raw_E_score_df)-1])
master <- list()
i <- 1

## main for loop subsets; lapply calculates t-statistics for all variables in the subset
for (group in unique(raw_E_score_df$CLUSTER)){
  # create a list of t-test result in a given "group" subset
  results <- lapply((1:length(varNames)), FUN = function(x, subset = raw_E_score_df[raw_E_score_df$CLUSTER == group,]) {
    t.test(subset[varNames[x]], mu = mean[x], alternative = "two.sided")
  })
  master[[group+1]] <- results
  i <- i + 1
  print(paste0("Cluster finished for ", group),quote=F)
}



for (i in 1:length(master)) {
  results_tmp=NULL
  for (j in 1:length(varNames)) {
    tmp=NULL
    tmp=cbind(tmp,as.double(master[[i]][[j]]$estimate))
    tmp=cbind(tmp,as.double(master[[i]][[j]]$conf.int[1]))
    tmp=cbind(tmp,as.double(master[[i]][[j]]$conf.int[2]))
    tmp=cbind(tmp,as.double(master[[i]][[j]]$p.value))
    rownames(tmp)=varNames[j]
    results_tmp<-rbind(results_tmp,tmp)
  }
  colnames(results_tmp)<-c("Mean","95%CI_low","95%CI_hi","p_value")
  results_tmp<-as.data.frame(results_tmp)
  results_tmp$FDR<-results_tmp$p_value
  results_tmp$FDR=p.adjust(results_tmp$p_value, method = 'BH', n = length(results_tmp$p_value))
  write.table(results_tmp,paste0(i-1,"_meanscores.txt"),sep="\t",row.names = T,col.names = T,quote=F)
}


for (i in 1:length(master)) {
  tmp<-read.csv(paste0(i-1,"_meanscores.txt"),sep="\t",header = T)
  tmp_sig<-tmp[tmp$FDR<0.01,]
  tmp_sig<-tmp_sig[order(tmp_sig$Mean,decreasing = T),]
  tmp_sig_top30<-as.data.frame(tmp_sig[1:50,])
  tmp_sig_top30$GSEA<-as.vector(rownames(tmp_sig_top30))
  ggplot(tmp_sig_top30) +
    geom_bar( aes(x=reorder(GSEA,Mean), y=Mean), stat="identity", fill="forestgreen", alpha=0.5) +
    geom_errorbar( aes(x=GSEA, ymin=X95.CI_low, ymax=tmp_sig_top30$X95.CI_hi), width=0.4, colour="grey", alpha=0.9, size=1.5) +
    ggtitle(paste0("Top50 ssGSEA ES for Cluster_",i-1))+coord_flip()
  ggsave(file=paste0("Top50_score_Cluster_",i-1,".pdf"),width = 40,height = 30,units = "cm")
  write.table(tmp_sig,paste0("Cluster_",i-1,"_Sigscores.txt"),col.names = T,row.names = T,sep="\t",quote=F)
}

# [1] "PA_PB1"    "PA_PB2"    "PA_PR5186" "PA_PR5196" "PA_PR5199" "AUG_PB1"   "PR5249_N" 
# [8] "PR5249_T"  "PR5251_N"  "PR5251_T"  "PR5254_N"  "PR5254_T"  "PR5261_N"  "PR5261_T" 
# [15] "PR5269"   

dge_E_id$refineID[dge_E_id$refineID=="BE"]="BE_PCA"
dge_E_id$refineID[dge_E_id$refineID=="LE"]="LE_PCA"
dge_E_id$refineID[dge_E_id$refineID=="Club"]="Club_PCA"
dge_E_id<-SetIdent(dge_E_id,value = as.vector(dge_E_id$refineID))
dge_BE_pairN<-SubsetData(dge_E_id,ident.use = "BE",cells = colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")])
dge_LE_pairN<-SubsetData(dge_E_id,ident.use = "LE",cells = colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")])
dge_Club_pairN<-SubsetData(dge_E_id,ident.use = "Club",cells = colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")])

dge_BE_T<-SubsetData(dge_E_id,ident.use = "BE",cells = colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T","PA_PB1","PA_PB2","PA_PR5186","PA_PR5196", "PA_PR5199", "AUG_PB1","PR5269")])
dge_LE_T<-SubsetData(dge_E_id,ident.use = "LE",cells = colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T","PA_PB1","PA_PB2","PA_PR5186","PA_PR5196", "PA_PR5199", "AUG_PB1","PR5269")])
dge_Club_T<-SubsetData(dge_E_id,ident.use = "Club",cells = colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T","PA_PB1","PA_PB2","PA_PR5186","PA_PR5196", "PA_PR5199", "AUG_PB1","PR5269")])

dge_E_pairN<-SubsetData(dge_E_id,ident.use = c("LE","BE","Club"),cells = colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")])
dge_E_T<-SubsetData(dge_E_id,ident.use = c("LE","BE","Club"),cells = colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T","PA_PB1","PA_PB2","PA_PR5186","PA_PR5196", "PA_PR5199", "AUG_PB1","PR5269")])

dge_E_Henry_subset500<-SubsetData(dge_E_Henry_subset,max.cells.per.ident = 500)
dge_E_pairN_subset<-SubsetData(dge_E_pairN,max.cells.per.ident = 500)
dge_E_T_subset<-SubsetData(dge_temp,max.cells.per.ident = 500) #Here dge_temp is the tumor samples in BE, LE and Club

dge_E_pairN_subset$ID<-as.vector(dge_E_pairN_subset@active.ident)
dge_E_pairN_subset$ID[dge_E_pairN_subset$ID=="Club"]="Club_pairN"
dge_E_pairN_subset$ID[dge_E_pairN_subset$ID=="BE"]="BE_pairN"
dge_E_pairN_subset$ID[dge_E_pairN_subset$ID=="LE"]="LE_pairN"
dge_E_pairN_subset<-SetIdent(dge_E_pairN_subset,value = as.vector(dge_E_pairN_subset$ID))
dge_E_T_subset$ID<-as.vector(dge_E_T_subset@active.ident)
dge_E_T_subset$ID[dge_E_T_subset$ID=="Club"]="Club_T"
dge_E_T_subset$ID[dge_E_T_subset$ID=="BE"]="BE_T"
dge_E_T_subset$ID[dge_E_T_subset$ID=="LE"]="LE_T"
dge_E_T_subset<-SetIdent(dge_E_T_subset,value = as.vector(dge_E_T_subset$ID))

dge_E_subset_500<-merge(dge_E_T_subset,y=dge_E_pairN_subset)
dge=list()
dge[[1]]=dge_E_subset_500
dge[[2]]=dge_E_Henry_subset500
for (i in 1:length(dge)) {
  dge[[i]] <- NormalizeData(dge[[i]], verbose = FALSE)
  dge[[i]] <- FindVariableFeatures(dge[[i]], selection.method = "vst", nfeatures = 10000, verbose = FALSE)
}
E.anchors_500 <- FindIntegrationAnchors(object.list = dge, dims = 1:30,anchor.features = 2000)
E_integrated_500_3group <- IntegrateData(anchorset =E.anchors_500, dims = 1:30)

DefaultAssay(E_integrated_500_3group) <- "integrated"
E_integrated_500_3group <- ScaleData(E_integrated_500_3group, verbose = FALSE)
E_integrated_500_3group <- RunPCA(E_integrated_500_3group, npcs = 50, verbose = FALSE)
E_integrated_500_3group <- RunUMAP(E_integrated_500_3group, reduction = "pca", dims = 1:30)
DimPlot(E_integrated_500_3group, reduction = "umap",label = T)
E_integrated_500_3group$source<-as.vector(E_integrated_500_3group@active.ident)
ggsave(file="E_integrated_500_3group_TandHenryUMAP.eps",width = 20,height = 20,units = "cm")

dge<-SubsetData(E_integrated_500_3group,ident.use = c("Club_normal","Club_pairN","Club_T"))
DefaultAssay(dge) <- "RNA"
name_type<-c("Club_normal","Club_pairN","Club_T")
##### This part is the signature score for all the E population
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/sum_subset.gmt.txt",header=FALSE,sep="\t")
for (i in 1:nrow(tmp)) {
  name_feature=as.character(tmp[i,1])
  feature_tmp<-na.omit(tmp[i,])
  Features<-as.vector(feature_tmp[3:ncol(feature_tmp)])
  
  dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature,assay = "RNA")
  #names(x = dge[[]])
  
  V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(notch = F) + #geom_jitter(shape=16,position = position_jitter(0.1))+
    # geom_signif(comparisons = list(c(name_type[2],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.05)+
    # geom_signif(comparisons = list(c(name_type[1],name_type[2])),map_signif_level=TRUE,y_position = max(V1)+0.25) +
    # geom_signif(comparisons = list(c(name_type[1],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.45) +
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    ggtitle(name_feature)+NoLegend()
  ggsave(file=paste0(name_feature,"_PCA_boxplot.eps"),width = 20,height = 20,units = "cm")
}



BE_integrated_500_TandHenry<-SubsetData(E_integrated_500_TandHenry,ident.use = c("BE_normal","BE_PCA"))
Data<-BE_integrated_500_TandHenry
name1<-"Ident"
name2<-"BE_integrated_500_TandHenry"
Prepare.for.GSEA(Data,name1,name2)

Club_integrated_500_TandHenry<-SubsetData(E_integrated_500_TandHenry,ident.use = c("Club_normal","Club_PCA"))
Data<-Club_integrated_500_TandHenry
name1<-"Ident"
name2<-"Club_integrated_500_TandHenry"
Prepare.for.GSEA(Data,name1,name2)

LE_integrated_500_TandHenry<-SubsetData(E_integrated_500_TandHenry,ident.use = c("LE_normal","LE_PCA"))
Data<-LE_integrated_500_TandHenry
name1<-"Ident"
name2<-"LE_integrated_500_TandHenry"
Prepare.for.GSEA(Data,name1,name2)




####### Divide into different ID and annotate by clusters!!!!!
DimPlot(dge_E_id,group.by = "seurat_clusters",label = T)
dge_E_id<-SetIdent(dge_E_id,value = as.vector(dge_E_id$seurat_clusters))
#Pay attention, the order of clusters are like this: 
#[1] "3"  "4"  "11" "0"  "8"  "10" "2"  "1"  "9"  "12" "6"  "16" "7"  "15" "14" "13" "5" 
new.cluster.ids <- c("Club","LE","ERGpos_Tumor","BE","ERGneg_Tumor","ERGneg_Tumor","LE","LE","ERGneg_Tumor","ERGneg_Tumor","ERGpos_Tumor","ERGneg_Tumor","ERGpos_Tumor","LE","LE","ERGpos_Tumor","LE")
names(new.cluster.ids) <- levels(dge_E_id)
dge_E_id <- RenameIdents(dge_E_id, new.cluster.ids)
DimPlot(dge_E_id)
ggsave(file="E_UMAP_IDbyclustersignature.eps",width = 20,height = 20,units = "cm")
DimPlot(dge_E_id,group.by="orig.ident",label = T)
ggsave(file="E_UMAP_IDbysample.eps",width = 20,height = 20,units = "cm")

dge_BE<-SubsetData(dge_E_id,ident.use = "BE_PCA")
dge_BE$ID[dge_BE$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")]="BE_pairedNormal"
dge_BE$ID[dge_BE$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T","PA_PB1","PA_PB2","PA_PR5186","PA_PR5196", "PA_PR5199", "AUG_PB1","PR5269")]="BE_Tumor"
dge_BE<-SetIdent(dge_BE,value = as.vector(dge_BE$ID))

dge_Club<-SubsetData(dge_E_id,ident.use = "Club")
dge_Club$ID[dge_Club$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")]="Club_pairedNormal"
dge_Club$ID[dge_Club$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T","PA_PB1","PA_PB2","PA_PR5186","PA_PR5196", "PA_PR5199", "AUG_PB1","PR5269")]="Club_Tumor"
dge_Club<-SetIdent(dge_Club,value = as.vector(dge_Club$ID))

dge_LE<-SubsetData(dge_E_id,ident.use = c("LE"))
dge_LE$ID[dge_LE$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")]="LE_pairedNormal"
dge_LE$ID[dge_LE$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T","PA_PB1","PA_PB2","PA_PR5186","PA_PR5196", "PA_PR5199", "AUG_PB1","PR5269")]="LE_Tumor"
dge_LE<-SetIdent(dge_LE,value = as.vector(dge_LE$ID))

dge_E_T<-SubsetData(dge_E_id,cells = colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T","PA_PB1","PA_PB2","PA_PR5186","PA_PR5196", "PA_PR5199", "AUG_PB1","PR5269")])

cell_KLK3_BE<-colnames(SubsetData(dge_BE,subset.name = "KLK3",low.threshold = 0))
cell_KRT5_BE<-colnames(SubsetData(dge_BE,subset.name = "KRT5",low.threshold = 0))
cell_KRT15_BE<-colnames(SubsetData(dge_BE,subset.name = "KRT15",low.threshold = 0))
cell_LTF_BE<-colnames(SubsetData(dge_BE,subset.name = "LTF",low.threshold = 0))

dge_KLK3_Club<-SubsetData(dge_Club,subset.name = "KLK3",low.threshold = 0)
dge_KRT5_Club<-SubsetData(dge_Club,subset.name = "KRT5",low.threshold = 0)
dge_KRT15_Club<-SubsetData(dge_Club,subset.name = "KRT15",low.threshold = 0)
dge_LTF_Club<-SubsetData(dge_Club,subset.name = "LTF",low.threshold = 0)



dge_infer_mtx <- as.matrix(GetAssayData(dge, slot = "counts"))
write.table(dge_infer_mtx,"LE_counts.txt",sep="\t",row.names = T,col.names = T,quote=F)
phenotype<-NULL
phenotype$cellname<-colnames(dge_LE)
phenotype$ID<-as.vector(dge_LE$orig.ident)
phenotype$ID<-as.vector(dge@active.ident)
phenotype<-as.data.frame(phenotype)
write.table(phenotype,"LE_phenotype_bycluster.txt",sep="\t",row.names = F,col.names = T,quote=F)

phenotype<-NULL
phenotype$cellname<-colnames(dge_LE)
phenotype$ID<-as.vector(dge_LE$seurat_clusters)
phenotype<-as.data.frame(phenotype)
#View(phenotype)
write.table(phenotype,"LE_phenotype_byrawcluster.txt",sep="\t",row.names = F,col.names = T,quote=F)

dge_LE$annotation<-as.vector(dge_LE$seurat_clusters)
dge_LE$annotation[dge_LE$seurat_clusters %in% c("1","2","4","5","14","15")]="Non-Malignant_LE"
dge_LE$annotation[dge_LE$seurat_clusters %in% c("6","7","11","13")]="ERGpos_Tumor"
dge_LE$annotation[dge_LE$seurat_clusters %in% c("8","9","10","12","16")]="ERGneg_Tumor"

phenotype<-NULL
phenotype$cellname<-colnames(dge_LE)
phenotype$ID<-as.vector(dge_LE$annotation)
phenotype<-as.data.frame(phenotype)
#View(phenotype)
write.table(phenotype,"LE_phenotype_byannotation.txt",sep="\t",row.names = F,col.names = T,quote=F)

AUG_PB1<-SubsetData(dge_E_id,cells = colnames(dge_E_id)[dge_E_id$orig.ident=="AUG_PB1"])
May_PB1<-SubsetData(dge_E_id,cells = colnames(dge_E_id)[dge_E_id$orig.ident=="PA_PB1"])
May_PB2<-SubsetData(dge_E_id,cells = colnames(dge_E_id)[dge_E_id$orig.ident=="PA_PB2"])
PR5186<-SubsetData(dge_E_id,cells=colnames(dge_E_id)[dge_E_id$orig.ident=="PA_PR5186"])
PR5196<-SubsetData(dge_E_id,cells=colnames(dge_E_id)[dge_E_id$orig.ident=="PA_PR5196"])
PR5199<-SubsetData(dge_E_id,cells=colnames(dge_E_id)[dge_E_id$orig.ident=="PA_PR5199"])
PR5269<-SubsetData(dge_E_id,cells=colnames(dge_E_id)[dge_E_id$orig.ident=="PR5269"])
PR5249<-SubsetData(dge_E_id,cells=colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5249_T","PR5249_N")])
PR5251<-SubsetData(dge_E_id,cells=colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5251_T","PR5251_N")])
PR5254<-SubsetData(dge_E_id,cells=colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5254_T","PR5254_N")])
PR5261<-SubsetData(dge_E_id,cells=colnames(dge_E_id)[dge_E_id$orig.ident %in% c("PR5261_T","PR5261_N")])


DimPlot(AUG_PB1)
ggsave(file="AUG_PB1.eps",width = 20,height = 20,units = "cm")

DimPlot(May_PB1)
ggsave(file="May_PB1.eps",width = 20,height = 20,units = "cm")

DimPlot(May_PB2)
ggsave(file="May_PB2.eps",width = 20,height = 20,units = "cm")

DimPlot(May_PB2)
ggsave(file="May_PB2.eps",width = 20,height = 20,units = "cm")

DimPlot(PR5186)
ggsave(file="PR5186.eps",width = 20,height = 20,units = "cm")

DimPlot(PR5196)
ggsave(file="PR5196.eps",width = 20,height = 20,units = "cm")

DimPlot(PR5199)
ggsave(file="PR5199.eps",width = 20,height = 20,units = "cm")

DimPlot(PR5269)
ggsave(file="PR5269.eps",width = 20,height = 20,units = "cm")

Tab<-table(PR5249$orig.ident,Idents(object = PR5249))
write.table(Tab,"PR5249_clustering.txt",sep="\t",row.names = T,col.names = T)
PR5249$ID<-paste0(as.vector(PR5249$orig.ident),"_",as.vector(PR5249@active.ident))
PR5249<-SetIdent(PR5249,value = as.vector(PR5249$ID))
DimPlot(PR5249,cols = c("gold1","blue1","green3","gray1","coral","darkcyan","red2","pink4","aquamarine2","blueviolet","gray1","gold4","darkolivegreen4"))
ggsave(file="PR5249.eps",width = 20,height = 20,units = "cm")

Tab<-table(PR5251$orig.ident,Idents(object = PR5251))
write.table(Tab,"PR5251_clustering.txt",sep="\t",row.names = T,col.names = T)
PR5251$ID<-paste0(as.vector(PR5251$orig.ident),"_",as.vector(PR5251@active.ident))
PR5251<-SetIdent(PR5251,value = as.vector(PR5251$ID))
DimPlot(PR5251,cols = c("gold1","blue1","green3","gray1","coral","darkcyan","red2","pink4","aquamarine2","blueviolet","gray1","gold4","darkolivegreen4"))
ggsave(file="PR5251.eps",width = 20,height = 20,units = "cm")

Tab<-table(PR5254$orig.ident,Idents(object = PR5254))
write.table(Tab,"PR5254_clustering.txt",sep="\t",row.names = T,col.names = T)
PR5254$ID<-paste0(as.vector(PR5254$orig.ident),"_",as.vector(PR5254@active.ident))
PR5254<-SetIdent(PR5254,value = as.vector(PR5254$ID))
DimPlot(PR5254,cols = c("gold1","blue1","green3","gray1","coral","darkcyan","red2","pink4","aquamarine2","blueviolet","gray1","gold4","darkolivegreen4"))
ggsave(file="PR5254.eps",width = 20,height = 20,units = "cm")

Tab<-table(PR5261$orig.ident,Idents(object = PR5261))
write.table(Tab,"PR5261_clustering.txt",sep="\t",row.names = T,col.names = T)
PR5261$ID<-paste0(as.vector(PR5261$orig.ident),"_",as.vector(PR5261@active.ident))
PR5261<-SetIdent(PR5261,value = as.vector(PR5261$ID))
DimPlot(PR5261,cols = c("gold1","blue1","green3","gray1","coral","darkcyan","red2","pink4","aquamarine2","blueviolet","gray1","gold4","darkolivegreen4"))
ggsave(file="PR5261.eps",width = 20,height = 20,units = "cm")

Tab_LE<-table(dge_E_id$orig.ident,Idents(dge_E_id))
write.table(Tab_LE,"All_E_clustering.txt",sep="\t",row.names = T,col.names = T)

dge_all<-readRDS("../dge_combined_annotated.rds")
DimPlot(dge_all)
##### This part is the signature score for all the E population
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/prostate.gmt",header=FALSE,sep="\t")
name_type<-c("Epithelial","Myeloid","Endothelial","Smooth_muscle","Fibroblast","T-cell")
for (i in 1:nrow(tmp)) {
  name_feature=as.character(tmp[i,1])
  feature_tmp<-na.omit(tmp[i,])
  Features<-as.vector(feature_tmp[3:ncol(feature_tmp)])
  
  dge_all <- AddModuleScore(object = dge_all, features = list(Features), name = name_feature)
  #names(x = dge[[]])
  
  V1<-as.vector(FetchData(object = dge_all,vars = paste0(name_feature,"1")))
  BOX_df<-NULL
  BOX_df$id<-dge_all@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(notch = F) + geom_jitter(shape=16,position = position_jitter(0.1))+
    geom_signif(comparisons = list(c(name_type[1],name_type[2])),map_signif_level=TRUE,y_position = max(V1)+0.05)+
    geom_signif(comparisons = list(c(name_type[1],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.15) +
    geom_signif(comparisons = list(c(name_type[1],name_type[4])),map_signif_level=TRUE,y_position = max(V1)+0.25) +
    geom_signif(comparisons = list(c(name_type[1],name_type[5])),map_signif_level=TRUE,y_position = max(V1)+0.35) +
    geom_signif(comparisons = list(c(name_type[1],name_type[6])),map_signif_level=TRUE,y_position = max(V1)+0.45) +
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    ggtitle(name_feature)
  ggsave(file=paste0(name_feature,"_PCA_boxplot.eps"),width = 30,height = 20,units = "cm")
}
DotPlot(dge_all,features=c("KRT5","KLK3","EPCAM","SELE","VWF","IFI27","RGS5","ACTA2","MYH11","DCN","C7","C1S","IGJ","MZB1","IGKC","TPSAB1","KIT","CPA3","IL7R","CD3D","CD8A","IL1B","LYZ","APOE","CD79A","CD22","MS4A1"),)
ggsave(file="../Combined_figure/Combined_dotplot.eps",width = 50,height = 20,units = "cm")

dge_temp<-SubsetData(KLK3_main,ident.use = c("KLK3","BE_Normal"))
all_markers <-FindAllMarkers(dge_temp,min.diff.pct=0.00001,min.pct =  0.00001,test.use = "wilcox",assay = "RNA")
mito.genes <- grep(pattern = "^MT-", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
RPL.genes <- grep(pattern = "^RPL", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)
all_markers<-all_markers[!(all_markers$gene %in% remove_gene),]
write.table(all_markers,"Club_vs_otherE.txt",sep="\t",row.names=F, col.names = T)
all_markers<-all_markers[1:(nrow(all_markers)/2),]
DE_select<-all_markers$p_val_adj<1
length(which(DE_select))
all_markers$threshold <- DE_select
ggplot(all_markers) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj), colour=threshold)) +
  ggtitle("Tumor vs Normal overexpression") +
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
all_markers_reorder$genelabels[1:30] <- rownames(all_markers_reorder)[1:30]
all_markers_reorder$genelabels[all_markers_reorder$gene =="KLK2"]="KLK2"
all_markers_reorder$genelabels[all_markers_reorder$gene =="KLK3"]="KLK3"
all_markers_reorder$genelabels[all_markers_reorder$gene =="PIGR"]="PIGR"
all_markers_reorder$genelabels[all_markers_reorder$gene =="MMP7"]="MMP7"
all_markers_reorder$genelabels[all_markers_reorder$gene =="CP"]="CP"
all_markers_reorder$genelabels[all_markers_reorder$gene =="SCGB1A1"]="SCGB1A1"
all_markers_reorder$genelabels[all_markers_reorder$gene =="CD74"]="CD74"
all_markers_reorder$genelabels[all_markers_reorder$gene =="LTF"]="LTF"
all_markers_reorder$genelabels[all_markers_reorder$gene =="KRT8"]="KRT8"
all_markers_reorder$genelabels[all_markers_reorder$gene =="KRT18"]="KRT18"
all_markers_reorder$genelabels[all_markers_reorder$gene =="KRT13"]="KRT13"
all_markers_reorder$genelabels[all_markers_reorder$gene =="IL8"]="IL8"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="CLDN4"]="CLDN4"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="CDKN2B"]="CDKN2B"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="SERPINB1"]="SERPINB1"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="IL1RN"]="IL1RN"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="MUC3A"]="MUC3A"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="APOBEC2"]="APOBEC2"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="CLCA4"]="CLCA4"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="PPARG"]="PPARG"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="DAPP1"]="DAPP1"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="LYPD3"]="LYPD3"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="TMPRSS4"]="TMPRSS4"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="C9orf16"]="C9orf16"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="CLDN7"]="CLDN7"
# all_markers_reorder$genelabels[all_markers_reorder$gene =="AQP3"]="AQP3"
ggplot(all_markers_reorder) +
  geom_point(aes(x = avg_logFC, y = -log10(p_val_adj), colour = threshold)) +
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = ifelse(genelabels !="", genelabels,""))) +
  ggtitle("Club vs Other Epithelial overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave(file="Club_vs_Other_Epithelial.eps",width = 20,height = 20,units = "cm")



dge=list()
dge[[1]]=dge_subset
dge[[2]]=dge_Henry_Club
for (i in 1:length(dge)) {
  dge[[i]] <- NormalizeData(dge[[i]], verbose = FALSE)
  dge[[i]] <- FindVariableFeatures(dge[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
club.anchors <- FindIntegrationAnchors(object.list = dge, dims = 1:30,anchor.features = 2000)
club <- IntegrateData(anchorset = club.anchors, dims = 1:30)

DefaultAssay(club) <- "integrated"
club <- ScaleData(club, verbose = FALSE)
club <- RunPCA(club, npcs = 50, verbose = FALSE)
club <- RunUMAP(club, reduction = "pca", dims = 1:30)
DimPlot(club)

dge<-club
DefaultAssay(dge) <- "RNA"
dge$ID<-as.vector(dge@active.ident)
dge$ID[dge$ID %in% c("KRT13pos_Club","Club")]="Club"
dge<-SetIdent(dge,value = as.vector(dge$ID))
Featurename<-c("BE","Club","Luminal","TMPRSS2ERG","LIU_PROSTATE","WALLACE_PROSTATE_RACE_UP","TOMLIN_PROSTATE","NE","Hillock","lungclub")
Features<-list()
Features[[1]] <- as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/BE.txt",header=FALSE,sep="\n")$V1)
Features[[2]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Club.txt",header=FALSE,sep="\n")$V1)
Features[[3]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/LE.txt",header=FALSE,sep="\n")$V1)
Features[[4]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/TMPRSS2ERG.txt",header=FALSE,sep="\n")$V1)
Features[[5]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/LIU_PROSTATE.txt",header=FALSE,sep="\n")$V1)
Features[[6]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/WALLACE_PROSTATE_RACE.txt",header=FALSE,sep="\n")$V1)
Features[[7]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/TOMLIN_PROSTATE.txt",header=FALSE,sep="\n")$V1)
Features[[8]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/NE.txt",header=FALSE,sep="\n")$V1)
Features[[9]]<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Hillock.txt",header=FALSE,sep="\n")$V1)
Features[[10]]<-c("lungclub","SCGB1A1","BPIFA1","SCGB3A2","SFTPD","ALOX15","PON1","CYP4A12B","WFDC1","MGST1","AKR1C18","SFTPA1","LYPD2","CES1F","GSTO1","TST","FAFR4","HP","CHAD","SELENBP1","MGAT3","GABRP")

dge<-dge_E
Featurename<-c("adult_SC","NhESC","PhESC")
Features<-list()

Features[[1]]<- Feature_ASC
Features[[2]]<- Feature_NhESC
Features[[3]]<- Feature_PhESC

for (i in 1:length(Features)) {
  gene.set <- Features[[i]][1:length(Features[[i]])]
  dge <- AddModuleScore(object = dge, features = list(gene.set), ctrl = 5, name = Featurename[i])
  V1<-as.vector(FetchData(object = dge,vars = paste0(Featurename[i],"1")))
  FeaturePlot(dge,features=paste0(c(Featurename[i]),"1"))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))
  ggsave(file=paste0("Modulescore_cluster_",Featurename[i],"_feature.eps"),width = 20,height = 20,units = "cm")
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() +#geom_boxplot(outlier.size = 0.1)+ #geom_jitter(shape=16,position = position_jitter(0.1))+
    #stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
    stat_compare_means(method = "anova",label.x = 3,label.y = max(V1)+0.05)+
    # geom_signif(comparisons = list(c(name_type[1],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.05) +
    # geom_signif(comparisons = list(c(name_type[2],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.25) +
    # geom_signif(comparisons = list(c(name_type[1],name_type[2])),map_signif_level=TRUE,y_position = max(V1)+0.45) +
    ggtitle(Featurename[i])
  ggsave(file=paste0("Modulescore_cluster_",Featurename[i],"_boxplot.eps"),width = 20,height = 15,units = "cm")
  
}
# quantile(dge$adult_SC1, c(0.90))
cell_ASC_organoid<-colnames(dge)[dge$adult_SC1>(mean(dge$adult_SC1)+sd(dge$adult_SC1))]
cell_NhESC_organoid<-colnames(dge)[dge$NhESC1>(mean(dge$NhESC1)+sd(dge$NhESC1)) ]
cell_PhESC_organoid<-colnames(dge)[dge$PhESC1>max((mean(dge$PhESC1)+sd(dge$PhESC1)),0)  ]

cell_ASC_PCA<-colnames(dge)[dge$adult_SC1>(mean(dge$adult_SC1)+sd(dge$adult_SC1))]
cell_NhESC_PCA<-colnames(dge)[dge$NhESC1>(mean(dge$NhESC1)+sd(dge$NhESC1)) ]
cell_PhESC_PCA<-colnames(dge)[dge$PhESC1>max((mean(dge$PhESC1)+sd(dge$PhESC1)),0)  ]

name="PhESC_organoid"
dge<-SubsetData(dge_organoid,ident.use = c("hMSC","BE"))
dge$SC<-"Others"
dge$SC[colnames(dge) %in% cell_PhESC_organoid]="PhESC_organoid"
dge<-SetIdent(dge,value = as.vector(dge$SC))
pbmc.markers <- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
top50 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
write.table(top50,paste0("E_markers.txt"),col.names = T,row.names = T,quote = F,sep="\t")



my_levels <- c("0", "1", "2", "3", "4", "5", "6","7","8","9","10","11","12","13","14","15","16")

# Relevel object@timepoint
dge@active.ident<-factor(dge@active.ident,levels=my_levels)

##### Compare between Gleason and LN involvement
dge_backup<-SubsetData(dge_E_id,ident.use = c("ERGpos_Tumor","ERGneg_Tumor"))
dge_backup<-dge_Club
dge_backup<-SetIdent(dge_backup,value = as.vector(dge_backup$orig.ident))
dge<-SubsetData(dge_backup,ident.use = c("PA_PR5186","PA_PR5196","PA_PR5199","PR5249_N","PR5249_T","PR5251_N","PR5251_T","PR5254_T","PR5254_N","PR5261_T","PR5261_N","PR5269"))
dge<-SubsetData(dge_backup,ident.use = c("PA_PR5186","PA_PR5196","PA_PR5199","PR5249_T","PR5251_T","PR5254_T","PR5261_T","PR5269"))

dge$gleason<-"Gleason7_T"
dge$gleason[dge$orig.ident %in% c("PR5249_N","PR5251_N","PR5261_N")]="Gleason7_N"
dge$gleason[dge$orig.ident %in% c("PR5254_N")]="Gleason9_N"
dge$gleason[dge$orig.ident %in% c("PA_PR5186","PR5254_T")]="Gleason9_T"
level_order <- c("Gleason7_N","Gleason9_N","Gleason7_T","Gleason9_T")
dge$LN="Neg_T"
dge$LN[dge$orig.ident %in% c("PR5249_N","PR5251_N","PR5261_N")]="Neg_N"
dge$LN[dge$orig.ident %in% c("PR5249_N")]="Pos_N"
dge$LN[dge$orig.ident %in% c("PA_PR5186","PR5249_T")]="Pos_T"

dge<-SetIdent(dge,value = as.vector(dge$gleason))
dge<-SetIdent(dge,value = as.vector(dge$LN))
mkdir("/Users/hsong/Desktop/Seqwell_combined/E_analysis/Modulescore_comparison/Tumor_LN")
setwd("/Users/hsong/Desktop/Seqwell_combined/E_analysis/Modulescore_comparison/Tumor_LN")
name_type=unique(dge$gleason)
name_type=unique(dge$LN)
dge<-SubsetData(dge,max.cells.per.ident = min(table(Idents(dge))))
tmp<-read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/sum_subset.gmt.txt",header=FALSE,sep="\t")
for (i in 1:nrow(tmp)) {
  name_feature=as.character(tmp[i,1])
  feature_tmp<-na.omit(tmp[i,])
  Features<-as.vector(feature_tmp[3:ncol(feature_tmp)])
  
  dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature,assay = "RNA")
  #names(x = dge[[]])
  
  V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  BOX_df<-BOX_df[order(match(BOX_df$id, level_order)),]
  ggplot(BOX_df, aes(reorder(id),value,fill=id)) + geom_boxplot(notch = F) + geom_jitter(shape=16,position = position_jitter(0.1))+
    # geom_signif(comparisons = list(c("Gleason7_T","Gleason9_T")),test = "wilcox.test",map_signif_level=TRUE,y_position = max(V1)+0.25) +
    # geom_signif(comparisons = list(c("Gleason7_N","Gleason9_N")),test = "wilcox.test",map_signif_level=TRUE,y_position = max(V1)+0.15) +
    # geom_signif(comparisons = list(c("Gleason9_T","Gleason9_N")),test = "wilcox.test",map_signif_level=TRUE,y_position = max(V1)+0.05)+
    # geom_signif(comparisons = list(c("Gleason7_T","Gleason7_N")),test = "wilcox.test",map_signif_level=TRUE,y_position = max(V1)+0.05) +
    # stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    geom_signif(comparisons = list(c("Pos_T","Neg_T")),test = "wilcox.test",map_signif_level=TRUE,y_position = max(V1)+0.35) +
    geom_signif(comparisons = list(c("Pos_N","Neg_N")),test = "wilcox.test",map_signif_level=TRUE,y_position = max(V1)+0.2) +
    geom_signif(comparisons = list(c("Pos_T","Pos_N")),test = "wilcox.test",map_signif_level=TRUE,y_position = max(V1)+0.05)+
    geom_signif(comparisons = list(c("Neg_T","Neg_N")),test = "wilcox.test",map_signif_level=TRUE,y_position = max(V1)+0.05) +
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    ggtitle(name_feature)
  
  ##### edit this for what you want to test and save as sig or insig. 
  #### sig means a significant difference between combined (matched normal)
  # value_first<-BOX_df$value[BOX_df$id %in% c("Gleason7_T","Gleason7_N")]
  # value_second<-BOX_df$value[BOX_df$id %in% c("Gleason9_T","Gleason9_N")]
  value_first<-BOX_df$value[BOX_df$id %in% c("Pos_N","Pos_T")]
  value_second<-BOX_df$value[BOX_df$id %in% c("Neg_N","Neg_T")]
  stat<-wilcox.test(value_first,value_second)
  if (stat$p.value>0.05) {
    ggsave(file=paste0("insig_",name_feature,"_boxplot.eps"),width = 20,height = 20,units = "cm")
  } else {
      ggsave(file=paste0("sig_",name_feature,"_boxplot.eps"),width = 20,height = 20,units = "cm")
  }
  
  
}
cell_tumor_ERGp<-colnames(dge_ERGpos)
dge_all<-RenameCells(dge_all,new.names = as.vector(dge_all$sample))
DimPlot(dge_all, label=T, cells.highlight= list(colnames(dge_ERGpos),intersect(colnames(SubsetData(dge_all,ident.use = "Epithelial")),colnames(dge_ERGneg))), cols.highlight = c("red","purple"), cols= "grey")

gene_NE<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/NE.txt",header=FALSE,sep="\n")$V1)
gene_NE<-gene_NE[2:length(gene_NE)]
gene_Hillock<-as.vector(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Hillock.txt",header=FALSE,sep="\n")$V1)
gene_Hillock<-gene_Hillock[2:length(gene_Hillock)]
dge<-dge_organoid
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x = dge)
dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = dge), 10)
top10
all.genes <- rownames(x = dge)
dge <- ScaleData(object=dge,features=gene_Hillock)
dge <- ScaleData(object=dge,features=gene_NE)
dge <- RunPCA(object = dge, features = gene_Hillock,npcs = 100)
dge <- RunPCA(object = dge, features = gene_NE, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
slot(dge[["pca"]], "misc")
print(x = dge[["pca"]], dims = 1:5, nfeatures = 5)
mat <- Seurat::GetAssayData(dge,assay="RNA",slot="scale.data")
pca <-dge[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)
#varExplained would be the % of variance explained by the PC
#If we choose PC1 as the signature score of a specific cell type, use dge@reductions$pca@cell.embeddings to select cells
df<-as.data.frame(dge@reductions$pca@cell.embeddings[,1])
View(df)
colnames(df)<-"Hillock_scores"
colnames(df)<-"NE_scores"

png("PCA_heatmap_NE.png")
hist(df, labels = T,ylim=c(0,2000),ylab = "Freq of cells",xlab="NE_scores")
dev.off()

DimPlot(object = dge, reduction = "pca")
ggsave(file="../PCA_NE.eps",width = 30,height = 30,units = "cm")

png("PCA_heatmap.png")
DimHeatmap(object = dge, dims = 1:10, cells = length(dge$nCount_RNA), balanced = TRUE)
dev.off()
df_NEscores<-df
NE_cutoff<-quantile(df$NE_scores, c(.99, .995, .999)) 
cell_NE<-rownames(df)[df$NE_scores>NE_cutoff[2]]
dge_NE_organoid<-SubsetData(dge,cells=cell_NE)



########### BE and LE integration
sample_tumor=c("AUG_PB1","PR5249_T","PR5251_T","PR5254_T","PR5261_T","PR5269","MAY_PB1","MAY_PB2","PR5186","PR5196","PR5199")
sample_pairedN<-c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")
dge_BEmain$ID[dge_BEmain$orig.ident %in% sample_pairedN]="dge_BE_Pairednormal"
dge_BEmain$ID[dge_BEmain$orig.ident %in% sample_tumor]="dge_BE_Tumor"
dge_BEmain<-SetIdent(dge_BEmain,value = as.vector(dge_BEmain$ID))
dge_BE_subset<-SubsetData(dge_BEmain,max.cells.per.ident = 500)

dge_Henry_BE_subset
dge_Henry_Hillock_subset

dge_KRT13pos_BE
dge_BE<-readRDS("../E_analysis/BE_analysis/After_BE.rds")


dge_Henry<-readRDS("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/Henry/Henry_clusteredandid.rds")
DimPlot(dge_Henry)
dge_Henry$ID<-paste0(dge_Henry$ID,"_Normal")
dge_Henry<-SetIdent(dge_Henry,value = as.vector(dge_Henry$ID))
dge_Henry<-AddModuleScore(object = dge_Henry,features=list(Features),ctrl = 2,name=Featurename)
FeaturePlot(dge_Henry,features=c("TMPRSS2_ACE21"))+scale_color_gradientn(colours = c("blue","green","yellow","red"),  limits = c(0, max(dge_Henry$TMPRSS2_ACE21)))



dge=list()
dge[[1]]=dge_E
dge[[2]]=dge_Henry
for (i in 1:length(dge)) {
  dge[[i]] <- NormalizeData(dge[[i]], verbose = FALSE)
  dge[[i]] <- FindVariableFeatures(dge[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
E_main.anchors <- FindIntegrationAnchors(object.list = dge, dims = 1:30,anchor.features = 2000)
E_main <- IntegrateData(anchorset = E_main.anchors, dims = 1:30)


DefaultAssay(E_main) <- "integrated"
E_main <- ScaleData(E_main, verbose = FALSE)
E_main <- RunPCA(E_main, npcs = 50, verbose = FALSE)
E_main <- RunUMAP(E_main, reduction = "pca", dims = 1:50)
E_main$ID<-as.vector(E_main@active.ident)
E_main$ID[E_main$ID %in% c("ERGneg_Tumor_PairN","ERGneg_Tumor_T")]="ERGneg_Tumor"
E_main$ID[E_main$ID %in% c("ERGpos_Tumor_PairN","ERGpos_Tumor_T")]="ERGpos_Tumor"
E_main<-SetIdent(E_main,value = as.vector(E_main$ID))
DimPlot(E_main)
my_levels <- c("BE_Normal","BE_PairN","BE_T","LE_Normal","LE_PairN","LE_T","ClubOE_Normal","Club_PairN","Club_T","ERGpos_Tumor","ERGneg_Tumor","HillockOE_Normal")
# Relevel object@ident
E_main@active.ident <- factor(x = E_main@active.ident, levels = my_levels)

VlnPlot(E_main,features=c("TMPRSS2"),assay = "RNA")
ggsave(file="TMPRSS2_vlnplot.eps",width = 60,height = 20,units = "cm")
VlnPlot(E_main,features=c("ACE2"),assay = "RNA")
ggsave(file="ACE2_vlnplot.eps",width = 60,height = 20,units = "cm")



LE_main$ID<-as.vector(LE_main$orig.ident)
LE_main$ID[LE_main$ID %in% samplelist_t]="LE_Tumor"
LE_main$ID[LE_main$ID %in% samplelist_n]="LE_Pairednormal"
LE_main$ID[LE_main$ID=="Henry"]="LE_Normal"
LE_main<-SetIdent(LE_main,value = as.vector(LE_main$ID))
DimPlot(LE_main, reduction = "umap",laLEl = T)
ggsave(file="LE_main_integrated_subset.eps",width = 20,height = 20,units = "cm")


dge<-SM
V1<-t(as.vector(FetchData(object = dge,vars = "CXCL8")))
V2<-t(as.vector(FetchData(object = dge,vars = "IL8")))
V=rep(0,length(V1))
for (i in 1:length(V1)) {
  V[i]=max(V1[i],V2[i])
}

dge$CXCL8_IL8<-as.vector(V)
name_type<-unique(dge@active.ident)
V1<-as.vector(FetchData(object = dge,vars = "CXCL8_IL8"))
BOX_df<-NULL
BOX_df$id<-dge@active.ident
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(notch = F) + geom_jitter(shape=16,position = position_jitter(0.1))+
  geom_signif(comparisons = list(c(name_type[2],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.05)+
  geom_signif(comparisons = list(c(name_type[1],name_type[2])),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  geom_signif(comparisons = list(c(name_type[1],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+1) +
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
  ggtitle("CXCL8_IL8")
ggsave(file="CXCL8_IL8_SM.eps",width = 20,height = 20,units = "cm")

dge$sample=as.vector(dge$orig.ident)
dge$sample[dge$orig.ident %in% c("PR5251_N","PR5251_T")]="PR5251"
dge$sample[dge$orig.ident %in% c("PR5254_N","PR5254_T")]="PR5254"
dge$sample[dge$orig.ident %in% c("PR5261_N","PR5261_T")]="PR5254"
dge<-SetIdent(dge,value = as.vector(dge$sample))
dge<-SubsetData(dge,ident.remove = c("PR"))


pbmc.markers <- FindAllMarkers(dge_pca, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.01)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge_pca@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge_pca@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge_pca@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
pbmc.markers<-pbmc.markers[pbmc.markers$avg_logFC>=0.5,]
pbmc.markers<-pbmc.markers[pbmc.markers$p_val_adj<0.05,]
write.table(pbmc.markers,"id_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")
DF=NULL
i=1
for (gene in as.vector(pbmc.markers$gene)) {
  DF$gene[i]=gene
  DF$cluster[i]=as.character(pbmc.markers$cluster[i])
  dge_temp=SubsetData(dge_pca,subset.name = gene,low.threshold = 0,ident.use = as.character(DF$cluster[i]))
  DF$percentage[i]=ncol(dge_temp)/length(dge_pca@active.ident[dge_pca@active.ident==DF$cluster[i]])
  i=i+1
}
DF<-as.data.frame(DF)

dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- ScaleData(object=dge,features=rownames(dge))
cluster.averages <- AverageExpression(dge,assays = "RNA",return.seurat = F)
mtx<-cluster.averages$RNA

Type<-Feature$type[order(Feature$type)]
Names<-Feature$names[order(Feature$type)]
Feature$names=as.vector(Names)
Feature$type=as.vector(Type)
Feature<-as.data.frame(Feature)

mtx<-cluster.averages$RNA
mtx<-mtx[rownames(mtx) %in% Feature$names,]
mtx<-log2(mtx+1)
mtx<-mtx[Feature$names,]

matrix<-dge_organoid@assays$RNA@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod[c("ACE2"),])
gene<-as.numeric(dge_organoid@meta.data$SARS_cov1)
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations<-sort(correlations,decreasing = T)
correlation_ACE2<-correlations
#correlation_TMPRSS2<-correlations
correlation_SARS_cov<-correlations
write.table(correlation_SARS_cov,"correlation_organoid.txt",col.names = T,row.names = T,sep = "\t")


library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


my_levels <- c("B-cell","T-cell","Mast","Plasma","Fibroblast","Smooth_muscle","Endothelial","Myeloid","Epithelial")

# Re-level object@ident
dge_combine@active.ident <- factor(x = dge_combine@active.ident, levels = my_levels)


features<- c("MS4A1","CD22","CD79A","CD8A","CD3D","IL7R","CPA3","KIT","TPSAB1","IGKC","MZB1","IGJ","C1S",
             "C7","DCN","MYH11","ACTA2","RGS5","IFI27","VWF","SELE","APOE","LYZ","IL1B","EPCAM","KLK3","KRT5")

StackedVlnPlot(obj = dge, features = features)

StackedVlnPlot(obj = dge_temp, features = c("KRT5","KRT15","KLK3","ACPP",
                                            "NKX3-1","PIGR","CP","ERG",
                                            "PCA3","SPON2"))

plot1 <- FeatureScatter(dge_pca, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "all")
plot2 <- FeatureScatter(dge_pca, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "all")
plot1 + plot2
ggsave(file="Scatter.eps",width = 20,height = 20,units = "cm")

dge<-SubsetData(dge_temp,ident.use = c("BE","BE_Normal","Hillock_Normal"))
dge<-SubsetData(dge_temp,ident.use = c("Club","Club_Normal","Hillock_Normal"))
dge<-SubsetData(dge_temp,ident.use = c("LE","LE_Normal","Club_Normal","Club"))

dge<-SubsetData(E_integrated,ident.use = c("Club","ERGpos_Tumor","ERGneg_Tumor"))
dge<-LE_PCA
dge$pair_type="others"
dge$pair_type[dge$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")]="PairN"
dge$pair_type[dge$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T")]="PairT"
DefaultAssay(dge)<-"RNA"
dge<-SetIdent(dge,value = as.vector(dge$pair_type))
dge<-SubsetData(dge,ident.remove = c("others"))
V1<-as.vector(FetchData(object = dge,vars = "AR"))
BOX_df<-NULL
BOX_df$id<-dge@active.ident
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
  geom_signif(comparisons = list(c("PairN","PairT")),map_signif_level=TRUE,y_position = max(V1)+0.15)+
  ggtitle("AR_expression")

dge<-SetIdent(dge,value = as.vector(dge@active.ident))
dge_1<-SubsetData(dge,ident.use = c("ERGneg_Tumor","ERGpos_Tumor"))
dge_1<-SetIdent(dge_1,value = as.vector(dge_1$ID))

pbmc.markers <- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.001, logfc.threshold = 0.001)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge_1_pca@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)
pbmc.markers<-pbmc.markers[!(rownames(pbmc.markers) %in% remove_gene),]
#dge_subset is 12 idents, 50 cells each.
top50_BE_Pair <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top50_Club_Tumor <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top20_Tumor <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

dge_temp<-SubsetData(dge,max.cells.per.ident = 100)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
heatmap_gene<-top50_Club_Pair$gene[top50_Club_Pair$p_val_adj<=0.1]
DoHeatmap(dge_temp, features = c(heatmap_gene),raster = F) + theme(axis.text.y = element_text(size = 6))#+ scale_fill_gradientn(colors = c("white", "red"))
ggsave(file="./Heatmap_Club_Pair_purple.pdf",width = 20,height = 30,units = "cm")
write.table(top50_Club_Pair,"Club_Pair_top50_DEG.txt",col.names = T,row.names = T)

levels(x = dge_temp) <- c("ERGneg_Tumor","ERGpos_Tumor","ERGneg_Stromal","ERGpos_Stromal","ERGneg_CD4","ERGpos_CD4")
levels(x = dge_temp)<-c("BE","BE_Normal","Club","Club_Normal","LE","LE_Normal","Hillock_Normal","ERGpos_Tumor","ERGneg_Tumor")
DoHeatmap(dge_temp, features = c(top20_CD4$gene,top20_Stromal$gene,top20_Tumor$gene),raster = F) + theme(axis.text.y = element_text(size = 5))#+ scale_fill_gradientn(colors = c("white", "red"))
DoHeatmap(dge_temp, features = c(top50_Club_Tumor$gene),raster = F) + theme(axis.text.y = element_text(size = 6))#+ scale_fill_gradientn(colors = c("white", "red"))
DoHeatmap(dge_temp, features = unique(c(top50_BE$gene,top50_Club$gene,top50_LE$gene)),raster = F) + theme(axis.text.y = element_text(size = 6))#+ scale_fill_gradientn(colors = c("white", "red"))

ggsave(file="./Heatmap_Club_Tumor_purple.pdf",width = 20,height = 30,units = "cm")
ggsave(file="./Heatmap_ERGstatus_purple_byERGtype_top20.eps",width = 20,height = 25,units = "cm")
write.table(top50_Club_Tumor,"Club_Tumor_top50_DEG.txt",col.names = T,row.names = T)

cluster.averages_club <- AverageExpression(dge,assays = "RNA",return.seurat = T)
genelist_6<-genelist
genelist_5<-genelist_6[1:50]
genelist_temp<-Response_features
DoHeatmap(cluster.averages_organoid,assay = "RNA",slot = "data", features = unlist(genelist_temp),raster = F, size = 3,draw.lines = FALSE,lines.width = 0.5)+ 
  scale_fill_gradientn(colors = c("blue", "red"))+ theme(axis.text.y = element_text(size = 4))
DoHeatmap(cluster.averages_club,assay = "RNA",features = unlist(top10$gene),raster = F, size = 3,draw.lines = FALSE,lines.width = 0.5)+ 
  scale_fill_gradientn(colors = c("white", "red"))+ theme(axis.text.y = element_text(size = 4))
ggsave(file="Heatmap_10_whitered.pdf",width = 20,height = 10,units = "cm",limitsize = FALSE)


dge<-SubsetData(dge_temp,max.cells.per.ident = 100)

DoHeatmap(dge,assay = "RNA",slot = "scale.data", features = c(Response_features),raster = F, size = 3,draw.lines = FALSE)#+ scale_fill_gradientn(colors = c("blue", "red"))
FeaturePlot(dge_temp,features=c("AR_signature1"))
ggsave(file="./heatmap_NE.pdf",width = 10,height = 15,units = "cm")

set_BE<-dot_df_BE_ERGneg[dot_df_BE_ERGneg$FDR.q.val<0.1,]
set_LE<-dot_df_LE_ERGneg[dot_df_LE_ERGneg$FDR.q.val<0.1,]
set_Club<-dot_df_Club_ERGneg[dot_df_Club_ERGneg$FDR.q.val<0.1,]



dge_temp<-SetIdent(dge_E,value = as.vector(dge_E$seurat_clusters))
levels(dge_temp)=c("0","1","2","3","4","5","6",
                   "7","8","9","10","11","12","13",
                   "14","15","16","17","18","19")
StackedVlnPlot(dge_temp,features=c("KRT5","KRT15","KLK3","ACPP","NKX3-1","PIGR","MMP7","CP","ERG","PCA3","SPON2"))
ggsave(file="./E_supervised_stacked.pdf",width = 25,height = 15,units = "cm")

VlnPlot(dge_temp,features=c("BE1"),pt.size = 0)
ggsave(file="./E_BE_signature.pdf",width = 20,height = 10,units = "cm")

VlnPlot(dge_temp,features=c("Luminal1"),pt.size = 0)
ggsave(file="./E_LE_signature.pdf",width = 20,height = 10,units = "cm")

ggsave(file="./E_25clusters.eps",width = 20,height =20,units = "cm")


dge<-Club_integrated
dge<-SetIdent(dge,value = as.vector(dge$ID))
dge<-SubsetData(dge_E,ident.use = c("LE","Club"))
dge<-SubsetData(dge_E,ident.use = c("BE","LE"))
dge$refinedID="LE"
dge$refinedID[colnames(dge) %in% cell_club_0]="Club_0"
dge$refinedID[colnames(dge) %in% cell_club_others]="Club_Others"
DefaultAssay(dge)="RNA"
dge$refinedID="LE"
dge$refinedID[colnames(dge) %in% cell_BE_6]="BE_6"
dge$refinedID[colnames(dge) %in% cell_BE_others]="BE_Others"
DefaultAssay(dge)="RNA"
dge<-SetIdent(dge,value = as.vector(dge$refinedID))
StackedVlnPlot(dge,features = c("SCGB3A1","LCN2","S100A6","EEF1A1","NEAT1","SAT1","B2M"))
StackedVlnPlot(dge,features = c("AR","KLK3","KLK2","NKX3-1","ACPP","CP","LTF"))
StackedVlnPlot(dge,features = c("KLK2","KLK3","ACPP","NKX3-1","KRT5","KRT15","TP63"))
StackedVlnPlot(dge,features = c("KLK2","KLK3","ACPP","NKX3-1","PIGR","MMP7","CP","SCGB1A1"))
ggsave(file="./BE_stacked.pdf",width = 10,height =25,units = "cm")


ggsave(file="./Club_normal_PCa_up.pdf",width = 10,height =25,units = "cm")
dge_temp<-SubsetData(dge_E,ident.remove = c("NE"))
StackedVlnPlot(dge_temp,features = c("KRT8","KRT18","SPINK1","NKX3-1","TACSTD2"))
ggsave(file="../luminal3_NG_check_2.pdf",width = 15,height =10,units = "cm")

StackedVlnPlot(dge_Henry,features=c("KRT8","KRT18","NKX3-1","TACSTD2","PSCA","PIGR","SPINK1","KRT4"))
ggsave(file="../luminal3_NG_check_1.pdf",width = 15,height =10,units = "cm")
cell_PIGR<-colnames(SubsetData(dge_tmp,subset.name = "PSCA",low.threshold = 0))
dge_temp$ID[colnames(dge_temp) %in% cell_PIGR]="Club_organoid"
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$ID))

dge<-SetIdent(dge,value = as.vector(dge$seurat_clusters))
my_levels <- c("0", "1", "2", "3", "4", "5", "6","7","8","9","10","11","12","13","14","15","16","17","18","19")

# Relevel object@timepoint
dge@active.ident<-factor(dge@active.ident,levels=my_levels)
StackedVlnPlot(dge,features=c("SPON2","TRPM8","SRGN","RGS1","PCAT4","TRGC1","TFF3","AZGP1","NCAPD3"))
dge<-SetIdent(dge,value = as.vector(dge$ID))
all_markers<-FindMarkers(object = dge,ident.1 = "")

library(SoupX)
library(Matrix)
rowSums = Matrix::rowSums
colSums = Matrix::colSums
tod=as.matrix(club_PCA@assays$RNA@counts)
toc=as.matrix(club_PCA@assays$RNA@counts)
sc = SoupChannel(tod = as.matrix(club_PCA@assays$RNA@counts),toc = as.matrix(club_PCA@assays$RNA@counts), calcSoupProfile = FALSE)
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), 
                      counts = rowSums(toc))
sc = setSoupProfile(sc, soupProf)




