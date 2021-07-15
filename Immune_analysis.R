setwd("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/")
library(Seurat)
library(cowplot)
library(dplyr)
library(icesTAF)
library(GEOquery)
library(ggplot2)
library(data.table)
library(cluster)
library(factoextra)
library(reshape2)
theme_set(theme_cowplot())

dge<-readRDS("../dge_Myeloid.rds")
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 500)
top20 <- head(x = VariableFeatures(object = dge), 20)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = dge)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE,xnudge = 0,ynudge = 0)
CombinePlots(plots = list(plot1, plot2))
ggsave(file="VariableFeature.eps",width = 40,height = 40,units = "cm")
all.genes <- rownames(x = dge)
dge <- ScaleData(object=dge,features=VariableFeatures(object = dge))
dge <- ScaleData(object=dge,features=rownames(dge))
dge <- RunPCA(object = dge, features = VariableFeatures(object = dge),npcs = 50)
VizDimLoadings(object = dge, dims = 1:20, reduction = "pca")
ggsave(file="PCA1_20.eps",width = 100,height = 100,units = "cm")
DimPlot(object = dge, reduction = "pca")
PCAPlot(object=dge,dim.1=1,dim.2=2)
ggsave(file="PCA.eps",width = 30,height = 30,units = "cm")

slot(dge[["pca"]], "misc")
print(x = dge[["pca"]], dims = 1:5, nfeatures = 5)
mat <- Seurat::GetAssayData(dge,assay="RNA",slot="scale.data")
pca <-dge[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)


png("PCA_heatmap.png")
DimHeatmap(object = dge, dims = 1:20, cells = length(dge$nCount_RNA), balanced = TRUE)
dev.off()
png("Elbowplot_dge.png")
ElbowPlot(dge,ndims = 50)
dev.off()

dge <- JackStraw(dge, num.replicate = 100,dims = 50)
dge <- ScoreJackStraw(dge, dims = 1:50)
JackStrawPlot(dge, dims = 1:50)
ggsave(file="Strawplot.eps",width = 20,height = 20,units = "cm")
ElbowPlot(dge)
ggsave(file="Elbow.eps",width = 10,height = 10,units = "cm")


#Myeloid
n_pc=27
resolut=1.5
dge_precluster <- dge
dge<-dge_precluster
dge <- FindNeighbors(dge, dims = 1:n_pc,k.param = 50)
dge <- FindClusters(dge, resolution = resolut)

dge=RunTSNE(dge,dims.use = 1:n_pc,perplexity=40,seed.use = 10)
TSNEPlot(dge,pt.size = 1,label=T)
ggsave(file="TSNE.eps",width = 20,height = 20,units = "cm")
dge <- RunUMAP(dge, dims = 1:n_pc)
DimPlot(dge, reduction = "umap",label=T)
ggsave(file="Umap.eps",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "tsne",label=T,group.by = "orig.ident")
ggsave(file="TSNE_group.eps",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "umap",label=T,group.by = "orig.ident")
ggsave(file="Umap_group.eps",width = 20,height = 20,units = "cm")

Tab<-table(dge$orig.ident,Idents(object = dge))
write.table(Tab,"table_IDclustering.txt",sep="\t",row.names = T,col.names = T)

dge.markers <- FindAllMarkers(dge, only.pos = TRUE, logfc.threshold = 0.5)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
dge.markers<-dge.markers[!(dge.markers$gene %in% remove_gene),]
write.table(dge.markers, file = "Myeloid_markers_bycluster.txt",sep="\t",col.names=T)
write.table(dge.markers, file = "Myeloid_markers_byID.txt",sep="\t",col.names=T)
write.table(dge.markers, file = "Myeloid_markers_bypatient.txt",sep="\t",col.names=T)
top10 <- dge.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
View(top10)

Myeloid_features <- unique(as.vector(top10$gene))
DotPlot(dge, features =Myeloid_features) + RotatedAxis()
ggsave(file="dot_plot_bycluster.eps",width = 70,height = 50,units = "cm")
ggsave(file="dot_plot_bypatient.eps",width = 70,height = 50,units = "cm")
ggsave(file="dot_plot_byID.eps",width = 50,height = 30,units = "cm")


dge <- ScaleData(object=dge,features=rownames(dge))
DoHeatmap(dge, features = as.character(top10$gene),combine = T,raster = T) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_bycluster.eps",width = 75,height = 75,units = "cm")
ggsave(file="Heatmap_bypatient.eps",width = 75,height = 75,units = "cm")
ggsave(file="Heatmap_byID.eps",width = 75,height = 75,units = "cm")

VlnPlot(dge,features = c("CD3E","CD3D","CD3G","CD8A","CD8B","CTLA4","PDCD1","LAG3","TIGIT","CD14","CD68","CD163","FCER1A","KIT","FCER2","ENPP3","ENPP3","KIT","CPA3","LYZ","VCAN","KLK2","KLK3","ACPP"))
ggsave(file="Markers.eps",width = 75,height = 75,units = "cm")

new.cluster.ids <- c("Monocyte","Dendritic_cell","Monocyte","Macrophage","Monocyte")
names(new.cluster.ids) <- levels(dge)
dge <- RenameIdents(dge, new.cluster.ids)
DimPlot(dge, reduction = "umap",label=T)
ggsave(file="Umap_ID.eps",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "tsne",label=T)
ggsave(file="TSNE_ID.eps",width = 20,height = 20,units = "cm")

saveRDS(dge,"myeloid_dge.RDS")
dge<-readRDS("Immune_dge.RDS")



VlnPlot(object = dge, features = c("CD14","CD68","CD163"))
ggsave(file="Macrophage_markers_feature.eps",width = 40,height = 40,units = "cm")
VlnPlot(object = dge, features = c("CCL2","TGFB1","VEGFA","PCNA","EGF","CCR2"<"CSF1R","CCL22","CCL17"))
ggsave(file="TAM_markers_feature.eps",width = 40,height = 40,units = "cm")

cell_CD4<-colnames(SubsetData(dge,subset.name="CD4",low.threshold=0))
cell_FOXP3<-colnames(SubsetData(dge,subset.name="FOXP3",low.threshold=0))
cell_CD25<-colnames(SubsetData(dge,subset.name="IL2RA",low.threshold=0))

cell_Treg<-union(intersect(cell_CD4,cell_FOXP3),intersect(cell_CD25,cell_FOXP3))
cell_CD4T<-setdiff(colnames(SubsetData(dge,subset.name = "CD4",low.threshold = 0)),cell_Treg)
cell_CD8T<-union(setdiff(colnames(SubsetData(dge,subset.name = "CD8A",low.threshold = 0)),cell_Treg),setdiff(colnames(SubsetData(dge,subset.name = "CD8B",low.threshold = 0)),cell_Treg))

dge$newid="Other T-cell"
dge$newid[colnames(dge) %in% cell_Treg]="Treg"
dge$newid[colnames(dge) %in% cell_CD4T]="CD4+ T-cell"
dge$newid[colnames(dge) %in% cell_CD8T]="CD8+ T-cell"

DimPlot(dge,reduction = "tsne",group.by = "newid")
ggsave(file="TSNE_id.eps",width = 20,height = 20,units = "cm")
DimPlot(dge,reduction = "umap",group.by = "newid")
ggsave(file="UMAP_id.eps",width = 20,height = 20,units = "cm")


dge_immune<-dge
dge_Macrophage<-SubsetData(dge,ident.use = c("0","7","8","2","14","13","14","16","10"))
dge_Tcell<-SubsetData(dge,ident.use = c("1","3","4","5","6","11","12"))


dge<-dge_Tcell
dge<-dge_immune
tmp<-read.csv("../../geneset/Tcell.gmt",header=FALSE,sep="\t")
for (i in 1:nrow(tmp)) {
  
  name_feature=as.character(tmp[i,1])
  feature_tmp<-na.omit(tmp[i,])
  
  Features<-as.vector(feature_tmp[2:ncol(feature_tmp)])
  dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature)
  names(x = dge[[]])
  
  #str<-paste0(name_feature,"1")
  #dge[[str]]<-range01(dge[[str]])
  
  FeaturePlot(object = dge, features =paste0(name_feature,"1") ,cols =c("grey","blue","pink","red"),label = T)
  ggsave(file=paste0("Modulescore_",name_feature,".eps"),width = 20,height = 20,units = "cm")
  
  VlnPlot(object = dge, features =paste0(name_feature,"1"))
  ggsave(file=paste0("Modulescore_",name_feature,"_vlnplot.eps"),width = 20,height = 20,units = "cm")
  
  
}

FeaturePlot(dge,features=c("CD4","CD3D","CD3E","CD3G","PDCD1","CD274"))
ggsave(file="Misc_markers_feature.eps",width = 40,height = 40,units = "cm")


dge_backup<-dge
new.cluster.ids <- c("Macrophage","T-cell","Macrophage","T-cell","T-cell","T-cell","T-cell","Macrophage","Macrophage","Mast","Leu","NK-cell","Exhausted_Tcell","Macrophage","Macrophage","Plasma","Macrophage")
names(new.cluster.ids) <- levels(dge)
dge <- RenameIdents(dge, new.cluster.ids)
DimPlot(dge, reduction = "umap",label=T)
ggsave(file="Immune_Umap_ID.eps",width = 20,height = 20,units = "cm")


dge_subset<-SubsetData(dge,max.cells.per.ident = 65)
pbmc.markers <- FindAllMarkers(dge_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(dge_subset, features = top10$gene,combine = T) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_subset.eps",width = 50,height = 50,units = "cm")



dge<-dge_combined
dge_backup<-dge
dge<-SubsetData(dge,ident.use = c("T-cell","NK-cell","Mast","Macrophage","Leu","Exhausted_Tcell"))
dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 2000)
dge<-subset(x=dge,features=VariableFeatures(object = dge))
#Correlation of average expression in each cluster
av.exp <- AverageExpression(dge,)$RNA
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
#
cor.df <- tidyr::gather(data = cor.exp, y, correlation, c("T-cell","Macrophage","Mast","NK-cell","Exhausted_Tcell","Leu"))

#cor.df <- tidyr::gather(data = cor.exp, y, correlation, c("T-cell","Macrophage","Mast","NK-cell","Plasma","Exhausted_Tcell","Leu"))
#cor.df <- tidyr::gather(data = cor.exp, y, correlation, c("T-cell","Macrophage","LE","Monocyte","Mast","NK","MZB","Cell-cycle"))
ggplot(cor.df, aes(x, y),lab_size = 5) +geom_tile(aes(fill = correlation))+scale_fill_gradientn(colours = c("blue","cyan","lightcyan","red"))
ggsave(file="Cluster_correlation_Immune.eps",width = 20,height = 20,units = "cm")

#Plot correlation in each cluster
dge_Tcell=subset(x=dge,ident = c("T-cell"))
dge_Macro=subset(x=dge,ident = c("Macrophage"))
dge_NK=subset(x=dge,ident = c("NK"))
dge_MZB=subset(x=dge,ident = c("Plasma"))
dge_Cellcycle=subset(x=dge,ident = c("Cell-cycle"))

dge_temp<-dge
av.exp <- Seurat::GetAssayData(dge_temp,assay="RNA")#,slot="scale.data")
dge_temp <- NormalizeData(object = dge_temp, normalization.method = "LogNormalize", scale.factor = 10000)
dge_temp <- FindVariableFeatures(object = dge_temp, selection.method = "vst", nfeatures = 200)
varied_list = c(VariableFeatures(object = dge_temp))
av.exp<-av.exp[rownames(av.exp) %in% varied_list,]
cell_list<-sprintf("Cell_%d",seq(1:ncol(av.exp)))
colnames(av.exp)<-cell_list
#av<-t(av.exp)
av<-av.exp
res.dist <- get_dist(av.exp, stand = T, method = "pearson")
fviz_dist(res.dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"),lab_size = 5)
ggsave(file=paste0(string,"_gene_Correlation_cluster.pdf"),width = 50,height = 50,units = "cm")


av.exp <- Seurat::GetAssayData(dge_temp,assay="RNA",slot="scale.data")
cell_list<-sprintf("Cell_%d",seq(1:ncol(av.exp)))
colnames(av.exp)<-cell_list
#av<-t(av.exp)
av<-av.exp
res.dist <- get_dist(t(av.exp), stand = F, method = "pearson")
fviz_dist(res.dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"),lab_size = 3)
ggsave(file=paste0(string,"_cell_Correlation_cluster.pdf"),width = 50,height = 50,units = "cm")

cell_immune<-colnames(dge)


# Try module score


#Macrophage analysis
dge<-readRDS("Immune_dge.RDS")
dge_backup<-dge
dge<-dge_Macro

# TAM markers
FeaturePlot(object = dge, features = c("CCL2","TGFB1","VEGFA","PCNA","EGF","CCR2","CSF1R","CCL22","CCL17"))
ggsave(file="TAM_markers.eps",width = 40,height = 40,units = "cm")

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "IFNG", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_IFNG <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "TNFRSF1A", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_TNFRSF1A <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "CXCL10", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CXCL10 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "IL12A", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_IL12A <- c(cell_names, genes[! genes %chin% cell_names])

cell_M1Macro<-unique(c(intersect(cell_CXCL10,cell_TNFRSF1A),intersect(cell_TNFRSF1A,cell_IFNG),intersect(cell_TNFRSF1A,cell_IL12A),intersect(cell_CXCL10,cell_IL12A),intersect(cell_CXCL10,cell_IFNG),intersect(cell_IFNG,cell_IL12A)))
dge_M1Macro<-SubsetData(dge,cells=cell_M1Macro)
#M2 markers: "MMP2","VTCN1","STAT3","CD163","APOE"
VlnPlot(dge,features=c("MMP2","VTCN1","STAT3","CD163","APOE"))
seurat_subset <- SubsetData(dge, subset.name = "MMP2", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_MMP2 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "VTCN1", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_VTCN1 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "STAT3", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_STAT3 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "CD163", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CD163 <- c(cell_names, genes[! genes %chin% cell_names])

cell_M2Macro<-unique(c(intersect(cell_VTCN1,cell_MMP2),intersect(cell_VTCN1,cell_STAT3),intersect(cell_VTCN1,cell_CD163),intersect(cell_MMP2,cell_STAT3),intersect(cell_MMP2,cell_CD163),intersect(cell_STAT3,cell_CD163)))
dge_M2Macro<-SubsetData(dge,cells=cell_M2Macro)

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "CCL2", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CCL2 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "VEGFA", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_VEGFA <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "TGFB1", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_TGFB1 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "CCR2", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CCR2 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "PCNA", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_PCNA <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "CSF1R", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CSF1R <- c(cell_names, genes[! genes %chin% cell_names])

cell_TAM_1<-unique(c(intersect(cell_CCL2,cell_CCR2),intersect(cell_CCR2,cell_VEGFA),intersect(cell_CCL2,cell_VEGFA),intersect(cell_VEGFA,cell_PCNA),intersect(cell_TGFB1,cell_PCNA),intersect(cell_VEGFA,cell_TGFB1)))
#cell_TAM_2<-unique(c(intersect(cell_CCL2,cell_CSF1R),intersect(cell_CCR2,cell_CSF1R),intersect(cell_CSF1R,cell_VEGFA),intersect(cell_CSF1R,cell_PCNA),intersect(cell_TGFB1,cell_CSF1R)))
cell_TAM<-unique(cell_TAM_1)
dge_TAM<-SubsetData(dge,cell=cell_TAM)
saveRDS(dge_TAM,"dge_TAM.RDS")
table(dge_TAM$orig.ident)

name_feature="TAM_Markers"
Features<-as.vector(c("CCL2","TGFB1","VEGFA","PCNA","EGF","CCR2","CSF1R","CCL22","CCL17"))
dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature)
names(x = dge[[]])
dge@active.ident<-as.character(dge$orig.ident)
FeaturePlot(object = dge, features =paste0(name_feature,"1") )
ggsave(file=paste0("Modulescore_",name_feature,".eps"),width = 20,height = 20,units = "cm")

dge$TAM_Markers1<-range01(dge$TAM_Markers1)
score<-dge$TAM_Markers1

# LE<-Score[dge@active.ident=="LE"]
# BE<-Score[dge@active.ident=="BE"]
# OE<-Score[dge@active.ident=="OE"]
# Tumor_LE<-Score[dge@active.ident=="Tumor_LE"]
# VlnPlot(object = dge, features =paste0(name_feature,"1") ,cols =c("grey","blue","red","pink") )
# ggsave(file=paste0("Modulescore_",name_feature,"_vlnplot.eps"),width = 20,height = 20,units = "cm")
VlnPlot(object = dge, features =paste0(name_feature,"1"),group.by = "orig.ident") #,cols =c("green","blue","red","cyan","purple","orange") )
ggsave(file=paste0("Modulescore_",name_feature,"_vlnplot.eps"),width = 20,height = 20,units = "cm")
dge_E<-dge
saveRDS(dge_E,"E_subcluster_Nov.rds")


#T-cell analysis
dge<-readRDS("Immune_dge.RDS")
dge_backup<-dge
  dge<-dge_Tcell
dge<-dge_backup


dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 500)
dge <- ScaleData(object=dge,features=VariableFeatures(object = dge))
dge <- RunPCA(object = dge, features = VariableFeatures(object = dge),npcs = 50)
VizDimLoadings(object = dge, dims = 1:10, reduction = "pca")
ggsave(file="Tcell_PCA1_10.eps",width = 100,height = 100,units = "cm")
DimPlot(object = dge, reduction = "pca")
PCAPlot(object=dge,dim.1=1,dim.2=2)
ggsave(file="Tcell_PCA.eps",width = 30,height = 30,units = "cm")
png("Tcell_PCA_heatmap.png")
DimHeatmap(object = dge, dims = 1:20, cells = length(dge$nCount_RNA), balanced = TRUE)
dev.off()
png("Elbowplot_dge.png")
ElbowPlot(dge,ndims = 50)
dev.off()
dge <- JackStraw(dge, num.replicate = 100,dims = 50)
dge <- ScoreJackStraw(dge, dims = 1:50)
JackStrawPlot(dge, dims = 1:50)
ggsave(file="Tcell_Strawplot.eps",width = 20,height = 20,units = "cm")
n_pc=50
resolut=1
dge_precluster <- dge
dge<-dge_precluster
dge <- FindNeighbors(dge, dims = 1:n_pc,k.param = 50)
dge <- FindClusters(dge, resolution = resolut)
dge <- RunUMAP(dge, dims = 1:50)
DimPlot(dge, reduction = "umap",label=T)
ggsave(file="Immune_Umap_rawid.eps",width = 20,height = 20,units = "cm")
dge=RunTSNE(dge,dims.use = 1:n_pc,perplexity=40,seed.use = 10)
TSNEPlot(dge,pt.size = 1,label=T)
ggsave(file="TSNE_rawid.eps",width = 20,height = 20,units = "cm")
dge.markers <- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(dge.markers, file = "Immune_dge_markers.txt",sep="\t",col.names=T)
top50 <- dge.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
View(top50)
saveRDS(dge,"Tcell_dge.RDS")
dge<-readRDS("Tcell_dge.RDS")
dge_cor <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 200)
dge_cor<-subset(x=dge_cor,features=VariableFeatures(object = dge))
av.exp <- AverageExpression(dge_cor)$RNA
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, c("0","1","2","3","4"))
ggplot(cor.df, aes(x, y),lab_size = 5) +geom_tile(aes(fill = correlation))+scale_fill_gradientn(colours = c("blue","cyan","lightcyan","red"))
ggsave(file="Tcell_correlation_200feature.eps",width = 20,height = 20,units = "cm")
DoHeatmap(dge, features = c("CD8A","CD8B","CTLA4","PDCD1","LAG3","TIGIT","CD14","CD68","CD163")) + theme(axis.text.y = element_text(size = 10))
#DoHeatmap(dge, features = top15$gene) + theme(axis.text.y = element_text(size = 10))
ggsave(file="Tcell_Heatmap.eps",width = 50,height = 50,units = "cm")

# CD8+  markers
FeaturePlot(object = dge, features = c("CD8A","CD8B","FOXP3"))
ggsave(file="CD8+_markers.eps",width = 40,height = 40,units = "cm")
VlnPlot(object = dge, features = c("CD8A","CD8B"))
ggsave(file="CD8+_markers_violin.eps",width = 40,height = 40,units = "cm")

# Treg  markers
FeaturePlot(object = dge, features = c("CD4","IL7R","IL2RA"))
ggsave(file="Treg_markers.eps",width = 40,height = 40,units = "cm")
VlnPlot(object = dge, features = c("CD4","IL7R","IL2RA"))
ggsave(file="Treg_markers_violin.eps",width = 40,height = 40,units = "cm")

# ExhaustT markers
FeaturePlot(object = dge, features = c("CTLA4","PDCD1","LAG3","TIGIT"))
ggsave(file="Exhaust_markers.eps",width = 40,height = 40,units = "cm")
VlnPlot(object = dge, features = c("CTLA4","PDCD1","LAG3","TIGIT"))
ggsave(file="Exhaust_markers_violin.eps",width = 40,height = 40,units = "cm")

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "CTLA4", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CTLA4 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "PDCD1", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_PDCD1 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "LAG3", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_LAG3 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "TIGIT", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_TIGIT <- c(cell_names, genes[! genes %chin% cell_names])


cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "CD8A", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CD8A <- c(cell_names, genes[! genes %chin% cell_names])
cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge, subset.name = "CD8B", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CD8B <- c(cell_names, genes[! genes %chin% cell_names])

cell_CD8<-unique(intersect(c(cell_CD8A,cell_CD8B),colnames(dge_Tcell)))
dge_CD8<-SubsetData(dge,cell=cell_CD8)
saveRDS(dge_CD8,"dge_CD8.RDS")
table(dge_CD8$orig.ident)

cell_exhaust<-unique(c(intersect(cell_CTLA4,cell_PDCD1),intersect(cell_CTLA4,cell_LAG3),intersect(cell_PDCD1,cell_LAG3),intersect(cell_TIGIT,cell_LAG3),intersect(cell_TIGIT,cell_CTLA4),intersect(cell_TIGIT,cell_PDCD1)))
cell_exhaust<-intersect(cell_exhaust,colnames(dge_Tcell))
cell_exhaust8<-intersect(cell_exhaust,cell_CD8)
dge_exhaust<-SubsetData(dge,cell=cell_exhaust)
saveRDS(dge_exhaust,"dge_exhaust.RDS")
table(dge_exhaust$orig.ident)

dge_Tcell<-SubsetData(dge,ident.use = "T-cell")
DimPlot(object=dge,cells.highlight = cell_exhaust,reduction = "umap",label = F, pt.size = 0.5)
ggsave(file="Exhaust_T-cell.eps")
DimPlot(object=dge,cells.highlight = cell_CD8,reduction = "umap",label = F, pt.size = 0.5)
ggsave(file="CD8T-cell.eps")
DimPlot(object=dge,cells.highlight = intersect(colnames(dge_Macro),colnames(dge_TAM)),reduction = "umap",label = F, pt.size = 0.5)
ggsave(file="TAM-cell.eps")

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge_Tcell, subset.name = "CD4", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CD4 <- c(cell_names, genes[! genes %chin% cell_names])
dge_CD4<-SubsetData(dge_Tcell,cells = cell_CD4)


#Treg: "CD4","IL7R","IL2RA"
cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge_Tcell, subset.name = "CD4", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_CD4 <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge_Tcell, subset.name = "IL7R", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_IL7R <- c(cell_names, genes[! genes %chin% cell_names])

cell_names <- vector(mode="character")
seurat_subset <- SubsetData(dge_Tcell, subset.name = "IL2RA", low.threshold  = 0)
genes <- colnames(seurat_subset)
cell_IL2RA <- c(cell_names, genes[! genes %chin% cell_names])

cell_Treg=intersect(intersect(cell_IL2RA,cell_IL7R),cell_CD4)
dge_Treg<-SubsetData(dge_Tcell, cell=cell_Treg, low.threshold  = 0)

dge_backup<-dge
cell_ERGpos<-colnames(dge_id)[dge_id$orig.ident %in% c("PB1A","PB1B","AUG_PB1A","AUG_PB1B","PA_PR5199")]
cell_ERGneg<-colnames(dge_id)[dge_id$orig.ident %in% c("PB2A","PB2B","PA_PR5186","PA_PR5196")]
dge_id$ID<-dge_id@active.ident

dge_id<-SetIdent(object = dge_id,cells = cell_ERGpos,value = "ERGpos")
dge_id<-SetIdent(object = dge_id,cells = cell_ERGneg,value = "ERGneg")
dge_id$ERG<-dge_id@active.ident

dge_id@active.ident=dge_id$ID

dge<-SubsetData(dge_id,ident.use=c("T-cell"))

dge<-SubsetData(dge_id,ident.use = c("Fib"))
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 5000)
top20 <- head(x = VariableFeatures(object = dge), 20)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = dge)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE,xnudge = 0,ynudge = 0)
CombinePlots(plots = list(plot1, plot2))
ggsave(file="VariableFeature.eps",width = 40,height = 40,units = "cm")
all.genes <- rownames(x = dge)
dge <- ScaleData(object=dge,features=VariableFeatures(object = dge))
dge <- ScaleData(object=dge,features=all.genes)
dge <- RunPCA(object = dge, features = VariableFeatures(object = dge),npcs = 50)
VizDimLoadings(object = dge, dims = 1:10, reduction = "pca")
ggsave(file="PCA1_10.eps",width = 100,height = 100,units = "cm")
DimPlot(object = dge, reduction = "pca")
PCAPlot(object=dge,dim.1=1,dim.2=2)
ggsave(file="PCA.eps",width = 30,height = 30,units = "cm")
png("PCA_heatmap.png")
DimHeatmap(object = dge, dims = 1:20, cells = length(dge$nCount_RNA), balanced = TRUE)
dev.off()
png("Elbowplot_dge.png")
ElbowPlot(dge,ndims = 50)
dev.off()
dge$oldid<-dge@active.ident
n_pc=50
resolut=1.2
dge_precluster <- dge
dge<-dge_precluster
dge <- FindNeighbors(dge, dims = 1:n_pc,k.param = 50)
dge <- FindClusters(dge, resolution = resolut)

dge=RunTSNE(dge,dims.use = 1:n_pc,perplexity=40,seed.use = 10)
TSNEPlot(dge,pt.size = 1,label=T)
ggsave(file="TSNE.eps",width = 20,height = 20,units = "cm")
dge <- RunUMAP(dge, dims = 1:n_pc)
DimPlot(dge, reduction = "umap",label=T)
ggsave(file="Umap.eps",width = 20,height = 20,units = "cm")
dge.markers <- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.05)
write.table(dge.markers, file = "dge_markers.txt",sep="\t",col.names=T)
top15 <- dge.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
View(top15)

#DoHeatmap(dge, features = c("CD3E","CD3D","CD3G","CD8A","CD8B","CTLA4","PDCD1","LAG3","TIGIT","CD14","CD68","CD163","FCER1A","KIT","FCER2","ENPP3","ENPP3","KIT","CPA3","LYZ","VCAN","KLK2","KLK3","ACPP")) + theme(axis.text.y = element_text(size = 10))
DoHeatmap(dge, features = top15$gene,combine = T,assay = "RNA",slot="scale.data") + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap.eps",width = 75,height = 75,units = "cm")

Exhaust_Features<- as.vector(read.csv("../../geneset/PD1.gmt",header=FALSE,sep="\n")$V1)
#FeaturePlot(object = dge, features = c("CCL2","TGFB1","VEGFA","PCNA","EGF","CCR2"<"CSF1R","CCL22","CCL17"),reduction = "tsne")
#ggsave(file="TAM_markers_feature.eps",width = 40,height = 40,units = "cm")
PD1_features=c("CD247",	"CD274",	"CD3D","CD3E",	"CD3G",	"CD4",	"CSK",	"HLA-DPA1",	"HLA-DPB1",	"HLA-DQA1",	"HLA-DQA2", 	"HLA-DRB1",	"HLA-DRB3",	"HLA-DRB5",	"LCK",	"PDCD1",	"PDCD1LG2",	"PTPN6")

M1_macro_features=c("CCR7","IL2RA","IL15RA","IL7R","CXCL11","CCL19","CXCL10","CXCL9",
                    "TNF","CCL5","IL12B","IL15","IL6","CCL20","BCL2A1","RAS","BIRC3",
                    "GADD45G","SLC7A5","SLC2A6","SLC31A2","PLA1A","OASL","CHI3L2","HSD11B1",
                    "AK3","SPHK1","PFKFB3","PSME2","PFKP","PSMB9","PSMA2","OAS2")
M2_macro_features=c("TGFBR2","HRH1","TLR5","MSR1","CXCR4","P2RY14","MS4A6A",
                    "CD36","MS4A4A","MRC1","IGF1","CCL23","CCL18","SLC4A7","SLC38A6",
                    "CTSC","HEXB","LIPA","ADK","HNMT","TPST2","CERK","HS3ST2",
                    "LTA4H","CA2","ALOX15","HS3ST1","TGFBI","SEPP1","CHN2",
                    "FN1","FGL2","GAS7","EGR2","MAF")
TAM_features=c("CCL5","IL10","CSF1R","TGFB1","CCL2","CCL17","CCL22","IL6",
               "CXCL10","CXCL12","CD274","PDCD1LG","NOS2","MRC1","VTCN1","TGFB2")

dge <- AddModuleScore(object = dge, features = list(M2_macro_features), name = "M2")
names(x = dge[[]])
dge@active.ident<-as.character(dge$orig.ident)
FeaturePlot(object = dge, features ="M21" ,label = F,reduction = "umap")
ggsave(file=paste0("Modulescore_","M2_macro",".eps"),width = 20,height = 20,units = "cm")
VlnPlot(object = dge, features ="TAM1" )

ggsave(file="PD1_score.eps",width = 20,height = 20,units = "cm")
VlnPlot(object = dge, features ="PD11",group.by = "ERG" ,cols =c("gold1","blue1","green3","purple4","coral","darkcyan","gray4","pink4","aquamarine2","blueviolet","red2","gold4","darkolivegreen4","cyan","orange","orchid4","chocolate","royalblue","deeppink4","indianred") )

ggsave(file="Exhaust_ERG_score.eps",width = 20,height = 20,units = "cm")
ERGpos_CAF<-dge$Cytotoxic1[dge$ERG=="ERGpos"]
ERGneg_CAF<-dge$Cytotoxic1[dge$ERG=="ERGneg"]
wilcox.test(ERGpos_CAF,ERGneg_CAF,alternative = "two.sided")

dge$TAM_Markers1<-range01(dge$TAM_Markers1)
score<-dge$TAM_Markers1

cell_Exhaust<-colnames(dge)[dge$Exhaust1>=0.415]
dge_Exhaust<-SubsetData(dge,cells=cell_Exhaust)
cell_TAM<-colnames(dge)[dge$TAM1>0]
dge_TAM<-SubsetData(dge,cells = cell_TAM)
tmp<-read.csv("../../geneset/PD1.gmt",header=FALSE,sep="\t")

for (i in 1:nrow(tmp)) {
  
  name_feature=as.character(tmp[i,1])
  feature_tmp<-na.omit(tmp[i,])
  
  Features<-as.vector(feature_tmp[3:ncol(feature_tmp)])
  dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature)
  names(x = dge[[]])
  
  #str<-paste0(name_feature,"1")
  #dge[[str]]<-range01(dge[[str]])
  
  #FeaturePlot(object = dge, features =paste0(name_feature,"1") ,cols =c("grey","blue","red","pink"),label = T)
  #ggsave(file=paste0("Modulescore_",name_feature,".eps"),width = 20,height = 20,units = "cm")
  
  VlnPlot(object = dge,group.by = "ERG", features =paste0(name_feature,"1") ,cols =c("gold1","blue1","green3","purple4","coral","darkcyan","gray4","pink4","aquamarine2","blueviolet","red2","gold4","darkolivegreen4","cyan","orange","orchid4","chocolate","royalblue","deeppink4","indianred") )
  ggsave(file=paste0("Modulescore_",name_feature,"_vlnplot.eps"),width = 20,height = 20,units = "cm")
  
  
}

FeaturePlot(dge,features=c("COL1A2","FAP","PDPN","DCN","COL3A1","COL6A1","ACTA2"))

cell_Tcell<-colnames(SubsetData(dge_immune,ident.use = c("CD4+ T_cell","CD8+ T_cell","Treg","T_cell")))
cell_Macrophage<-colnames(SubsetData(dge_immune,ident.use = c("Macrophage","DC")))
dge_immune<-SetIdent(dge_immune,cells = cell_Tcell,value = "T-cell")
dge_immune<-SetIdent(dge_immune,cells = cell_Macrophage,value = "Macrophage")
immune_markers<- FindAllMarkers(dge_immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(dge.markers, file = "dge_markers.txt",sep="\t",col.names=T)
top10 <- immune_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
View(top15)
dge_subset<-SubsetData(dge_immune,max.cells.per.ident = 50)
all.genes <- rownames(x = dge_subset)

dge_subset <- ScaleData(object=dge_subset,features=all.genes)
DoHeatmap(dge_subset, features = top10$gene,combine = T,assay = "RNA",slot="scale.data") + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap.eps",width = 75,height = 75,units = "cm")


immune_markers<- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(immune_markers, file = "dge_markers.txt",sep="\t",col.names=T)
top10 <- immune_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
View(top15)
dge_subset<-SubsetData(dge,max.cells.per.ident = 50)
all.genes <- rownames(x = dge_subset)
dge_subset <- ScaleData(object=dge_subset,features=all.genes)
DoHeatmap(dge_subset, features = top10$gene,combine = T) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_0109.eps",width = 75,height = 75,units = "cm")



dge_backup<-dge
new.cluster.ids <- c("T-cell","Monocyte","T-cell","Monocyte","T-cell","T-cell","Macrophage","Monocyte","Mast","T-cell","Monocyte","11","Macrophage","NK-cell","Plasma","Neutrophil","Monocyte")
names(new.cluster.ids) <- levels(dge)
dge <- RenameIdents(dge, new.cluster.ids)
DimPlot(dge, reduction = "umap",label=T)
ggsave(file="Immune_Umap_ID_0110_label.eps",width = 20,height = 20,units = "cm",)

new.cluster.ids <- c("T-cell","Myeloid","T-cell","Myeloid","T-cell","Myeloid","T-cell","Mast","Myeloid","Myeloid","PTPRC+ Luminal","NK-cell","Plasma","Neutrophil","Myeloid")
names(new.cluster.ids) <- levels(dge)
dge <- RenameIdents(dge, new.cluster.ids)
DimPlot(dge,label = T)

immune_markers<- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.table(immune_markers, file = "myeloid_markers_id.txt",sep="\t",col.names=T)
top25 <- immune_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)
top5 <- immune_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
View(top15)
dge_subset<-SubsetData(dge,max.cells.per.ident = 50)
dge_subset <- FindVariableFeatures(object = dge_subset, selection.method = "vst", nfeatures = 5000)
dge_subset <- ScaleData(object=dge_subset,features=VariableFeatures(object = dge_subset))
DoHeatmap(dge_subset, features = top5$gene,combine = T,raster = F) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_myeloid_id.eps",width = 75,height = 75,units = "cm")

dge_T<-SubsetData(dge,ident.use = c("T-cell"))

Tcell_markers<- FindAllMarkers(dge_T, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.table(Tcell_markers, file = "immune_markers.txt",sep="\t",col.names=T)
top25 <- Tcell_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)
top5 <- Tcell_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
View(top15)
dge_T <- FindVariableFeatures(object = dge_T, selection.method = "vst", nfeatures = 5000)
dge_T <- ScaleData(object=dge_T,features=VariableFeatures(object = dge_T))
DoHeatmap(dge_T, features = top5$gene,combine = T,raster = F) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_Tcell_new.eps",width = 75,height = 75,units = "cm")

dge_myeloid<-SubsetData(dge,ident.use = "Myeloid")
dge_myeloid<-SetIdent(dge_myeloid,value = dge_myeloid$orig.ident)
#MB21D1 is CGAS
#TMEM173 is STING
VlnPlot(dge_myeloid,features=c("MB21D1","TMEM173"))
ggsave(file="Innate_inflam.eps",width = 20,height = 20,units = "cm")

dge$ID_2<-dge$ID
dge$ID_2[colnames(dge) %in% cell_M1 ]="M1_Macrophage"
dge$ID_2[colnames(dge) %in% cell_M2 ]="M2_Macrophage"
dge$ID_2[dge$ID_2=="Macrophage"]="M0_Macrophage"
dge<-SetIdent(dge,value = as.vector(dge$ID_2))
DimPlot(dge)
ggsave(file="Myeloid_refinedID.eps",width = 20,height = 20,units = "cm")
dge$ERGtype<-"ERGneg"
dge$ERGtype[dge$orig.ident %in% c("AUG_PB1","MAY_PB1","PR5249","PR5251","PR5199","PR5269")]="ERGpos"

DimPlot(dge,group.by = "ERGtype")
ggsave(file="Myeloid_byERGtype.eps",width = 20,height = 20,units = "cm")
table(dge$ERGtype,dge@active.ident)



#####Checking PD-1 genes:
Features=c("CD28","CD274","CD247","IFNG","CTLA4","INPP5D","INPPL1","CD58","CD27","CD70","HLA-A","CD74")
for (i in 1:length(Features)) {
V1<-as.vector(FetchData(object = dge,vars = Features[i]))
BOX_df<-NULL
BOX_df$id<-dge@active.ident
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
  geom_signif(comparisons = list(c("CD4_cluster2","CD4_cluster1")),map_signif_level=TRUE,y_position = max(V1)+0.15)+
  ggtitle(Features[i])
ggsave(file=paste0(Features[i],"_vlnplot.eps"),width = 20,height = 20,units = "cm")
}


StackedVlnPlot(obj = dge, features = Features,plot.margin = 2)
ggsave(file="CD4_PD1_correlation_stackedviolin.eps",width = 15,height = 20,units = "cm")


pbmc.markers <- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.001, logfc.threshold = 0.001)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge_1_pca@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)
pbmc.markers<-pbmc.markers[!(rownames(pbmc.markers) %in% remove_gene),]
#dge_subset is 12 idents, 50 cells each.
top50_BE_Pair <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top50_Club_Tumor <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top20_Tumor <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

dge_temp<-SubsetData(dge,max.cells.per.ident = 100)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
heatmap_gene<-top50_Club_Pair$gene[top50_Club_Pair$p_val_adj<=0.1]
DoHeatmap(dge_temp, features = c(heatmap_gene),raster = F) + theme(axis.text.y = element_text(size = 6))#+ scale_fill_gradientn(colors = c("white", "red"))
ggsave(file="./Heatmap_Club_Pair_purple.pdf",width = 20,height = 30,units = "cm")
write.table(top50_Club_Pair,"Club_Pair_top50_DEG.txt",col.names = T,row.names = T)

M1_macro_features=c("CCR7","IL2RA","IL15RA","IL7R","CXCL11","CCL19","CXCL10","CXCL9",
                    "TNF","CCL5","IL12B","IL15","IL6","CCL20","BCL2A1","RAS","BIRC3",
                    "GADD45G","SLC7A5","SLC2A6","SLC31A2","PLA1A","OASL","CHI3L2","HSD11B1",
                    "AK3","SPHK1","PFKFB3","PSME2","PFKP","PSMB9","PSMA2","OAS2")
M2_macro_features=c("TGFBR2","HRH1","TLR5","MSR1","CXCR4","P2RY14","MS4A6A",
                    "CD36","MS4A4A","MRC1","IGF1","CCL23","CCL18","SLC4A7","SLC38A6",
                    "CTSC","HEXB","LIPA","ADK","HNMT","TPST2","CERK","HS3ST2",
                    "LTA4H","CA2","ALOX15","HS3ST1","TGFBI","SEPP1","CHN2",
                    "FN1","FGL2","GAS7","EGR2","MAF")
TAM_features=c("CCL5","IL10","CSF1R","TGFB1","CCL2","CCL17","CCL22","IL6",
               "CXCL10","CXCL12","CD274","PDCD1LG","NOS2","MRC1","VTCN1","TGFB2")

DotPlot(dge, features =c("CD68","CD163","MRC1","IL1B","NLRP3","C1QC","SPP1","PLTP","LYVE1","IL10","HLA-DRA")) + RotatedAxis()
DotPlot(dge, features =c("SLAMF7","DUSP4","SDS","MT1G","MT2A","IRF1","CXCR2","TNFRSF1A","CXCL10","IFNG","THBS1","LYVE1","FCN1","VCAN","SEPP1")) + RotatedAxis()
DotPlot(dge, features = c("IL1A","CXCL3","PTGS2","ARG1","CCL22","FLT1")) +  RotatedAxis()

ggsave(file="dot_plot_bycluster.eps",width = 70,height = 50,units = "cm")

dge <- ScaleData(object=dge,features=rownames(dge))
DoHeatmap(dge, features = top5$gene,raster = F) + theme(axis.text.y = element_text(size = 6))#+ scale_fill_gradientn(colors = c("white", "red"))

DoHeatmap(dge, features = c("CXCL10","CCR7","IL6","TGFBR2","CXCR4","CD36","C1QC","SPP1","FCN1","IL10","CCL2","NLRP3","PLTP","IL1B"),raster = F) + theme(axis.text.y = element_text(size = 6))#+ scale_fill_gradientn(colors = c("white", "red"))

cell_IL1B<-colnames(SubsetData(dge,subset.name = "IL1B",low.threshold = 0))
cell_NLRP3<-colnames(SubsetData(dge,subset.name = "NLRP3",low.threshold = 0))
cell_PLTP<-colnames(SubsetData(dge,subset.name = "PLTP",low.threshold = 0))
cell_C1QC<-colnames(SubsetData(dge,subset.name = "C1QC",low.threshold = 0))
cell_SPP1<-colnames(SubsetData(dge,subset.name = "SPP1",low.threshold = 0))
cell_FCN1<-colnames(SubsetData(dge,subset.name = "FCN1",low.threshold = 0))
cell_CSF1R<-colnames(SubsetData(dge,subset.name = "CSF1R",low.threshold = 0))

DimPlot(dge_myeloid,cells.highlight = cell_NLRP3)
ggsave(file="NLRP3_RPM.eps",width = 20,height = 20,units = "cm")
DimPlot(dge_myeloid,cells.highlight = intersect(cell_IL1B,cell_PLTP))
ggsave(file="PLTP_RPM.eps",width = 20,height = 20,units = "cm")
DimPlot(dge_myeloid,cells.highlight = cell_C1QC)
ggsave(file="TAM_threemarkers.eps",width = 20,height = 20,units = "cm")

DimPlot(dge_myeloid, cells.highlight= list(cell_C1QC, cell_SPP1,cell_FCN1), cols.highlight = c("red", "blue","green"), cols= "grey")
ggsave(file="TAM_highlight_2019.eps",width = 20,height = 20,units = "cm")

FeaturePlot(dge_myeloid,features=c("IL1B","PLTP"),blend = T)
FeaturePlot(dge_myeloid,features=c("IL1B","NLRP3"),blend = T)
ggsave(file="Monocyte_RTM.eps",width = 20,height = 20,units = "cm")
FeaturePlot(dge_myeloid,features=c("MKI67","TOP2A"))
ggsave(file="Proliferating_Myeloid.eps",width = 40,height = 20,units = "cm")
FeaturePlot(dge_myeloid,features=c("SPP1","C1QC"),blend = T)

FeaturePlot(dge_myeloid,features=c("FCN1","C1QC"),blend = T)



dge_Tcell$ERGtype<-"ERGneg"
dge_Tcell$ERGtype[dge_Tcell$orig.ident %in% c("AUG_PB1","MAY_PB1","PR5249_T","PR5249_N","PR5251_T","PR5251_N","PR5199","PR5269")]="ERGpos"

dge_temp<-SubsetData(dge_Tcell,ident.use = c("CD8","CD8_cluster3","CD8_cluster2"))
Featurelist<-c("CD247","ENTPD1","GZMA","GZMB","GZMK")
Featurelist<-c("TIGIT","HAVCR2","CTLA4","LAG3","ITGA4")
Featurelist<-c("TGFB1","TGFB2","CXCR4","IFNG","IFNGR1","IFNGR2","TGFBR1","TGFBR2")
Featurelist<-c("PDCD1","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1","HLA-DRB5","STAT1","MT2A","HLA-A","HLA-B","HLA-C",
               "HLA-E","HLA-F","HLA-G","ICAM1","CD44","B2M","GBP2","GBP4","GBP5")

###cytotoxic
Featurelist<-c("CST7","GZMA","GZMB","IFNG","NKG7","PRF1","TNFSF10")
###Exhausted
Featurelist<-c("BTLA","CTLA4","HAVCR2","LAG3","PDCD1","TIGIT")
###Regulatory
Featurelist<-c("IL2RA","IL4R","IL7","TGFB1","TGFB3","TGFBI","TGFBR1")
###Naive
Featurelist<-c("CCR7","LEF1","SELL","TCF7")
###Costimulatory
Featurelist<-c("ICOS","CD226","SLAMF1","TNFRSF14","TNFRSF25","TNFRSF9")


Pos_cells <- colnames(dge_temp)[dge_temp$ERGtype == "ERGpos"]
Pos_data <- FetchData(dge_temp,
                  vars = Featurelist,
                  cells = Pos_cells ,
                  slot = "data")
Neg_cells <- colnames(dge_temp)[dge_temp$ERGtype == "ERGneg"]
Neg_data <- FetchData(dge_temp,
                  vars = Featurelist,
                  cells = Neg_cells ,
                  slot = "data")
Pos_data<-as.data.frame(Pos_data)
Neg_data<-as.data.frame(Neg_data)
colnames(Pos_data)<-paste0(colnames(Pos_data),"_Pos")
colnames(Neg_data)<-paste0(colnames(Neg_data),"_Neg")
long_data_Pos <- melt(Pos_data)
long_data_Neg <- melt(Neg_data)
long_data<-rbind(long_data_Pos,long_data_Neg)
long_data<-as.data.frame(long_data)
long_data <- long_data[order(as.vector(long_data$variable)),]

ggplot(long_data,
       aes(x = reorder(variable), y = value)) +
  geom_violin() +
  geom_jitter(size = 0.1)
ggsave(file="CD4_Cytotoxicmarkers_violinplot.eps",width = 45,height =15,units = "cm")


Featurelist<-c("CST7","GZMA","GZMB","IFNG","NKG7","PRF1","TNFSF10","BTLA","CTLA4","HAVCR2","LAG3","TIGIT","PDCD1")
###Exhausted

dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$ERGtype))
StackedVlnPlot(dge_temp,features=Featurelist)+geom_jitter(size=0.1)+geom_boxplot()
ggsave(file="CD8_Exhaustedmarkers_stack.eps",width = 10,height =20,units = "cm")

dge_CD4$pair="Others"
dge_CD4$pair[dge_CD4$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")]="pairN"
dge_CD4$pair[dge_CD4$orig.ident %in% c("PR5249_T","PR5251_T","PR5254_T","PR5261_T")]="pairT"
dge_temp<-SetIdent(dge_CD4,value = as.vector(dge_CD4$pair))
dge_temp<-SubsetData(dge_temp,ident.remove = "Others")
VlnPlot(dge_temp,features = Featurelist,ncol = 3)
ggsave(file="~/Desktop/Seqwell_combined/Tcell_analysis/CD8_pair_cytotoxicmarkers_violinplot.eps",width = 30,height =20,units = "cm")

dge_temp<-SubsetData(dge_ERGstatus,ident.use = c("ERGpos_CD4","ERGneg_CD4"))
FC_df<-NULL
for (i in 1:length(Featurelist)) {
  name_feature=Featurelist[i]
  #FeaturePlot(LE_Henry,features=Featurelist[i])+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, 5))
  #ggsave(file=paste0(name_feature,"_Henry_LE_Cluster.eps"),width = 20,height = 20,units = "cm")
  V1<-as.vector(FetchData(object = dge_temp,vars = name_feature))
  BOX_df<-NULL
  #BOX_df$id<-dge_temp$target
  #BOX_df$id<-as.vector(dge_temp$ERGtype)
  BOX_df$id<-as.vector(dge_temp$pair)
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  # Pos_mean=mean(BOX_df$value[BOX_df$id=="ERGpos"])
  # Neg_mean=mean(BOX_df$value[BOX_df$id=="ERGneg"])
  Pos_mean=mean(BOX_df$value[BOX_df$id=="pairN"])
  Neg_mean=mean(BOX_df$value[BOX_df$id=="pairT"])
  FC=log(Neg_mean/Pos_mean)/log(2)
  FC_temp<-c(Pos_mean,Neg_mean,FC)
  FC_df<-rbind(FC_df,FC_temp)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    #geom_signif(comparisons = list(c("ERGneg","ERGpos")),map_signif_level=TRUE,y_position = max(V1)+0.25)
    geom_signif(comparisons = list(c("pairN","pairT")),map_signif_level=TRUE,y_position = max(V1)+0.25)
  #geom_signif(comparisons = list(c("Target","Others")),map_signif_level=TRUE,y_position = max(V1)+0.25) +
  ggtitle(name_feature)+ RotatedAxis()
  ggsave(file=paste0(name_feature,"_in_CD4.pdf"),width = 20,height = 15,units = "cm")
} 

FC_df<-as.data.frame(FC_df)
#colnames(FC_df)=c("Pos_mean","Neg_mean","log2FC")
colnames(FC_df)=c("PairN_mean","PairT_mean","log2FC")
rownames(FC_df)=as.vector(Featurelist)
write.table(FC_df,"CD4_DEG_Exhaust_cyto.txt",sep="\t",quote = F,col.names = T,row.names = T)

Featurelist<-c("CCR7","IL2RA","TNFRSF4","FOXP3","CD69","GZMA","GZMB","GZMK","PRF1","IL17A","CXCL13","HSPA1A","TIGIT","STMN1","MKI67","IFNG","TOX")
DF=NULL
for (i in 12:length(Featurelist)) {
  dge_temp<-SubsetData(dge_CD8,subset.name = Featurelist[i],low.threshold = 0)
  tab_temp<-table(dge_temp$ERGtype)
  mtx_temp<-matrix(c(as.numeric(tab_temp[1]),as.numeric(tab_temp[2]),146-as.numeric(tab_temp[1]),323-as.numeric(tab_temp[2])),nrow = 2)
  fisher_temp<-fisher.test(mtx_temp)
  DF_temp<-c(tab_temp,fisher_temp$p.value)
  DF<-rbind(DF,DF_temp)
}
DF<-as.data.frame(DF)
rownames(DF)<-Featurelist
write.table(DF,"CD4_subtype_markers.txt",row.names = T,col.names = T,sep = "\t")

dge_temp<-dge_CD8
cell_CCR7<-SubsetData(dge_temp,subset.name = "CCR7",low.threshold = 0)
cell_TOX<-SubsetData(dge_temp,subset.name = "TOX",low.threshold = 0)
cell_IL17A<-SubsetData(dge_temp,subset.name = "IL17A",low.threshold = 0)
cell_CD69<-SubsetData(dge_temp,subset.name = "CD69",low.threshold = 0)
cell_GZMB<-SubsetData(dge_temp,subset.name = "GZMB",low.threshold = 0)
cell_GZMK<-SubsetData(dge_temp,subset.name = "GZMK",low.threshold = 0)
cell_IL2RA<-SubsetData(dge_temp,subset.name = "IL2RA",low.threshold = 0)
cell_MKI67<-SubsetData(dge_temp,subset.name = "MKI67",low.threshold = 0)
cell_CXCL13<-SubsetData(dge_temp,subset.name = "CXCL13",low.threshold = 0)
cell_HSP<-SubsetData(dge_temp,subset.name = "HSPA1A",low.threshold = 0)

dge_temp$ID="CD4_IL2RAlow"
dge_temp$ID[colnames(dge_temp) %in% colnames(cell_HSP)]="CD4_HSP"
dge_temp$ID[colnames(dge_temp) %in% colnames(cell_CD69)]="CD4_Activated"
dge_temp$ID[colnames(dge_temp) %in% colnames(cell_CCR7)]="CD4_CM"
dge_temp$ID[colnames(dge_temp) %in% colnames(cell_TOX)]="CD4_Exhausted"
dge_temp$ID[colnames(dge_temp) %in% colnames(cell_IL17A)]="CD4_Th17"
dge_temp$ID[colnames(dge_temp) %in% colnames(cell_IL2RA)]="CD4_IL2RAhi"
dge_temp$ID[colnames(dge_temp) %in% c(colnames(cell_GZMB),colnames(cell_GZMK))]="CD4_Cytotoxic"
dge_temp$ID[colnames(dge_temp) %in% colnames(cell_MKI67)]="CD4_Proliferation"
dge_temp$ID[colnames(dge_temp) %in% colnames(cell_CXCL13)]="CD4_CXCL13"
DimPlot(dge_temp,group.by = "ID")
cytotoxic_pathway<-c("CST7","GZMA","GZMB","GZMK","IFNG","NKG7","PRF1","TNFSF10")
Exhausted_pathway<-c("BTLA","CTLA4","HAVCR2","LAG3","PDCD1","TIGIT","TOX")

dge_temp <- AddModuleScore(object = dge_temp, features = list(pd1_feature_filter), name = "PD1_pathway",assay = "RNA")
dge_temp <- AddModuleScore(object = dge_temp, features = list(interferon_gamma_feature_filter), name = "Interferon_gamma",assay = "RNA")
dge_temp <- AddModuleScore(object = dge_temp, features = list(cytotoxic_pathway), name = "cytotoxic",assay = "RNA")
dge_temp <- AddModuleScore(object = dge_temp, features = list(Exhausted_pathway), name = "exhausted",assay = "RNA")

VlnPlot(dge_temp,features=c("PD1_pathway1"),group.by = "ID")
VlnPlot(dge_temp,features=c("Interferon_gamma1"),group.by = "ID")
VlnPlot(dge_temp,features=c("cytotoxic1"),group.by = "ID")
VlnPlot(dge_temp,features=c("exhausted1"),group.by = "ID")

VlnPlot(dge_temp,features=c("PD1_pathway1"),group.by = "ERGtype")
VlnPlot(dge_temp,features=c("Interferon_gamma1"),group.by = "ERGtype")
VlnPlot(dge_temp,features=c("exhausted1"),group.by = "ERGtype")
VlnPlot(dge_temp,features=c("cytotoxic1"),group.by = "ERGtype")

# mean(dge_temp$exhausted1[dge_temp$ERGtype=="Pos"])
# mean(dge_temp$exhausted1[dge_temp$ERGtype=="Neg"])
# 
# mean(dge_temp$cytotoxic1[dge_temp$ERGtype=="Pos"])
# mean(dge_temp$cytotoxic1[dge_temp$ERGtype=="Neg"])
# 
# mean(dge_temp$PD1_pathway1[dge_temp$ERGtype=="Pos"])
# mean(dge_temp$PD1_pathway1[dge_temp$ERGtype=="Neg"])
# 
# mean(dge_temp$Interferon_gamma_pathway1[dge_temp$ERGtype=="Pos"])
# mean(dge_temp$Interferon_gamma_pathway1[dge_temp$ERGtype=="Neg"])

mean(dge_temp$exhausted1[dge_temp$ERGtype=="ERGpos"])
mean(dge_temp$exhausted1[dge_temp$ERGtype=="ERGneg"])

mean(dge_temp$cytotoxic1[dge_temp$ERGtype=="ERGpos"])
mean(dge_temp$cytotoxic1[dge_temp$ERGtype=="ERGneg"])

mean(dge_temp$PD1_pathway1[dge_temp$ERGtype=="ERGpos"])
mean(dge_temp$PD1_pathway1[dge_temp$ERGtype=="ERGneg"])

mean(dge_temp$Interferon_gamma1[dge_temp$ERGtype=="ERGpos"])
mean(dge_temp$Interferon_gamma1[dge_temp$ERGtype=="ERGneg"])

interferon_gamma_signature<-as.vector(dge_temp$Interferon_gamma_pathway1)
expression_GZMB<-FetchData(dge_temp,vars = "GZMB")
cor.test(as.vector(expression_GZMB$GZMB),interferon_gamma_signature)
expression_GZMA<-FetchData(dge_temp,vars = "GZMA")
cor.test(as.vector(expression_GZMA$GZMA),interferon_gamma_signature)
expression_GZMK<-FetchData(dge_temp,vars = "GZMK")
cor.test(as.vector(expression_GZMK$GZMK),interferon_gamma_signature)

dge<-dge_CD8

dge$ERGtype_new<-"ERGneg"
dge$ERGtype_new[dge$orig.ident %in% c("PR5199","PR5269")]="ERGpos"
dge$ERGtype_new[dge$orig.ident %in% c("AUG_PB1","MAY_PB1","PR5249_T","PR5249_N","PR5251_T","PR5251_N")]="Mix"
DimPlot(dge,group.by = "ERGtype_new")
dge_CD8<-dge


dge<-dge_CD8
dge<-SetIdent(dge,value = as.vector(dge$ERGtype))

dge_temp<-SubsetData(dge,ident.remove = c("Mix"))
source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
#dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$ERGtype_new))
Data<-dge
name1<-"Ident"
name2<-"HBO_epithelial"
Prepare.for.GSEA(Data,name1,name2)

dge_temp<-SubsetData(dge,ident.remove = c("ERGpos"))
source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
dge_temp<-SetIdent(dge_temp,value=as.vector(dge_temp$ERGtype_new))
Data<-dge_temp
name1<-"Ident"
name2<-"ERGneg_Mix_CD8"
Prepare.for.GSEA(Data,name1,name2)

dge_temp<-SubsetData(dge,ident.remove = c("ERGneg"))
source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
dge_temp<-SetIdent(dge_temp,value=as.vector(dge_temp$ERGtype_new))
Data<-dge_temp
name1<-"Ident"
name2<-"Mix_ERGpos_CD8"
Prepare.for.GSEA(Data,name1,name2)





