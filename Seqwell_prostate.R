#HS 2019

library(Seurat)
library(cowplot)
library(dplyr)
library(icesTAF)
library(GEOquery)
library(ggplot2)
library(data.table)
theme_set(theme_cowplot())
#samples <- as.vector(read.csv("../Filename.txt",header=FALSE,sep="\n")$V1)
#l=length(samples)

setwd("/Users/hsong/Desktop/Seqwell_0913/human_output/")
mkdir("prostate")
#Kidney samples
samples <- as.vector(read.csv("Filename_prostate.txt",header=FALSE,sep="\n")$V1)
for (i in 2:length(samples)) {
zz=gunzip(paste0("dge/",samples[i],"_dge.txt.gz")) 
}
dge<-list()

dge_summary=NULL
for (i in 1:length(samples)) {
dg_temp <- read.table(file = paste0("dge/",samples[i],"_dge.txt"), header = TRUE,row.names = 1)
dge[[i]] <- CreateSeuratObject(counts = dg_temp, project = samples[i], min.cells = 3,min.features = 0)
saveRDS(dge[[i]], file = paste0("prostate/",samples[i],"_Seurat.rds"))
dge_summary$sample[i]=samples[i]
dge_summary$median_gene[i]=median(dge[[i]]$nFeature_RNA)
dge_summary$mean_gene[i]=mean(dge[[i]]$nFeature_RNA)
dge_summary$median_UMI[i]=median(dge[[i]]$nCount_RNA)
dge_summary$mean_UMI[i]=mean(dge[[i]]$nCount_RNA)
dge[[i]][["percent.mt"]] <- PercentageFeatureSet(object = dge[[i]], pattern = "^MT-")
dge_summary$mean_mt[i]=mean(na.omit(dge[[i]]$percent.mt))
}

write.table(dge_summary, file = "prostate_summary.txt", sep = "\t",row.names = FALSE,col.names = TRUE)
dg_temp <- read.table(file = "./2wk_count.txt", header = TRUE,row.names = 1)
setwd("./prostate/")
dge<-list()
for (i in 1:length(samples)) {

  dge[[i]]<-readRDS(paste0(samples[i],"_Seurat.rds"))
  dge[[i]][["percent.mt"]] <- PercentageFeatureSet(object = dge[[i]], pattern = "^MT-")
  VlnPlot(object = dge[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file=paste0(samples[i],"_VlnPlot.pdf"),width = 20,height = 20,units = "cm")

  plot1 <- FeatureScatter(object = dge[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object = dge[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  ggsave(file=paste0(samples[i],"Scatter.eps"),width = 30,height = 30,units = "cm")
}


#merge 2 objects
dge_combine<-dge[[1]]
for (i in 1:(length(samples)-1)) {
#dge_combine <- merge(dge_combine,y=c(dge[[2]],dge[[3]],dge[[4]],dge[[5]],),add.cell.ids=c("11603T","12050T","11734","12050N","12041"))
  dge_combine <- merge(dge_combine,y=dge[[i+1]],add.cell.ids=c(samples[i],samples[i+1]))
}
saveRDS(dge_combine, file = "Combined_Seurat.rds")
dge_combine<-readRDS("Combined_Seurat.rds")
dge_combine[["percent.mt"]] <- PercentageFeatureSet(object = dge_combine, pattern = "^MT-")

head(colnames(dge_combine))
table(dge_combine$orig.ident)
# Set up control object
median(dge_combine$nFeature_RNA)
mean(dge_combine$nFeature_RNA)
median(dge_combine$nCount_RNA)
mean(dge_combine$nCount_RNA)
mean(na.omit(dge_combine$percent.mt))

dge <- dge_combine
dge_backup <- dge


VlnPlot(object = dge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")
ggsave(file="VlnPlot_afterfilter.eps",width = 20,height = 20,units = "cm")

print(paste("Median nGene:", median(dge$nFeature_RNA), sep = " "))
print(paste("Median nUMI:", median(dge$nCount_RNA), sep = " "))
dge_backup <- dge

#dge<-SCTransform(dge, verbose = FALSE)

dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
head(VariableFeatures(object=dge))
head(HVFInfo(object = dge))
dge_normalized <-dge

dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = dge), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = dge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge = 0,ynudge = 0)
CombinePlots(plots = list(plot1, plot2))
ggsave(file="VariableFeature.eps",width = 20,height = 20,units = "cm")

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

VizDimLoadings(object = dge, dims = 1:20, reduction = "pca")
ggsave(file="PCA1_20.eps",width = 100,height = 100,units = "cm")
DimPlot(object = dge, reduction = "pca")
PCAPlot(object=dge,dim.1=1,dim.2=2)
ggsave(file="PCA.eps",width = 20,height = 20,units = "cm")

png("PCA_heatmap.png")
DimHeatmap(object = dge, dims = 1:20, cells = length(dge$nCount_RNA), balanced = TRUE)
dev.off()
png("Elbowplot_dge.png")
ElbowPlot(dge)
dev.off()

dge <- JackStraw(dge, num.replicate = 100,dims = 50)
dge <- ScoreJackStraw(dge, dims = 1:50)
JackStrawPlot(dge, dims = 1:50)
ggsave(file="Strawplot.eps",width = 20,height = 20,units = "cm")
ElbowPlot(dge)
ggsave(file="Elbow.eps",width = 10,height = 10,units = "cm")

n_pc=43
resolut=0.8
dge_precluster <- dge
dge<-tenx
dge <- FindNeighbors(dge, dims = 1:n_pc,k.param = 50)
dge <- FindClusters(dge, resolution = resolut)


dge <- FindClusters(dge, resolution =c(0.25,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,2))

dge <- FindClusters(dge, resolution =1)

pdf("Cluster_stability.pdf",width = 20,height = 15)
clustree(dge, prefix = "RNA_snn_res.")
dev.off()
#head(Idents(dge), 10)

dge=RunTSNE(dge,dims.use = 1:n_pc,perplexity=40,seed.use = 10)
TSNEPlot(dge,pt.size = 1,label=T)
ggsave(file="TSNE.eps",width = 20,height = 20,units = "cm")

dge <- RunUMAP(dge, dims = 1:n_pc)
DimPlot(dge, reduction = "umap",label=T)
ggsave(file="Umap.eps",width = 20,height = 20,units = "cm")

DimPlot(dge, reduction = "umap",group.by="orig.ident",label=F)
ggsave(file="Umap_group.eps",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "tsne",group.by="orig.ident",label=F)
ggsave(file="TSNE_group.eps",width = 20,height = 20,units = "cm")
# Print number of cells per cluster



saveRDS(dge, file = "Merged.rds")
dge<-SetIdent(dge,value = as.vector(dge$orig.ident))

pbmc.markers <- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)


cluster1.markers <- FindMarkers(dge, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers,10)
cluster2.markers <- FindMarkers(dge, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster2.markers,10)
cluster3.markers <- FindMarkers(dge, ident.1 = 2, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster3.markers,10)
cluster4.markers <- FindMarkers(dge, ident.1 = 3, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster4.markers,10)
cluster5.markers <- FindMarkers(dge, ident.1 = 4, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster5.markers,10)

cluster6.markers <- FindMarkers(dge, ident.1 = 5, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster6.markers,10)
cluster7.markers <- FindMarkers(dge, ident.1 = 6, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster7.markers,10)
cluster8.markers <- FindMarkers(dge, ident.1 = 7, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster8.markers,10)
cluster9.markers <- FindMarkers(dge, ident.1 = 8, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster9.markers,10)
cluster10.markers <- FindMarkers(dge, ident.1 = 9, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster10.markers,10)

cluster11.markers <- FindMarkers(dge, ident.1 = 10, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster11.markers,10)
cluster12.markers <- FindMarkers(dge, ident.1 = 11, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster12.markers,10)
cluster13.markers <- FindMarkers(dge, ident.1 = 12, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster13.markers,10)
cluster14.markers <- FindMarkers(dge, ident.1 = 13, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster14.markers,10)
cluster15.markers <- FindMarkers(dge, ident.1 = 14, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster15.markers,10)

cluster16.markers <- FindMarkers(dge, ident.1 = 15, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster16.markers,10)
cluster17.markers <- FindMarkers(dge, ident.1 = 16, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster17.markers,10)
cluster18.markers <- FindMarkers(dge, ident.1 = 17, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster18.markers,10)



VlnPlot(dge, features = c(rownames(head(cluster1.markers,5)),rownames(head(cluster2.markers,5)),rownames(head(cluster3.markers,5)),rownames(head(cluster4.markers,5)),rownames(head(cluster5.markers,5))),slot = "counts", log = TRUE)
ggsave(file="VlnPlot_features_count.eps",width = 40,height = 40,units = "cm")
VlnPlot(object = dge, features = c("KRT14","DST","KRT15","KRT5","RGCC","MSMB","KLK3","ACPP","PLA2G2A","KLK2","SCGB3A1","LCN2","PIGR","PTEN","ERG","TMPRSS2","CCND1","WFDC2","FCGBP","KRT13","APOBEC3A","CSTB","LYPD3","LY6D","CHGA","GRP","CALCA","SCG2","TPH1","APOD","FBLN1","PTGDS","CFD","DCN","TPM2","ACTA2","RGF5","MT1A","MYH11","IFI27","ACKR1","SELE","HMOX1","CLDN5","RGS1","C1QA","C1QB","TYROBP","C1QC","CD45","CD200","CD326","PDPN","CD26","PSCA"))
ggsave(file="VlnPlot_features_marker.png",width = 100,height = 100,units = "cm")
VlnPlot(object = dge, features = c("P63"))
ggsave(file="AR.eps",width = 40,height = 40,units = "cm")

FeaturePlot(dge, features = c(rownames(head(cluster1.markers,5)),rownames(head(cluster2.markers,5)),rownames(head(cluster3.markers,5)),rownames(head(cluster4.markers,5)),rownames(head(cluster5.markers,5))))
ggsave(file="FeaturePlot.eps",width = 40,height = 40,units = "cm")
FeaturePlot(dge, features = c("KRT14","DST","KRT15","KRT5","RGCC","MSMB","KLK3","ACPP","PLA2G2A","KLK2","SCGB3A1","LCN2","PIGR","PTEN","ERG","TMPRSS2","CCND1","WFDC2","FCGBP","KRT13","APOBEC3A","CSTB","LYPD3","LY6D","CHGA","GRP","CALCA","SCG2","TPH1","APOD","FBLN1","PTGDS","CFD","DCN","TPM2","ACTA2","RGF5","MT1A","MYH11","IFI27","ACKR1","SELE","HMOX1","CLDN5","RGS1","C1QA","C1QB","TYROBP","C1QC","CD45","CD200","CD326","PDPN","CD26","PSCA"))
ggsave(file="FeaturePlot_marker.png",width = 100,height = 100,units = "cm")

VlnPlot(object = dge, features = c("EPCAM"))
ggsave(file="CD326.eps",width = 40,height = 40,units = "cm")


av.exp <- AverageExpression(dge,)$RNA
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
#
cor.df <- tidyr::gather(data = cor.exp, y, correlation, as.vector(unique(dge$seurat_clusters)))

#cor.df <- tidyr::gather(data = cor.exp, y, correlation, c("T-cell","Macrophage","Mast","NK-cell","Plasma","Exhausted_Tcell","Leu"))
#cor.df <- tidyr::gather(data = cor.exp, y, correlation, c("T-cell","Macrophage","LE","Monocyte","Mast","NK","MZB","Cell-cycle"))
ggplot(cor.df, aes(x, y),lab_size = 5) +geom_tile(aes(fill = correlation))+scale_fill_gradientn(colours = c("blue","cyan","lightcyan","red"))
ggsave(file="Cluster_correlation.eps",width = 20,height = 20,units = "cm")


library(ggplot2)
library(reshape2)
library(RColorBrewer)
#Get table listing number of cells from each mouse per cluster
#Note: subs is my Seurat object
cells_per_mouse <- table(dge@active.ident, dge@meta.data$orig.ident)

#Get percent contribution per mouse for each cluster 
#Divide cells per mouse by total cells from that cluster
percent_cells_mouse <- apply(cells_per_mouse,1,function(x){x/sum(x)})

#Reformat the table for it to work with ggplot2
percent_cells_mouse <- as.data.frame(t(percent_cells_mouse))
percent_cells_mouse$cluster <- rownames(percent_cells_mouse)
percent_cells_mouse <- melt(percent_cells_mouse, id.var="cluster")

percent_cells_mouse$variable =factor(percent_cells_mouse$variable, levels = c("concat_HB008", "concat_HB012", "concat_HB016"))

#Store graph

bargraph <- ggplot(percent_cells_mouse, aes(x = cluster, y = value, fill = variable)) + geom_bar(stat = "identity") + xlab("Cluster") + ylab("Percent of cells in cluster") + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +  theme_classic()

#Plot graph
bargraph
ggsave(file="Barchart.eps",width = 15,height = 15,units = "cm")


top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(top10,"top10_markers_bypatient.txt",col.names = T,row.names = T,quote = F,sep="\t")
dge <- ScaleData(object=dge,features=rownames(dge))
DoHeatmap(dge, features = top10$gene,combine = T,raster = F) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_10_bypatient.pdf",width = 15,height = 15,units = "cm")

FeaturePlot(dge,features=c("EPCAM","PTPRC"))
ggsave(file="CD326_CD45.eps",width = 20,height = 20,units = "cm")

CombinePlots(plots = list(plot1, plot2,plot3,plot4, plot5,plot6,plot7, plot8,plot9))
#DimPlot(object=dge,cells.highlight = WhichCells(object=dge,expression=orig.ident =="PA_PB1A_Pool_1_3_S50_L002"),reduction = "umap",label = TRUE, pt.size = 0.5)
ggsave(file="Array_Cluster.eps",width = 100,height = 100,units = "cm")

VlnPlot(object = dge, features = c("MSMB","KLK3","ACPP","PLA2G2A","KLK2","DBI","EFNA1","IDH2","FDPS","IL6ST"))
ggsave(file="LE_markers.eps",width = 40,height = 40,units = "cm")

VlnPlot(object = dge, features = c("DST","ZFP36L2","CAV1","CAV2","CLU","DEFB1","HBEGF","EDN1","F3","FHL2","KRT5","KRT14","KRT15","KRT17","S100A2","STK17A","PHLDA1","RGCC","HSPA14"))
ggsave(file="BE_markers.eps",width = 40,height = 40,units = "cm")

VlnPlot(object = dge, features = c("BIK","CEACAM5","CP","CRABP2","CYP1B1","FOXO3","FTH1","GLUL","LCN2","MMP7"))
ggsave(file="Club_markers.eps",width = 40,height = 40,units = "cm")

VlnPlot(object = dge, features = c("AOC1","ADARB2","ADCY2","ANK1","AOX1","ASCL1","HCN2","CBLB","CLGN","CHN1","CGA"))
ggsave(file="Hillock_markers.eps",width = 40,height = 40,units = "cm")

#VlnPlot(object = dge, features = c("SCGB3A1","LCN2","PIGR","WFDC2","FCGBP"))
#ggsave(file="OE1_markers.eps",width = 40,height = 40,units = "cm")

#VlnPlot(object = dge, features = c("KRT13","APOBEC3A","CSTB","LYPD3","LY6D"))
#ggsave(file="OE2_markers.eps",width = 40,height = 40,units = "cm")

VlnPlot(object = dge, features = c("AQP3","CDKN2B","CSTB","KRT13","IL1RN","TRIM3","LYPD3","TMPRSS4","CLCA4","CCL20"))
ggsave(file="NE_markers.eps",width = 40,height = 40,units = "cm")

VlnPlot(object = dge, features = c("APOD","FBLN1","PTGDS","CFD","DCN","TPM2","RSPO3","FGF2","FBLN2","GSN","GSTM3","CFH","CXCL1"))
ggsave(file="Fib_markers.eps",width = 40,height = 40,units = "cm")

VlnPlot(object = dge, features = c("TPM2","ACTA2","RGS5","MT1A","MYH11","APP","RHOC","ASPH","CRIP2","CD59","CXCL2","HLA-DRB5","CD200","PNP","MSN"))
ggsave(file="SM_markers.eps",width = 40,height = 40,units = "cm")

VlnPlot(object = dge, features = c("IFI27","ACKR1","SELE","VWF"),ncol=2)
ggsave(file="Endo_markers.eps",width = 20,height = 20,units = "cm")

VlnPlot(object = dge, features = c("RGS1","C1QA","C1QB","TYROBP","C1QC"))
ggsave(file="Leu_markers.eps",width = 40,height = 40,units = "cm")

VlnPlot(object = dge, features = c("CD74"))
ggsave(file="CD74.eps",width = 40,height = 40,units = "cm")

#CD326: EPCAM
#CD26: DPP4, DPP8, DPP9, FAP
#CD45: PTPRC, IL7R, FYN, LYN



for(j in 1:ncol(Pliver@meta.data)){
  if(is.factor(Pliver@meta.data[,j]) == T){
    Pliver@meta.data[,j] = as.character(Pliver@meta.data[,j]) # Force the variable type to be character
    Pliver@meta.data[,j][is.na(Pliver@meta.data[,j])] <- "N.A"
  }
  if(is.character(Pliver@meta.data[,j]) == T){
    Pliver@meta.data[,j][is.na(Pliver@meta.data[,j])] <- "N.A"
  }
}
Pliver.loom <- as.loom(Pliver, filename = "~/Desktop/hepatoblastoma/scRNA/Pliver.loom", verbose = T)
Pliver.loom



