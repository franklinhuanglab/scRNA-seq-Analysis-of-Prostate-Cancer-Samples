setwd("/Users/hsong/Desktop/Seqwell_combined/Tcell_analysis/")
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
library(ggpubr)
library(ggsignif)
theme_set(theme_cowplot())

dge<-readRDS("~/Desktop/Seqwell_combined/dge_")
dge<-SubsetData(dge_Tcell,ident.remove = "KLK3+")
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 2000)
top20 <- head(x = VariableFeatures(object = dge), 20)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = dge)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE,xnudge = 0,ynudge = 0)
CombinePlots(plots = list(plot1, plot2))
ggsave(file="VariableFeature.pdf",width = 20,height = 20,units = "cm")
all.genes <- rownames(x = dge)
dge <- ScaleData(object=dge,features=VariableFeatures(object = dge))
dge <- RunPCA(object = dge, features = VariableFeatures(object = dge),npcs = 100)

DimPlot(object = dge, reduction = "pca",group.by = "orig.ident")
PCAPlot(object=dge,dim.1=1,dim.2=2)
ggsave(file="PCA_id.pdf",width = 20,height = 20,units = "cm")

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

dge <- JackStraw(dge, num.replicate = 100,dims = 100)
dge <- ScoreJackStraw(dge, dims = 1:100)
JackStrawPlot(dge, dims = 1:100)
ggsave(file="Strawplot.pdf",width = 30,height = 30,units = "cm")
ElbowPlot(dge)
ggsave(file="Elbow.eps",width = 10,height = 10,units = "cm")



n_pc=28
resolut=0.3
dge_precluster <- dge
dge<-dge_precluster
dge <- FindNeighbors(dge, dims = 1:n_pc,k.param = 50)
dge <- FindClusters(dge, resolution = resolut)

dge=RunTSNE(dge,dims.use = 1:n_pc,perplexity=40,seed.use = 10)
TSNEPlot(dge,pt.size = 1,label=T)
ggsave(file="TSNE.pdf",width = 20,height = 20,units = "cm")
dge <- RunUMAP(dge, dims = 1:n_pc)
DimPlot(dge, reduction = "umap",label=T)
ggsave(file="Umap.pdf",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "tsne",label=T,group.by = "orig.ident")
ggsave(file="TSNE_group.eps",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "umap",label=T,group.by = "orig.ident")
ggsave(file="Umap_group.eps",width = 20,height = 20,units = "cm")

Tab<-table(dge_Tcell$orig.ident,Idents(object = dge_Tcell))
write.table(Tab,"../table_IDclustering.txt",sep="\t",row.names = T,col.names = T)
dge<-SetIdent(dge,value = as.vector(dge$ID))
Tcell.markers <- FindAllMarkers(dge_CD4, only.pos = TRUE, logfc.threshold = 0.001,min.pct = 0.001)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
#RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
Tcell.markers<-Tcell.markers[!(Tcell.markers$gene %in% remove_gene),]
write.table(Tcell.markers, file = "Tcell_dge_markers_bycluster.txt",sep="\t",col.names=T)
write.table(Tcell.markers, file = "Tcell_dge_markers_byID.txt",sep="\t",col.names=T)
top50 <- Tcell.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
View(top50)
write.table(top50, file = "CD4_markers_byID.txt",sep="\t",col.names=T)
T_features <- unique(as.vector(top10$gene))
DotPlot(dge, features =T_features) + RotatedAxis()
ggsave(file="dot_plot_bypatient.eps",width = 70,height = 50,units = "cm")
ggsave(file="dot_plot_byID.eps",width = 50,height = 30,units = "cm")


dge <- ScaleData(object=dge,features=rownames(dge))
DoHeatmap(dge, features = c("CD3E","CD3D","CD3G","CD8A","CD8B","CTLA4","PDCD1","LAG3","TIGIT","CD14","CD68","CD163","FCER1A","KIT","FCER2","ENPP3","ENPP3","KIT","CPA3","LYZ","VCAN","KLK2","KLK3","ACPP"),raster = F) + theme(axis.text.y = element_text(size = 10))
ggsave(file="supervised_Heatmap.eps",width = 75,height = 75,units = "cm")
DoHeatmap(dge, features = as.character(top10$gene),combine = T,raster = T) + theme(axis.text.y = element_text(size = 5))
ggsave(file="Heatmap_bypatient.eps",width = 75,height = 75,units = "cm")
ggsave(file="Heatmap_byID.eps",width = 75,height = 75,units = "cm")

VlnPlot(dge,features = c("CD3E","CD3D","CD3G","CD8A","CD8B","CTLA4","PDCD1","LAG3","TIGIT","CD14","CD68","CD163","FCER1A","KIT","FCER2","ENPP3","ENPP3","KIT","CPA3","LYZ","VCAN","KLK2","KLK3","ACPP"))
ggsave(file="Markers.eps",width = 75,height = 75,units = "cm")

new.cluster.ids <- c("CD4_cluster1","CD8","Tumor_specific_CD8","CD8","CD8","CD8","CD8","CD4_cluster2","CD8_cluster2","Treg","KLK3+","NK-cell_cluster1","NKcell_cluster2")
names(new.cluster.ids) <- levels(dge)
dge <- RenameIdents(dge, new.cluster.ids)
DimPlot(dge, reduction = "umap",group.by = "ID",label=T)
ggsave(file="Umap_ID.eps",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "tsne",group.by = "ID",label=T)
ggsave(file="TSNE_ID.eps",width = 20,height = 20,units = "cm")


saveRDS(dge,"T-cell_dge.RDS")
dge$ID<-as.vector(dge@active.ident)

##### Let's look at CD8+ T-cells 
dge_Tcell_all<-dge
cell_CD8T<-colnames(dge_Tcell_all)[dge_Tcell_all$ID=="CD8"]
dge<-SubsetData(dge_Tcell_all,cells=cell_CD8T)
patient<-c("MAY_PB1","MAY_PB2","PR5186","PR5196","PR5199","AUG_PB1","PR5251_N","PR5251_T")
report<-NULL
for (i in 1:length(patient)) {
dge_temp<-SubsetData(dge,ident.use = patient[i])
Featurename<-c("Cytotoxic","Exhausted","PD1","EGFR")
Features<-list()
#Features[[1]]<-as.vector(as.character(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Cytotoxic.txt",header=FALSE,sep="\t"))))
#Features[[2]]<-as.vector(as.character(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Exhausted_Tcell.txt",header=FALSE,sep="\t"))))
Features[[1]] <- as.vector(as.character(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Cytotoxic.txt",header=FALSE,sep="\t"))))
Features[[2]]<-c("LAG3", "CTLA4", "CD244", "CD160", "ENTPD1", "TIGIT","PDCD1","HAVCR2")
Features[[3]]<-c("CD247","CD274","CD3D","CD3E","CD3G","CD4","CSK","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DRB1","HLA-DRB3","HLA-DRB5","LCK","PDCD1","PDCD1LG2","PTPN6")
Features[[4]]<-c("HSP90AA1","AREG","GAB1","CBL","HBEGF","SOS1","PIK3CA","PLCG1","EREG","KRAS","EGF","RPS27A","PIK3R1","EGFR","UBC","SHC1","TGFA","UBB","HRAS","BTC","GRB2","EPGN","NRAS","UBA52","STAT3")
for (j in 1:length(Features)) {
  gene.set <- Features[[j]][1:length(Features[[j]])]
  dge_temp <- AddModuleScore(object = dge_temp, features = list(gene.set), ctrl = 5, name = Featurename[j])
}
df<-as.data.frame(cbind(dge_temp$Cytotoxic1,dge_temp$Exhausted1))
colnames(df)<-c("Cytotoxic1","Exhausted1")

####linear regression
# reg1 <- lm(Cytotoxic1~Exhausted1,data=df) 
# summary(reg1)
# 
# with(df,plot(Cytotoxic1, Exhausted1))
# abline(reg1)

#### Lowess regression with outlier detection
setEPS()
postscript(paste0("Regression_",patient[i],".eps"))
reg1 <- lm(Cytotoxic1~Exhausted1,data=df) 
plot(Cytotoxic1~Exhausted1,data=df,pch=19,col="grey")#+abline(reg1,lwd=2)
lo <- loess.smooth(df$Cytotoxic1, df$Exhausted1)
lines(lo$x, lo$y, lwd = 2)
lines(lo$x, lo$y * 2, lwd = 1, col = 2)
lines(lo$x, lo$y/2, lwd = 1, col = 2)

f1 <- approxfun(lo$x, lo$y * 2)
(wh1 <- which(df$Cytotoxic1 >2*f1(df$Exhausted1)))
f2 <- approxfun(lo$x, lo$y /2)
(wh2 <- which(df$Cytotoxic1 <f2(df$Exhausted1)/2))

mt_exhausted <- df[c(wh1), ]
points(mt_exhausted$Exhausted1 ,mt_exhausted$Cytotoxic1, pch = 19, col = "red")

mt_cytotoxic <- df[c(wh2), ]
points(mt_cytotoxic$Exhausted1 ,mt_cytotoxic$Cytotoxic1, pch = 19, col = "green")
dev.off()

report_tmp<-c(reg1$coefficients[1],reg1$coefficients[2],length(wh1),length(wh2))
report<-rbind(report,report_tmp)
}
report<-cbind(patient,report)
colnames(report)<-c("Patient","intercept","K_exhausted","NCyto","NExhausted")
write.table(report,"regression_report.txt",col.names = T,sep="\t")
####### Finished regression#####



###### 
dge_backup<-dge
dge<-SubsetData(dge_backup,ident.use = c("PA_PR5186","PA_PR5196","PA_PR5199","PR5249_N","PR5251_N","PR5251_T","PR5254_N","PR5261_T","PR5269"))
dge$gleason<-"Gleason_7"
dge$gleason[dge$orig.ident %in% c("PA_PR5186","PA_PR5254_N","PA_PR5254_T")]="Gleason_9"
dge$LN="Neg"
dge$LN[dge$orig.ident %in% c("PA_PR5186","PR5249_N","PR5249_T")]="Pos"

dge$source<-as.vector(dge$orig.ident)
dge$source[dge$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")]="Paired_normal"
dge$source[dge$orig.ident %in% c("AUG_PB1", "PA_PB1", "PA_PR5199", "PR5249_T", "PR5251_T", "PR5269")]="ERGpos_Tumor"
dge$source[dge$orig.ident %in% c("PA_PB2","PA_PR5186","PA_PR5196","PR5254_T","PR5261_T")]="ERGneg_Tumor"

dge<-SetIdent(dge,value = as.vector(dge$source))
dge<-SetIdent(dge,value = as.vector(dge$ID))
dge<-SetIdent(dge,value=as.vector(dge$TorN))
dge<-SetIdent(dge,value = as.vector(dge$gleason))
dge<-SetIdent(dge,value = as.vector(dge$LN))

name_type=c("Paired_normal","ERGpos_Tumor","ERGneg_Tumor")
name_type=c("CD4","CD8","NK-cell","Treg")
name_type=c("Tumor","Paired_normal")
name_type=unique(dge$gleason)
name_type=unique(dge$LN)

Featurename<-c("Cytotoxic","Exhausted","Naive","PD1","EGFR","AR_independence","AR","Stem","TILN","Memory")
Features<-list()
Features[[1]] <- as.vector(as.character(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Cytotoxic.txt",header=FALSE,sep="\t"))))
Features[[2]]<-c("Exhausted","LAG3", "CTLA4", "CD244", "CD160", "ENTPD1", "TIGIT","PDCD1","HAVCR2")
Features[[3]]<-c("Naive","CCR7","TCF7","LEF1","SELL","IL7R","IL2RG")
Features[[4]]<-c("PD1","CD247","CD274","CD3D","CD3E","CD3G","CD4","CSK","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DRB1","HLA-DRB3","HLA-DRB5","LCK","PDCD1","PDCD1LG2","PTPN6")
Features[[5]]<-c("EGFR","HSP90AA1","AREG","GAB1","CBL","HBEGF","SOS1","PIK3CA","PLCG1","EREG","KRAS","EGF","RPS27A","PIK3R1","EGFR","UBC","SHC1","TGFA","UBB","HRAS","BTC","GRB2","EPGN","NRAS","UBA52","STAT3")
Features[[6]]<-as.vector(as.character(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/AR_independence.txt",header=FALSE,sep="\t"))))
Features[[7]]<-as.vector(as.character(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/AR.txt",header=FALSE,sep="\t"))))
Features[[8]]<-as.vector(as.character(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/Stemness.txt",header=FALSE,sep="\t"))))
Features[[9]]<-as.vector(as.character(unlist(read.csv("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/geneset/TILN.txt",header=FALSE,sep="\t"))))
Features[[10]]<-c("Memory","ITGB1","ITGA4","ITGA5","ITGAL","ITGAM","ITGB2","CD2","CD44","ICAM1","CD58")
for (i in 1:length(Features)) {
  gene.set <- Features[[i]][2:length(Features[[i]])]
  dge <- AddModuleScore(object = dge, features = list(gene.set), ctrl = 5, name = Featurename[i])
  V1<-as.vector(FetchData(object = dge,vars = paste0(Featurename[i],"1")))
  FeaturePlot(dge,features=paste0(c(Featurename[i]),"1"))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))
  ggsave(file=paste0("Modulescore_",Featurename[i],"_featureplot.eps"),width = 20,height = 20,units = "cm")
  
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(notch = T) + geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    #geom_signif(comparisons = list(c(name_type[1],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.05) +
    #geom_signif(comparisons = list(c(name_type[2],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.25) +
    #geom_signif(comparisons = list(c(name_type[1],name_type[2])),map_signif_level=TRUE,y_position = max(V1)+0.45) +
    ggtitle(Featurename[i])
  ggsave(file=paste0("Modulescore_",Featurename[i],"_boxplot.eps"),width = 50,height = 20,units = "cm")
  
}
dge$Cytotoxic1[dge$ID %in% c("NK-cell_cluster1","NKcell_cluster2")]=+1
dge<-SubsetData(dge_Tcell,ident.use = c("Tumor_specific_CD8"))
dge$TorN="Tumor"
dge$TorN[dge$orig.ident %in% c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")]="Paired_normal"
dge<-SetIdent(dge,value = as.vector(dge$TorN))
dge<-SetIdent(dge,value = as.vector(dge$source))
dge<-SetIdent(dge,value = as.vector(dge$ID))
marker_list=c("CD274","PDCD1","CD8A","PRF1","FOXP3","GNLY","TIGIT","LAG3","CTLA4","MKI67","IDO1")
name_type=c("Paired_normal","ERGpos_Tumor","ERGneg_Tumor")
name_type=c("Tumor","Paired_normal")

for (i in 1:length(marker_list)) {
  V1<-as.vector(FetchData(object = dge,vars = marker_list[i],slot = "data"))
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  # ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(notch = F) + geom_jitter(shape=16,position = position_jitter(0.1))+geom_signif(comparisons = list(c("Club","BE")),map_signif_level=TRUE,y_position = max(V1)+0.5)+
  #    geom_signif(comparisons = list(c("Club","LE")),map_signif_level=TRUE,y_position = max(V1)+1) +stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
  #    ggtitle(marker_list[i])
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(notch = F) + geom_jitter(shape=16,position = position_jitter(0.1))+
    # geom_signif(comparisons = list(c(name_type[1],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.05) +
    # geom_signif(comparisons = list(c(name_type[2],name_type[3])),map_signif_level=TRUE,y_position = max(V1)+0.25) +
     geom_signif(comparisons = list(c(name_type[1],name_type[2])),map_signif_level=TRUE,y_position = max(V1)+0.45) +
    # 
    stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
    ggtitle(marker_list[i])
  ggsave(file=paste0(marker_list[i],"_marker_boxplot.eps"),width = 20,height = 20,units = "cm")
  
}


cell_exhausted<-colnames(dge)[dge$Exhausted1>mean(dge$Exhausted1)+sd(dge$Exhausted1)]
dge_exhausted<-SubsetData(dge,cells = cell_exhausted)
dge$exhaust="Others"
dge$exhaust[colnames(dge) %in% cell_exhausted]="Exhausted"
cell_Cytotoxic<-colnames(dge)[dge$Cytotoxic1>mean(dge$Cytotoxic1)+sd(dge$Cytotoxic1)]
dge_Cytotoxic<-SubsetData(dge,cells = cell_Cytotoxic)
dge$cytotoxic="Others"
dge$cytotoxic[colnames(dge) %in% cell_Cytotoxic]="Cytotoxic"
DimPlot(dge,cells.highlight = cell_exhausted)
DimPlot(dge,cells.highlight = cell_Cytotoxic)

dge$ID<-as.vector(dge@active.ident)
dge<-SetIdent(dge,value = as.vector(dge$exhaust))
dge<-SetIdent(dge,value = as.vector(dge$cytotoxic))

tmp<-read.csv("../../geneset/androgen.gmt",header=FALSE,sep="\t")
for (i in 1:nrow(tmp)) {
  
  name_feature=as.character(tmp[i,1])
  feature_tmp<-na.omit(tmp[i,])
  
  Features<-as.vector(feature_tmp[2:ncol(feature_tmp)])
  dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature)
  #names(x = dge[[]])
  
  #str<-paste0(name_feature,"1")
  #dge[[str]]<-range01(dge[[str]])
  
  #FeaturePlot(object = dge, features =paste0(name_feature,"1") ,label = F)
  #ggsave(file=paste0("BE_",name_feature,".eps"),width = 20,height = 20,units = "cm")
  V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  ggplot(BOX_df, aes(id,value)) + geom_boxplot() + geom_signif(comparisons = split(t(combn(levels(BOX_df$id), 2)), seq(nrow(t(combn(levels(BOX_df$id), 2))))),textsize=6,map_signif_level=TRUE)
  ggsave(file=paste0("BE_",name_feature,"_boxplot.eps"),width = 20,height = 20,units = "cm")
  #VlnPlot(object = dge, features =paste0(name_feature,"1"))+stat_summary(fun = mean, geom='point', size = 20, colour = "blue",shape=95)
  #ggsave(file=paste0("BE_",name_feature,"_vlnplot.eps"),width = 20,height = 20,units = "cm")
  
}

name_feature="PD1"
PD1_features=c("CD247","CD274","CD3D","CD3E","CD3G","CD4","CSK","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DRB1","HLA-DRB3","HLA-DRB5","LCK","PDCD1","PDCD1LG2","PTPN6")
Features<-PD1_features
dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature)
V1<-as.vector(FetchData(object = dge,vars = "PCAT4"))
BOX_df<-NULL
BOX_df$id<-dge@active.ident
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value)) + geom_boxplot() + geom_signif(comparisons = split(t(combn(levels(BOX_df$id), 2)), seq(nrow(t(combn(levels(BOX_df$id), 2))))),textsize=6,map_signif_level=TRUE)
ggsave(file=paste0("PD1_",name_feature,"_boxplot.eps"),width = 20,height = 20,units = "cm")

name_feature="TILN"
PD1_features=c("CD247","CD274","CD3D","CD3E","CD3G","CD4","CSK","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DRB1","HLA-DRB3","HLA-DRB5","LCK","PDCD1","PDCD1LG2","PTPN6")
Features<-PD1_features
dge <- AddModuleScore(object = dge, features = list(Features), name = name_feature)
dge<-SetIdent(object = dge,value = dge_Tcell_withexhaust$Idents)
V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
BOX_df<-NULL
BOX_df$id<-dge@active.ident
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value)) + geom_boxplot() + geom_signif(comparisons = split(t(combn(levels(BOX_df$id), 2)), seq(nrow(t(combn(levels(BOX_df$id), 2))))),textsize=6,map_signif_level=TRUE)
ggsave(file=paste0("Tcell_exhausted_",name_feature,"_boxplot.eps"),width = 20,height = 20,units = "cm")



dge<-SubsetData(dge_Tcell,ident.use = c("CD8","Tumor_specific_CD8"))

source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
Data<-dge_temp
name1<-"Ident"
name2<-"Fibroblast_ERGtype"
Prepare.for.GSEA(Data,name1,name2)


dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 2000)
dge<- ScaleData(object=dge,features=VariableFeatures(object = dge))

all_markers <-FindAllMarkers(dge,min.diff.pct=0.00001,min.pct =  0.00001,test.use = "wilcox")

write.table(all_markers,"NKcluster2_NKcluster2.txt",sep="\t",row.names=F, col.names = T)
all_markers<-all_markers[1:(nrow(all_markers)/2),]
DE_select<-all_markers$p_val_adj<=10e-05
length(which(DE_select))
all_markers$threshold <- DE_select
ggplot(all_markers) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj), colour=threshold)) +
  ggtitle("Club vs LE overexpression") +
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

ggplot(all_markers_reorder) +
  geom_point(aes(x = avg_logFC, y = -log10(p_val_adj), colour = threshold)) +
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = ifelse(genelabels !="", genelabels,""))) +
  ggtitle("NKcluster2_NKcluster1 T-cell") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave(file="NKcluster2_NKcluster2_Volcano.eps",width = 20,height = 20,units = "cm")


Featurename<-c("Cytotoxic","Exhausted","Naive","PD1","EGFR","AR_independence")
Featurename<-c("PDCD1","TIGIT","LAG3","HAVCR2","CTLA4","CD274","MKI67")
Featurename<-c("Exhausted1","Cytotoxic1","Naive1","PD11","EGFR1","AR_independence1","AR1","Stem1")
for (i in 1:length(Featurename)) {
  V1<-as.vector(FetchData(object = dge,vars = Featurename[i],slot = "data"))
  FeaturePlot(dge,features=c(Featurename[i]))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(min(V1), max(V1)))
  ggsave(file=paste0("../Tcell_analysis/Featureplot_",Featurename[i],".eps"),width = 20,height = 20,units = "cm")
}



dge<-dge_Tcell
dge<-SetIdent(dge,value = as.vector(dge$orig.ident))
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- ScaleData(object=dge,features=rownames(dge))
cluster.averages <- AverageExpression(dge,assays = "RNA",return.seurat = F,)

#CellScatter(cluster.averages, cell1 = "NK", cell2 = "CD8_T")
Feature_immune<-NULL
Feature_immune$names=c("IDO1","NKG7","NGLY","CD8B","B2M","C11orf45","LGMN","GPNMB","HS3ST2","CYBB","CD68","HAVCR2","LAIR1","LGALS9","TNFSF8","CD58","CLEC5A","FUCA1","AHR","TNFSF15","PVRL3","TNFSF10","C10orf54","ICOSLG","IRF8","LILRA4","GPR146","SELP","PLD4","PHEX","IFIT2","IFIT1","IFIT3","RSAD2","MX1","MX2","C15orf53","CD226","CD40LG","CLEC4C","GPR15","CD28","TNFSF18","FOXP3","CD27","CD2","BTLA","SLAMF1","CTLA4","ICOS","TIGIT","IL32","CXCR3","PDCD1","GZMA","CD8A","LAG3","TAP1","GZMB","PRF1","CD70","TNFRSF9","MMP9","ISG20","TNFSF14","TNFSF4","CD160","CD244","CD274","TNFSF9","CD40","HLA-A","PTCRA","TM4SF19","IL3RA","TNFRSF8","IRF7","TNFRSF4","TNFRSF18","TNFRSF25","CCR5","CCR7","CCR9","CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR7","CXCL8")
Feature_immune$type=c("Inhibition","Cyto","Cyto","CD8","DC_HLA","DC_HLA","Macro","DC_HLA","Macro","Macro","Macro","IFN","IFN","IFN","IFN","IFN","Macro","Macro","Inhibition","IFN","IFN","Inhibition","IFN","IFN","DC_HLA","DC_HLA","Inhibition","Inhibition","DC_HLA","DC_HLA","Inhibition","Inhibition","Inhibition","Inhibition","Inhibition","Inhibition","Treg","IFN","IFN","DC_HLA","Treg","IFN","IFN","Treg","IFN","IFN","IFN","IFN","IFN","IFN","IFN","Treg","DC_HLA","IFN","Cyto","CD8","IFN","DC_HLA","DC_HLA","Cyto","IFN","IFN","Macro","Inhibition","IFN","IFN","IFN","IFN","IFN","IFN","IFN","DC_HLA","DC_HLA","Macro","DC_HLA","IFN","DC_HLA","IFN","IFN","IFN","Chemokine","Chemokine","Chemokine","Chemokine","Chemokine","Chemokine","Chemokine","Chemokine","Chemokine","Chemokine")
Type<-Feature_immune$type[order(Feature_immune$type)]
Names<-Feature_immune$names[order(Feature_immune$type)]
Feature_immune$names=as.vector(Names)
Feature_immune$type=as.vector(Type)
Feature_immune<-as.data.frame(Feature_immune)
Feature_immune<-Feature_immune[Feature_immune$type %in% c("NK-cell","CD8","Cyto","Chemokine","IFN","Inhibition","Treg"),]
mtx<-cluster.averages$RNA
mtx<-mtx[rownames(mtx) %in% Feature_immune$names,]
#mtx<-log2(mtx+1)
mtx<-mtx[Feature_immune$names,]

mtx_backup<-mtx
mtx<-na.omit(mtx[order(match(rownames(mtx),Feature_immune$names)),])
mtx<-mtx[,c(7,8,10,9,11,12,1,2,3,4,5,6,13)]
mtx_norm<-na.omit(t(apply(mtx, 1, function(x)(x-min(x))/(max(x)-min(x)))))
metadata<-data.frame(as.character(Feature_immune$type),row.names = Feature_immune$names)
colnames(metadata)<-"Marker_Type"
require(pheatmap)
pheatmap(mtx_norm,cluster_rows = F,cluster_cols = F,filename = "T-cell_supervisedheatmap.pdf",annotation_row = metadata,width = 10,height = 10)

saveRDS(dge,"T-cell_dge.RDS")

dge_tam$ID<-as.vector(dge_tam$orig.ident)
dge_tam$ID[dge_tam$orig.ident %in% sample_tumor]="Tumor_tam"
dge_tam$ID[dge_tam$orig.ident %in% sample_pairedN]="Tumor_pairedN"
dge_tam<-SetIdent(dge_tam,value = as.vector(dge_tam$ID))

source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
Data<-dge_tam
name1<-"Ident"
name2<-"tam_compare"
Prepare.for.GSEA(Data,name1,name2)



library(dplyr)
library(tidyr)
library(broom)
raw_E_score<-read.csv("Immune_dge_ssGSEA.PROJ.gct",sep="\t",header = F)
cellname<-na.omit(raw_E_score[1,3:ncol(raw_E_score)])
raw_E_score<-raw_E_score[-c(1),]
raw_E_score<-raw_E_score[,-c(2)]

raw_E_score<-t(raw_E_score)
raw_E_cluster<-read.csv("Immune_dge_phenotypes.cls",sep=" ",header=F)
raw_E_score_df<-as.data.frame(raw_E_score[2:2536,1:ncol(raw_E_score)])

rownames(raw_E_score_df)<-as.vector(t(cellname))
colnames(raw_E_score_df)<-raw_E_score[1,1:ncol(raw_E_score)]
raw_E_score_df <- as.data.frame(sapply(raw_E_score_df, as.numeric))
raw_E_score_df$CLUSTER<-t(raw_E_cluster)


varNames <- names(raw_E_score_df)[-length(raw_E_score_df)]
mean <- colMeans(raw_E_score_df[,1:length(raw_E_score_df)-1])
master <- list()
i <- 1
groupname=unique(raw_E_score_df$CLUSTER)
raw_E_score_df[is.na(raw_E_score_df)] = 0
## main for loop subsets; lapply calculates t-statistics for all variables in the subset
for (group in unique(raw_E_score_df$CLUSTER)){
  # create a list of t-test result in a given "group" subset
  results <- lapply((1:length(varNames)), FUN = function(x, subset = raw_E_score_df[raw_E_score_df$CLUSTER == group,]) {
    t.test(subset[varNames[x]], mu = mean[x], alternative = "two.sided")
  })
  master[[i]] <- results
  i <- i + 1
  print(paste0("Cluster finished for ", groupname[i-1]),quote=F)
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
  write.table(results_tmp,paste0(groupname[i],"_meanscores.txt"),sep="\t",row.names = T,col.names = T,quote=F)
}


for (i in 1:length(master)) {
  tmp<-read.csv(paste0(groupname[i],"_meanscores.txt"),sep="\t",header = T)
  tmp_sig<-tmp[tmp$FDR<0.01,]
  tmp_sig<-tmp_sig[order(tmp_sig$Mean,decreasing = T),]
  tmp_sig_top30<-as.data.frame(tmp_sig[1:50,])
  tmp_sig_top30$GSEA<-as.vector(rownames(tmp_sig_top30))
  ggplot(tmp_sig_top30) +
    geom_bar( aes(x=reorder(GSEA,Mean), y=Mean), stat="identity", fill="forestgreen", alpha=0.5) +
    geom_errorbar( aes(x=GSEA, ymin=X95.CI_low, ymax=tmp_sig_top30$X95.CI_hi), width=0.4, colour="grey", alpha=0.9, size=1.5) +
    ggtitle(paste0("Top50 ssGSEA ES for Cluster_",groupname[i]))+coord_flip()
  ggsave(file=paste0("Top50_score_Cluster_",groupname[i],".pdf"),width = 40,height = 30,units = "cm")
  write.table(tmp_sig,paste0("Cluster_",groupname[i],"_Sigscores.txt"),col.names = T,row.names = T,sep="\t",quote=F)
}



dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- ScaleData(object=dge,features=rownames(dge))
cluster.averages <- AverageExpression(dge,assays = "RNA",return.seurat = F)
mtx<-cluster.averages$RNA

Feature_exhaustion=c("TIGIT","CTLA4","CD244","CD160","PDCD1","TBC1D4","ENTPD1") #ENTPD1 is CD39
Feature_cytotoxic<-c("ITGAL","ITGB2","CD247","GZMH")

Feature<-NULL
Feature$names=c(Feature_exhaustion,Feature_cytotoxic)
Feature$type=c(rep("Exhaustion",length(Feature_exhaustion)),rep("Cytotoxicity",length(Feature_cytotoxic)))

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
mtx<-mtx[,c(2,1,3,9,5,7,4,8)]
#for club, use this:
metadata<-data.frame(as.character(Feature$type),row.names = Feature$names)
colnames(metadata)<-"Marker_Type"
require(pheatmap)
pheatmap(mtx,cluster_rows = F,cluster_cols = F,filename = "Tcell_supervisedexhaustion.pdf",annotation_row = metadata,width = 10,height = 10)


###ERGpos_patient "PR5249","PR5251","PR5269","PR5199"
###ERGneg_patient "PR5186","PR5196","PR5254","PR5261"
DimPlot(dge_Tcell)
dge<-SubsetData(dge_Tcell,cells = colnames(dge_Tcell)[dge_Tcell$orig.ident %in% c("PR5249_N","PR5251_N",
                                                                                  "PR5186","PR5196",
                                                                                  "PR5199","PR5251_T",
                                                                                  "PR5254_N","PR5261_T",
                                                                                  "PR5269","PR5249_T",
                                                                                  "PR5254_T","PR5261_N")])
dge$ID<-as.vector(dge@active.ident)
dge$ID[dge$ID %in% c("CD4_cluster1","CD4_cluster2")]="CD4"
dge$ID[dge$ID %in% c("CD8_cluster2","CD8_cluster3","CD8")]="CD8"
dge$ID[dge$ID %in% c("NK-cell_cluster1","NKcell_cluster2")]="NKcell"

dge$ERGexpression="POS"
dge$ERGexpression[dge$orig.ident %in% c("PR5186","PR5196","PR5254_N","PR5254_T","PR5261_T","PR5261_N")]="NEG"

dge<-SetIdent(dge,value = as.vector(dge$ID))
group=c("CD4","CD8","NKcell","Treg","All_Tcell")


i=5
dge_temp<-SubsetData(dge,ident.use = group[i])
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$ERGexpression))

source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
Data<-SubsetData(dge_E,ident.use = c("ERGpos_Tumor","ERGneg_Tumor"))
name1<-"Ident"
name2<-"Compare_ERG"
Prepare.for.GSEA(Data,name1,name2)

#NKcell<-SubsetData(dge_Tcell,ident.use = c("NKcell_cluster2","NK-cell_cluster1"))
all_markers <-FindAllMarkers(dge_temp,min.diff.pct=0.001,min.pct =  0.001,test.use = "wilcox")
write.table(all_markers,paste0(group[i],"_markers.txt"),sep="\t",row.names=F, col.names = T)
all_markers<-all_markers[1:(nrow(all_markers)/2),]
DE_select<-all_markers$p_val_adj<=10e-05
length(which(DE_select))
all_markers$threshold <- DE_select
ggplot(all_markers) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj), colour=threshold)) +
  ggtitle("ERGpos ERGneg overexpression") +
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

ggplot(all_markers_reorder) +
  geom_point(aes(x = avg_logFC, y = -log10(p_val_adj), colour = threshold)) +
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = ifelse(genelabels !="", genelabels,""))) +
  ggtitle("ERGpos ERGneg overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave(file=paste0(group[i],"_Volcano.eps"),width = 20,height = 20,units = "cm")


dge$ERGtype<-"ERGneg"
dge$ERGtype[dge$orig.ident %in% c("AUG_PB1","MAY_PB1","PR5249","PR5251","PR5199","PR5269")]="ERGpos"

DimPlot(dge,group.by = "ERGtype")
ggsave(file="Myeloid_byERGtype.eps",width = 20,height = 20,units = "cm")
table(dge$ERGtype,dge@active.ident)




cell_M_MDSC<-intersect(colnames(dge_myeloid),colnames(dge_M_MDSC))
DimPlot(dge_myeloid,cells.highlight = intersect(colnames(dge_myeloid)[dge_myeloid$ID %in% c("CD14+ Monocyte","CD14- Monocyte")],cell_M_MDSC))
ggsave(file="../Myeloid_analysis/M_MDSC_highlight.eps",width = 20,height = 20,units = "cm")
cell_PMC_MDSC<-intersect(colnames(dge_myeloid),colnames(dge_PMN_MDSC))
DimPlot(dge_myeloid,cells.highlight = intersect(colnames(dge_myeloid)[dge_myeloid$ID %in% c("CD14+ Monocyte","CD14- Monocyte")],cell_PMC_MDSC))
ggsave(file="../Myeloid_analysis/PMN_MDSC_highlight.eps",width = 20,height = 20,units = "cm")



###Generate signature
for (i in 1:length(unique(dge@active.ident))) {
  markers_temp <- FindMarkers(dge, ident.1 = unique(dge@active.ident)[i], ident.2 = NULL, only.pos = TRUE)
  markers_temp$gene<-as.vector(rownames(markers_temp))
  markers_temp<-markers_temp[order(markers_temp$avg_logFC,decreasing = T),]
  sig_cluster<-markers_temp$gene[markers_temp$p_val_adj<0.05]
  sig_cluster<-sig_cluster[1:min(length(sig_cluster),50)]
  write.table(sig_cluster,paste0("",unique(dge@active.ident)[i],"_Signature.txt"),sep = "\t",quote = F,col.names = F,row.names = F)
}
Another_Exhausted_Signature<-c("PDCD1","CTLA4","HAVCR2","LAG3", "BTLA", "CD244", "CD160", "ENTPD1", "VSIR", "TIGIT")
dge<-AddModuleScore(dge,features = list(Another_Exhausted_Signature),name = "Another_Exhausted")
VlnPlot(dge,features=c("Another_Exhausted1"),group.by = "ERGtype")
ggsave(file="CD8_Anotherexhausted.pdf",width = 20,height = 20,units = "cm")
