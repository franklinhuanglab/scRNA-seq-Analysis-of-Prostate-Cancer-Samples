setwd("~/Desktop/hepatoblastoma/Epithelial/")
mkdir("./Tumor")
setwd("./Tumor")

dge<-SubsetData(HBO_epithelial,cells = colnames(HBO_epithelial)[HBO_epithelial$ID_final %in% c("Tumor_Hepatocyte","DCN_high_Tumor","Neuroendocrine")])
dge$sampletype<-"Normal"
dge$sampletype[dge$patient %in% c("Patient1_T","Patient2_T","Patient3_T",
                                  "Patient4_T","Patient5_T","Patient6_T",
                                  "Patient7_T","Patient8_T","Patient9_T")]="Tumor"

dge<-HBO_final
unique(dge$patient)
dge$patient[dge$patient %in% c("Patient1_H","Patient1_N")]="Patient1_N"
dge$patient[dge$patient %in% c("Patient2_H","Patient2_N")]="Patient2_N"
dge$patient[dge$patient %in% c("Patient3_H","Patient3_N")]="Patient3_N"
dge$patient[dge$patient %in% c("Patient4_H","Patient4_N")]="Patient4_N"
dge$patient[dge$patient %in% c("Patient5_H","Patient5_N")]="Patient5_N"
dge$patient[dge$patient %in% c("Patient6_H","Patient6_N")]="Patient6_N"
dge$patient[dge$patient %in% c("Patient7_H","Patient7_N")]="Patient7_N"
dge$patient[dge$patient %in% c("Patient8_H","Patient8_N")]="Patient8_N"
dge$patient[dge$patient %in% c("Patient9_H","Patient9_N")]="Patient9_N"

write.table(dfr,"./Table_Annotation_new.txt",sep="\t",quote = F,col.names = T,row.names = T)
dge = SubsetData(HBO_final,ident.use = c("Normal_Hepatocyte","Tumor_Hepatocyte","BECs","Neuroendocrine",
                                         "DCN_high_Tumor"))
dge<-SubsetData(dge,ident.use = c("Club"))
dge<-SubsetData(dge,ident.use = c("BE_Intermediate"))
dge<-SubsetData(dge,ident.use = c("LE_Club"))
dge <- CreateSeuratObject(counts = dg_temp, project = "2wk_data", min.features = 0)
dge<-dge_CD8
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = dge)
dge <- ScaleData(object = dge, features = VariableFeatures(dge))#, vars.to.regress = "percent.mt")
dge <- RunPCA(dge, features = VariableFeatures(object = dge),npcs = 100)


slot(dge[["pca"]], "misc")
print(x = dge[["pca"]], dims = 1:5, nfeatures = 5)
mat <- Seurat::GetAssayData(dge,assay="RNA",slot="scale.data")
pca <-dge[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)

dge$ID_2[dge$ID %in% c("BE_0","BE_1","BE_2","BE_3","BE_5","BE_6","BE_7")]="BE"
dge$ID_2[dge$ID %in% c("LE_0","LE_1","LE_2","LE_4","LE_5","LE_6")]="LE"
dge$ID_2[dge$ID %in% c("Club_0","Club_1","Club_2","Club_3","Club_4","Club_5")]="Club"

DimPlot(object = dge, reduction = "pca")
ggsave(file="PCA.pdf",width = 30,height = 30,units = "cm")

DimPlot(object = dge, reduction = "pca",group.by = "ID_2")
ggsave(file="PCA_byID.pdf",width = 30,height = 30,units = "cm")

pdf("PCA_heatmap.png")
DimHeatmap(object = dge, dims = 1:10, cells = nrow(dge), balanced = TRUE)
dev.off()
png("Elbowplot_dge.png")
ElbowPlot(dge,ndims = 100)
dev.off()

dge <- JackStraw(dge, num.replicate = 100,dims = 50)
dge <- ScoreJackStraw(dge, dims = 1:50)
JackStrawPlot(dge, dims = 1:50)
ggsave(file="Strawplot.pdf",width = 30,height = 30,units = "cm")

#TODO: Set resolution
n_pc <- 46
#Sawyer data Club cell
n_pc <- 75
resolut=0.2
#Chen data Club cell
n_pc <- 21
resolut=0.2
#Dong data Club cell
n_pc <- 14
resolut=0.4
#CD8 Tcell
n_pc <- 13
resolut=0.3


dge <- FindNeighbors(dge, reduction = "pca", dims = 1:n_pc)
dge <- FindClusters(dge, resolution =c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5))


pdf("Cluster_stability.pdf",width = 20,height = 25)
clustree(dge, prefix = "RNA_snn_res.")
dev.off()

dge <- FindNeighbors(dge, dims = 1:n_pc)
###for tumor cell use 0.6, all epithelial 0.3
dge <- FindClusters(dge, resolution = resolut, random.seed = 10)
dge <- RunUMAP(dge, dims = 1:n_pc, seed.use = 10)

DimPlot(dge,label = T)+ NoLegend()
ggsave(file="Umap_raw.pdf",width = 20,height = 20,units = "cm")
DimPlot(dge,label = T,group.by = "ID")+ NoLegend()
ggsave(file="Umap_initial_annotation.pdf",width = 20,height = 20,units = "cm")
DimPlot(dge,group.by = "ID_2",label = F)
ggsave(file="Umap_byID_2.pdf",width = 30,height = 20,units = "cm")


VlnPlot(dge,features=c("MKI67","TOP2A","PTPRC","EPCAM"),ncol=2)
ggsave(file="Vlnplot_MKI67_EPCAM_annotataed.pdf",width = 40,height = 20,units = "cm")
library(SDMTools)
library(Seurat)
plot = DimPlot(dge, label = T, repel = T, label.size = 3)
select.cells <- CellSelector(plot = plot)
DimPlot(dge, label = T, cells.highlight = select.cells,repel = T, label.size = 3)
dge@meta.data$annotation_789 = as.character(dge@meta.data$annotation_789)
dge@meta.data[select.cells, "annotation_789"] = "Stellate/smooth muscle"
dge = SetIdent(dge, value = "annotation_789")

# Cleaning up some mislabelled stellate cells
stellate = rownames(dge@meta.data[which(dge@active.ident == "Tumor cells" & dge@meta.data$type == "Normal"),])

dge <- SetIdent(object = dge, cells = stellate, value = 'Stellate/smooth muscle')

DotPlot(dge, cols = c("grey","red"), features = c("FCGR3B","REG1A","GPC3", "DLK1","CD3E","KLRB1",
                                                  'CLEC10A',"TF","LYZ","C1QA","KRT19", "FLT1","CD79A",
                                                  "MS4A2","ACTA2", "HBB", "LILRA4")) + 
  scale_y_discrete(limits= c("pDCs","Erythrocytes", "Fibroblasts", "Basophils", "B cells", 
                             "Endothelial cells", "BECs", "Macrophages", "Monocytes", "Hepatocytes",  
                             "Dendritic cells", "NK/T/NKT", "Tumor cells", "Neutrophils")) + 
  theme(axis.text.x = element_text(angle = 325, size = 10, hjust=0.1,vjust=0.1),
        axis.text.y = element_text(size = 10), axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + coord_flip()


library(RColorBrewer)
colors = c(brewer.pal(12, "Paired"), "aquamarine", "violet")
DimPlot(dge, label = T, repel = T, label.size = 3, cols = colors) + 
  annotate(geom="text", x=10, y=13, size = 5, label="26533 cells") +
  annotate("segment", x = -13, xend = -13, y = -14, yend = -10, colour = "black",
           size=0.5, alpha=1, arrow=arrow(type = "open", angle = 18,length = unit(0.1, "inches"))) +
  annotate("segment", x = -13, xend = -10, y = -14, yend = -14, colour = "black",
           size=0.5, alpha=1, arrow=arrow(type = "open", angle = 18,length = unit(0.1, "inches"))) +
  annotate("text", x = -14, y = -12.5, label = "UMAP_2", size = 3, angle = 90) +
  annotate("text", x = -11.5, y = -15, label = "UMAP_1", size = 3) & 
  NoAxes() 

library(dplyr)
markers = FindAllMarkers(dge_old)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)  %>% 
  arrange(desc(avg_logFC), .by_group = TRUE)
top10
write.table(top10,"./Previous_annotation_top10.txt",sep="\t",col.names = T)


dge_temp<-dge
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$patient))
dge_temp<-SubsetData(dge_temp,max.cells.per.ident = 100)
av.exp <- AverageExpression(dge_temp)$RNA
av.exp<-av.exp[rownames(av.exp) %in% VariableFeatures(dge),]
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, as.vector(unique(dge_temp@active.ident)))
ggplot(cor.df, aes(x, y),lab_size = 5) +geom_tile(aes(fill = correlation))+
  scale_fill_gradientn(colours = c("blue","cyan","lightcyan","red")) + RotatedAxis()+geom_text(aes(label = sprintf("%1.2f",correlation)))
ggsave(file="Correlation_bycluster.pdf",width = 30,height = 30,units = "cm")


library(clustifyr)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
# Matrix of normalized single-cell RNA-seq counts
pbmc_matrix <- dge_allpatient@assays$RNA@data

# meta.data table containing cluster assignments for each cell
# The table that we are using also contains the known cell identities in the "classified" column
pbmc_meta <- dge_allpatient@meta.data
vargenes <- VariableFeatures(dge_allpatient)

ref_matrix<-Pliver@assays$RNA@data
ref_meta<-Pliver@meta.data
Pliver$ID<-as.vector(Pliver@active.ident)
new_ref_matrix <- average_clusters(
  mat = ref_matrix,
  metadata = ref_meta$ID, # or use metadata = pbmc_meta, cluster_col = "classified"
  if_log = TRUE                    # whether the expression matrix is already log transformed
)

head(new_ref_matrix)


res <- clustify(
  input = pbmc_matrix, # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
  metadata = pbmc_meta, # meta.data table containing cell clusters
  ref_mat = new_ref_matrix,
  cluster_col = "ID", # name of column in meta.data containing cell clusters
  query_genes = vargenes # list of highly varible genes identified with Seurat
)

# Peek at correlation matrix
res[1:5, 1:5]
res2 <- cor_to_call(
  cor_mat = res,                  # matrix correlation coefficients
  cluster_col = "ID" # name of column in meta.data containing cell clusters
)
res2[1:5, ]

pbmc_meta2 <- call_to_metadata(
  res = res2,                     # data.frame of called cell type for each cluster
  metadata = pbmc_meta,           # original meta.data table containing cell clusters
  cluster_col = "ID" # name of column in meta.data containing cell clusters
)

plot_cor_heatmap(cor_mat = res,col = c("blue", "white", "red")) 

dge$Sampletype<-dge$patient
dge$Sampletype[!(dge$patient %in% c("Patient1_T","Patient2_T","Patient3_T","Patient4_T","Patient5_T","Patient6_T","Patient7_T","Patient8_T","Patient9_T"))]="Normal"

DimPlot(dge,label = T,group.by = "Sampletype")
ggsave(file="Umap_initial_bysampletype.pdf",width = 30,height = 30,units = "cm")


dge_temp<-SubsetData(dge,cells = colnames(dge)[dge$ID %in% c("BECs","Hepatocyte_2")])
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$ID))
dge_temp$patientID<-dge_temp$patient
dge_temp$patientID[dge_temp$patient %in% c("Patient1_H","Patient1_N","Patient1_T")]="Patient1"
dge_temp$patientID[dge_temp$patient %in% c("Patient2_H","Patient2_N","Patient2_T")]="Patient2"
dge_temp$patientID[dge_temp$patient %in% c("Patient3_H","Patient3_N","Patient3_T")]="Patient3"
dge_temp$patientID[dge_temp$patient %in% c("Patient4_H","Patient4_N","Patient4_T")]="Patient4"
dge_temp$patientID[dge_temp$patient %in% c("Patient5_H","Patient5_N","Patient5_T")]="Patient5"
dge_temp$patientID[dge_temp$patient %in% c("Patient6_H","Patient6_N","Patient6_T")]="Patient6"
dge_temp$patientID[dge_temp$patient %in% c("Patient7_H","Patient7_N","Patient7_T")]="Patient7"
dge_temp$patientID[dge_temp$patient %in% c("Patient8_H","Patient8_N","Patient8_T")]="Patient8"
dge_temp$patientID[dge_temp$patient %in% c("Patient9_H","Patient9_N","Patient9_T")]="Patient9"
dge_temp$ID_new<-dge_temp$patientID
dge_temp$ID_new[dge_temp$ID =="BECs"]="BECs"
dge_temp$ID_new[dge_temp$ID =="Hepatocyte_2"]="Hepatocyte_2"
table(dge_temp$ID_new)

dge_1<-club_PCA
dge_1$ID<-paste0("Club_",dge_1$seurat_clusters)
dge_1<-SetIdent(dge_1,value = as.vector(dge_1$ID))
dge_2<-BE_PCA
dge_2$ID<-paste0("BE_",dge_2$seurat_clusters)
dge_2<-SetIdent(dge_2,value = as.vector(dge_2$ID))
dge_2<-SubsetData(dge_2,ident.remove = c("BE_4","BE_8"))
dge_2<-SetIdent(dge_2,value = as.vector(dge_2@active.ident))
dge_3<-LE_PCA
dge_3$ID<-paste0("LE_",dge_3$seurat_clusters)
dge_3<-SetIdent(dge_3,value = as.vector(dge_3$ID))
dge_3<-SubsetData(dge_3,ident.remove = c("LE_3","LE_7"))
dge_3<-SetIdent(dge_3,value = as.vector(dge_3@active.ident))
dge_temp<-merge(dge_1,y=dge_2)
dge_temp<-merge(dge_temp,y=dge_3)

dge_infer_mtx <- as.matrix(GetAssayData(dge_temp, slot = "counts"))
write.table(dge_infer_mtx,"Nontumor_Epithelial_counts.txt",sep="\t",row.names = T,col.names = T,quote=F)
phenotype<-NULL
phenotype$cellname<-colnames(dge_temp)
phenotype$ID<-as.vector(dge_temp@active.ident)
phenotype<-as.data.frame(phenotype)
write.table(phenotype,"Nontumoe_Epithelial_phenotype.txt",sep="\t",row.names = F,col.names = T,quote=F)


dge<-readRDS("./dge_epithelial.rds")
dge_infer_mtx <- as.matrix(GetAssayData(dge, slot = "counts"))
write.table(dge_infer_mtx,"Epithelial_counts.txt",sep="\t",row.names = T,col.names = T,quote=F)

source("/Users/hsong/Desktop/Seqwell_0913/human_output/prostate/GSEA_prepare.R")
Data<-SubsetData(HBO_epithelial,ident.use = c("Hepatocyte_1","Hepatocyte_2","Hepatocyte_3","Hepatocyte_4","Hepatocyte_5"))
Data<-SetIdent(Data,value = as.vector(Data@active.ident))
name1<-"Ident"
name2<-"HBO_Hepatocyte"
Prepare.for.GSEA(Data,name1,name2)


Data$Sampletype<-Data$patient
Data$Sampletype[!(Data$patient %in% c("Patient1_T","Patient2_T","Patient3_T","Patient4_T","Patient5_T","Patient6_T","Patient7_T","Patient8_T","Patient9_T"))]="Normal"
Data<-SetIdent(Data,value = as.vector(Data$Sampletype))
name1<-"Ident"
name2<-"HBO_Hepatocyte"
Prepare.for.GSEA(Data,name1,name2)


####Supervised heatmap
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

cell_Erythroid_like_Tumor<-colnames(HBO_ery_tumor)[HBO_ery_tumor$ID %in% c("Patient2_Erythroid-like-tumor","Patient5_H_Erythroid-like-tumor")]
dge$ID<-paste0(dge$patientID,"_",dge$ID_final)
dge$ID[colnames(dge) %in% cell_Erythroid_like_Tumor]="Erythroid_like_Tumor"
dge$ID[dge$ID_final=="BECs"]="BECs"
dge$ID[dge$ID_final=="Normal_Hepatocyte"]="Hepatocyte"
dge$ID[dge$ID_final=="Neuroendocrine"]="Neuroendocrine"
dge$ID[dge$ID_final=="DCN_high_Tumor"]="DCN_high_Tumor"

####Erythrocytes
WNT_markers<-read.table("~/Desktop/Seqwell_combined/combined_analysis/geneset/WNT_SIGNALING.txt",sep="\n")
Ery_markers<-c(early_ery,mid_ery,late_ery)
Ery_type<-c(rep("Early",length(early_ery)),rep("Mid",length(mid_ery)),rep("Late",length(late_ery)))
#for (i in 73:(length(WNT_markers$V1)-2)) {
Ery_markers<-read.table("~/Desktop/hepatoblastoma/Signature_geneset/Tumor_20genes.txt",sep="\n")
Ery_markers<-Ery_markers$V1
Ery_type<-c(rep("Tumor",length(Ery_markers)))
for (i in 1:(length(Ery_markers))) {
  Feature<-Ery_markers[i]
  #Feature<-"H19"
  V1<-as.vector(FetchData(object = dge,vars = Feature))
  FeaturePlot(dge,features=Feature)+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(V1)))
  ggsave(file=paste0(Ery_type[i],"_",Feature,"_featureplot.pdf"),width = 20,height = 20,units = "cm")
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1[,1])
  BOX_df<-data.frame(BOX_df)
  # ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin()+ #geom_jitter(shape=16,position = position_jitter(0.1))+
  #   stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
  #   stat_compare_means(method = "anova",label.x = 3,label.y = max(V1)+0.05)+
  #   ggtitle(Feature)+RotatedAxis()+NoLegend()
  VlnPlot(dge,features=Feature)+NoLegend()
  ggsave(file=paste0(Ery_type[i],"_",Feature,"_boxplot.pdf"),width = 20,height = 20,units = "cm")
}



####DCN_high_Tumor vs Fibroblast

dge<-SubsetData(HBO_final,cells = colnames(HBO_final)[HBO_final$patient %in% c("Patient3_T","Patient5_T")])
dge<-SubsetData(dge,ident.use = c("WNT5A_int_Fibroblast","WNT5A_high_Fibroblast","DCN_high_Tumor"))
dge<-SetIdent(dge,value = as.vector(dge@active.ident))
dge$ID<-as.vector(dge@active.ident)
dge$ID[dge$ID %in% c("WNT5A_high_Fibroblast","WNT5A_int_Fibroblast")]="Fibroblast"
dge<-SetIdent(dge,value = as.vector(dge$ID))
DimPlot(dge,group.by = "patientID")
dge$ID_BP<-paste0(dge$patientID,"_",dge$ID)
dge$ID_BP[dge$ID_BP %in% c("Patient5_Fibroblast","Patient3_Fibroblast")]="Fibroblast"
dge$ID_BP[dge$ID_BP =="Patient3_DCN_high_Tumor"]="Patient3_Tumor"
dge$ID_BP[dge$ID_BP =="Patient5_DCN_high_Tumor"]="Patient5_Tumor"
dge<-SetIdent(dge,value = as.vector(dge$ID_BP))
dge<-SetIdent(dge,value = as.vector(dge$ID))
pbmc.markers <- FindAllMarkers(dge, only.pos = TRUE, min.pct = 0.001, logfc.threshold = 0.001)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
write.table(pbmc.markers,"all_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")
top10_annotated <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(top10_annotated,"top10_markers_byID.txt",col.names = T,row.names = T,quote = F,sep="\t")
dge_temp<-SubsetData(dge,max.cells.per.ident = 100)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
DoHeatmap(dge_temp, features = top10_annotated$gene,raster = F,size = 3) + theme(axis.text.y = element_text(size = 5))+NoLegend()
ggsave(file="Heatmap_10_byID.pdf",width = 20,height = 20,units = "cm",limitsize = FALSE)
DotPlot(dge, cols = c("grey","red"), features =unique(top10_annotated$gene) ) + 
  coord_flip()+ RotatedAxis() +
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 6))
ggsave(file="Dotplot_10_byID.pdf",width = 15,height = 20,units = "cm",limitsize = FALSE)

unique(dge$ERGtype)

for (i in 1:10) {
  dge_temp<-SetIdent(dge,value = as.vector(dge$ID))
  dge_temp<-SubsetData(dge_temp,max.cells.per.ident = 46)
  all_markers <-FindAllMarkers(dge_temp,min.diff.pct=0.001,min.pct =  0.001,test.use = "wilcox",assay = "RNA")
  mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
  RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
  #RPL.genes <- grep(pattern = "^RPL", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
  remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
  all_markers<-all_markers[!(all_markers$gene %in% remove_gene),]
  write.table(all_markers,paste0(i,"_Tumor_Normal_downsampled_marker.txt"),sep="\t",row.names=F, col.names = T)
  all_markers<-all_markers[1:(nrow(all_markers)/2),]
  DE_select<-all_markers$p_val_adj<0.05
  length(which(DE_select))
  all_markers$threshold <- DE_select
  ggplot(all_markers) +
    geom_point(aes(x=avg_logFC, y=-log10(p_val_adj), colour=threshold)) +
    ggtitle("Tumor_Fibroblast_Volcano") +
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
    ggtitle("Tumor_Fibroblast_Volcano") +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) 
  ggsave(file=paste0(i,"_Tumor_Fibroblast_downsampled_Volcano.pdf"),width = 20,height = 10,units = "cm")
}





