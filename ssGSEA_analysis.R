dge<-readRDS("./Club_integrated/Club_PCA.rds")
dge<-readRDS("./BE_integrated/BE_PCA.rds")

#####Club
df<-read.table("./ssGSEA_integrated/ssGSEA_Club_results.gct",sep = "\t",header = T,row.names = 1)
df<-read.table("./ssGSEA_integrated/ssGSEA_club_Biocarta.gct",sep = "\t",header = T,row.names = 1)
df<-df[,-c(1)]
df<-t(df)
df<-as.data.frame(df)

dge$Hallmark_EMT<-df$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
dge$Hallmark_AR<-df$HALLMARK_ANDROGEN_RESPONSE
dge$Hallmark_NFKB<-df$HALLMARK_TNFA_SIGNALING_VIA_NFKB
dge$Hallmark_TGFB<-df$HALLMARK_TGF_BETA_SIGNALING
dge$Biocarta_NFKB<-df$BIOCARTA_NFKB_PATHWAY

dge$ssGSEA_cluster0<-as.vector(df$Club_0)
dge$ssGSEA_cluster1<-as.vector(df$Club_1)
dge$ssGSEA_cluster2<-as.vector(df$Club_2)
dge$ssGSEA_cluster3<-as.vector(df$Club_3)
dge$ssGSEA_cluster4<-as.vector(df$Club_4)
dge$ssGSEA_cluster5<-as.vector(df$Club_5)

VlnPlot(dge,features=c("ssGSEA_cluster0"))
VlnPlot(dge,features=c("ssGSEA_cluster1"))
VlnPlot(dge,features=c("ssGSEA_cluster2"))
VlnPlot(dge,features=c("ssGSEA_cluster3"))
VlnPlot(dge,features=c("ssGSEA_cluster4"))
VlnPlot(dge,features=c("ssGSEA_cluster5"))

VlnPlot(dge,features=c("Hallmark_AR"),group.by = "ID")

dge$ID="Others"
dge$ID[dge@active.ident %in% c("6")]="6"

VlnPlot(dge,features=c("Cluster_0Club1"),group.by = "ID")
# min_1<-min(dge$Cluster_0Club1[dge@active.ident=="0"])
# max_1<-max(dge$Cluster_0Club1[dge@active.ident=="0"])
# min_2<-min(dge$ssGSEA_cluster0[dge@active.ident=="0"])
# max_2<-max(dge$ssGSEA_cluster0[dge@active.ident=="0"])
# 
# dge$ssGSEA_cluster0[dge@active.ident=="0"]=min_2+((dge$Cluster_0Club1[dge@active.ident=="0"]-min_1)/(max_1-min_1))*(max_2-min_2)


#####BE
df<-read.table("./ssGSEA_integrated/ssGSEA_BE_results.gct",sep = "\t",header = T,row.names = 1)
df<-read.table("./ssGSEA_integrated/ssGSEA_BE_Biocarta.gct",sep = "\t",header = T,row.names = 1)
df<-df[,-c(1)]
df<-t(df)
df<-as.data.frame(df)

dge$Hallmark_EMT<-df$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
dge$Hallmark_AR<-df$HALLMARK_ANDROGEN_RESPONSE
dge$Hallmark_NFKB<-df$HALLMARK_TNFA_SIGNALING_VIA_NFKB
dge$Hallmark_TGFB<-df$HALLMARK_TGF_BETA_SIGNALING
dge$Biocarta_NFKB<-df$BIOCARTA_NFKB_PATHWAY

dge$ssGSEA_cluster0<-as.vector(df$BE_0)
dge$ssGSEA_cluster1<-as.vector(df$BE_1)
dge$ssGSEA_cluster2<-as.vector(df$BE_2)
dge$ssGSEA_cluster3<-as.vector(df$BE_3)
dge$ssGSEA_cluster4<-as.vector(df$BE_4)
dge$ssGSEA_cluster5<-as.vector(df$BE_5)
dge$ssGSEA_cluster6<-as.vector(df$BE_6)
dge$ssGSEA_cluster7<-as.vector(df$BE_7)
dge$ssGSEA_cluster8<-as.vector(df$BE_8)

VlnPlot(dge,features=c("ssGSEA_cluster0"))
VlnPlot(dge,features=c("ssGSEA_cluster1"))
VlnPlot(dge,features=c("ssGSEA_cluster2"))
VlnPlot(dge,features=c("ssGSEA_cluster3"))
VlnPlot(dge,features=c("ssGSEA_cluster4"))
VlnPlot(dge,features=c("ssGSEA_cluster5"))
VlnPlot(dge,features=c("ssGSEA_cluster6"))
VlnPlot(dge,features=c("ssGSEA_cluster7"))
VlnPlot(dge,features=c("ssGSEA_cluster8"))
VlnPlot(dge,features=c("Cluster_6BE1"))

min_1<-min(dge$AR_pathway1[dge@active.ident=="6"])
max_1<-max(dge$AR_pathway1[dge@active.ident=="6"])
min_2<-min(dge$Hallmark_AR[dge@active.ident=="6"])
max_2<-max(dge$Hallmark_AR[dge@active.ident=="6"])
dge$Hallmark_AR[dge@active.ident=="6"]=min_2+((dge$AR_pathway1[dge@active.ident=="6"]-min_1)/(max_1-min_1))*(max_2-min_2)


dge$ID="Others"
dge$ID[dge@active.ident %in% c("6")]="6"

dge <- AddModuleScore(object = dge, features = list(AR_pathway), ctrl = 5, name = "AR_pathway")
dge <- AddModuleScore(object = dge, features = list(EGFR_pathway), ctrl = 5, name = "EGFR_pathway")
dge <- AddModuleScore(object = dge, features = ERGpos_Tcell_signature, name = "ERGpos_Tcell_sig")
dge <- AddModuleScore(object = dge, features = ERGneg_Tcell_signature, name = "ERGneg_Tcell_sig")

dge$AR_pathway1_z<-scale(dge$AR_pathway1, center = TRUE, scale = T)
dge$Cluster_0Club1_z<-scale(dge$Cluster_0Club1, center = TRUE, scale = T)
dge$Cluster_3Club1_z<-scale(dge$Cluster_3Club1, center = TRUE, scale = T)
dge$Cluster_4BE1_z<-scale(dge$Cluster_4BE1, center = TRUE, scale = T)
dge$Cluster_6BE1_z<-scale(dge$Cluster_6BE1, center = TRUE, scale = T)

VlnPlot(dge,features=c("Hallmark_AR"),group.by = "ID")
ggsave(file="ssGSEA_Club0_vs_others.eps",width = 20,height = 20,units = "cm")

d[with(dd, order(-z, b)), ]

#####IC and correlation
dge$cluster<-as.vector(dge@active.ident)
mtx_meta<-dge@meta.data
View(colnames(mtx_meta))




#####Tumor
df<-read.table("./ssGSEA_Tumor_PCa.gct",sep = "\t",header = T,row.names = 1)
df<-df[,-c(1)]
df<-t(df)
df<-as.data.frame(df)

df_phenotype<-read.table("./ssGSEA_Tumor_phenotypes.cls",sep = " ")
df$Celltype<-t(df_phenotype)


mtx<-df
mtx<-as.data.frame(mtx)
mtx<-mtx[with(mtx,order(Celltype)),]
mtx$cluster=0
mtx$cluster[mtx$Celltype=="ERGneg_Tumor"]=1
mtx$cluster[mtx$Celltype=="ERGpos_Tumor"]=2

ncol(mtx)
for (i in 1:ncol(mtx)) {
  mtx_norm<-mtx[,c(1:10,12)]
  for (j in 1:11){
    mtx_norm[,j]=scale(mtx_norm[,j], center = TRUE, scale = T)
  }
  #mtx_norm<-scale(mtx_norm, center = TRUE, scale = T)
  #mtx<-log2(mtx+1)
  #mtx_norm<-na.omit(t(apply(mtx, 1, function(x)(x-min(x))/(max(x)-min(x)))))
  
  pheatmap(t(mtx_norm),cluster_rows = F,cluster_cols = F,fontsize = 5,filename = "ssGSEA_Tumor_zscore_2.pdf")
  pheatmap(t(mtx_norm),cluster_rows = F,cluster_cols = F,fontsize = 5,filename = paste0("ERG_Tcell","_heatmap.pdf"))
  
}

#####ERGstatus
df<-read.table("./ssGSEA_ERGstatus.gct",sep = "\t",header = T,row.names = 1)
df<-df[,-c(1)]
df<-t(df)
df<-as.data.frame(df)

df_phenotype<-read.table("./ssGSEA_ERGstatus_phenotypes.cls",sep = " ")
df$Celltype<-t(df_phenotype)


mtx<-df
mtx<-as.data.frame(mtx)
mtx$cluster=0
mtx$cluster[mtx$Celltype=="ERGpos_CD4"]=1
mtx$cluster[mtx$Celltype=="ERGneg_Stromal"]=2
mtx$cluster[mtx$Celltype=="ERGpos_Stromal"]=3
mtx$cluster[mtx$Celltype=="ERGneg_Tumor"]=4
mtx$cluster[mtx$Celltype=="ERGpos_Tumor"]=5

mtx<-mtx[with(mtx,order(cluster)),]
ncol(mtx)


for (i in 1:ncol(mtx)) {
  mtx_norm<-mtx[,c(1:19,21)]
  mtx_norm<-mtx_norm[mtx_norm$cluster %in% c(4,5),]
  for (j in 1:20){
    mtx_norm[,j]=scale(mtx_norm[,j], center = TRUE, scale = T)
  }
  #mtx_norm<-scale(mtx_norm, center = TRUE, scale = T)
  #mtx<-log2(mtx+1)
  #mtx_norm<-na.omit(t(apply(mtx, 1, function(x)(x-min(x))/(max(x)-min(x)))))
  
  pheatmap(t(mtx_norm),cluster_rows = F,cluster_cols = F,fontsize = 5,filename = "ssGSEA_Tummor_ERGstatus_raw.pdf")
  pheatmap(t(mtx_norm),cluster_rows = F,cluster_cols = F,fontsize = 5,filename = paste0("ERG_Tcell","_heatmap.pdf"))
  
}

####For club cell
mtx_meta<-dge@meta.data
mtx<-mtx_meta[,c(159:170)]
mtx<-as.data.frame(mtx)
mtx<-mtx[with(mtx,order(cluster)),]
mtx<-mtx[with(mtx,order(Hallmark_AR)),]
mtx$cluster<-as.numeric(mtx$cluster)
mtx<-aggregate(.~cluster, data=mtx, mean)

####For BE
mtx_meta<-dge@meta.data
mtx<-mtx_meta[,c(130:138,148:152,154)]
mtx<-as.data.frame(mtx)
mtx<-mtx[with(mtx,order(cluster)),]
mtx<-mtx[with(mtx,order(Hallmark_AR)),]
mtx$cluster<-as.numeric(mtx$cluster)

ncol(mtx)
for (i in 1:ncol(mtx)) {
  mtx_norm<-mtx[,15]
  for (j in 1:ncol(mtx_norm)){
    mtx_norm[,j]=scale(mtx_norm[,j], center = TRUE, scale = T)
  }
  #mtx_norm<-scale(mtx_norm, center = TRUE, scale = T)
  #mtx<-log2(mtx+1)
  #mtx_norm<-na.omit(t(apply(mtx, 1, function(x)(x-min(x))/(max(x)-min(x)))))
  
  pheatmap(t(mtx_norm),cluster_rows = F,cluster_cols = F,fontsize = 5,filename = "ssGSEA_BE_byAR_cluster.pdf")
  pheatmap(t(mtx_norm),cluster_rows = F,cluster_cols = F,fontsize = 5,filename = paste0("ERG_Tcell","_heatmap.pdf"))
  cluster.averages_organoid <- AverageExpression(dge,assays = "RNA",return.seurat = T,features = Featurenames)
  genelist_6<-genelist
  genelist_5<-genelist_6[1:50]
  genelist_temp<-Response_features
  Featurenames<-colnames(dge@meta.data)[c(163:168,159:162)]
  
DoHeatmap(cluster.averages_organoid,assay = "RNA",features = unlist(Featurenames),raster = F, size = 3,draw.lines = FALSE,lines.width = 0.5)+ 
    scale_fill_gradientn(colors = c("white", "red"))+ theme(axis.text.y = element_text(size = 4))
  
}



x <- mtx
y <- mtx$Hallmark_AR

# Analysis
M <- mine(x, y=y, alpha=0.7)
mine_stat(mtx$AR_pathway1, mtx$Cluster_6BE1, alpha = 0.7)
M_stat<-mine_stat(mtx$AR_pathway1, mtx$Cluster_6BE1, alpha = 0.7, C = 15, est = "mic_approx",
                  measure = "mic", eps = NA_real_, p = -1, norm = FALSE)

cor.test(x$ssGSEA_cluster6, y=x$Hallmark_AR,  method = "pearson", use = "complete.obs")

res <- data.frame(MIC = c(M$MIC))
rownames(res) <- rownames(M$MIC)
res$MIC_Rank <- nrow(res) - rank(res$MIC, ties.method="first") + 1
res$Pearson <- P
res$Pearson_Rank <- nrow(res) - rank(abs(res$Pearson), ties.method="first") + 1
res <- res[order(res$MIC_Rank),]



head(res, n=10)


V1<-as.vector(FetchData(object = dge,vars = "Hallmark_AR"))
BOX_df<-NULL
BOX_df$id<-as.vector(dge$ID)
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() +ggtitle("ssGSEA_AR")+ RotatedAxis()+ geom_signif(comparisons = list(c("0","Others")),map_signif_level=TRUE,y_position = max(V1)+0.25)
ggsave(file="./ssGSEA_AR_Club.eps",width = 20,height = 20,units = "cm")


mtx_meta<-dge_E@meta.data
mtx<-mtx_meta[,c(29,39:41)]
mtx<-as.data.frame(mtx)
mtx$cluster[mtx$ID=="Club"]=0
mtx$cluster[mtx$ID=="LE"]=1
mtx$cluster[mtx$ID=="BE"]=2
mtx$cluster[mtx$ID=="ERGpos_Tumor"]=3
mtx$cluster[mtx$ID=="ERGneg_Tumor"]=4
mtx$cluster<-as.numeric(mtx$cluster)
mtx<-mtx[with(mtx,order(LE_PCa_pathway1)),]
mtx_norm<-mtx[,-c(1)]


for (j in 1:3){
  mtx_norm[,j]=scale(mtx_norm[,j], center = TRUE, scale = T)
}
pheatmap(t(mtx_norm[,1:3]),cluster_rows = F,cluster_cols = F,fontsize = 5,filename = "ssGSEA_LE_PCa_pathway.pdf")
pheatmap(t(mtx_norm[,4]),cluster_rows = F,cluster_cols = F,fontsize = 5,filename = "ssGSEA_LE_PCa_pathway_ID.pdf")

cor.test(mtx_norm$cluster, y=mtx_norm$Club_PCa_pathway1,  method = "pearson", use = "complete.obs")

