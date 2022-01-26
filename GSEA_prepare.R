#BiocManager::install("ArrayTools")
library(ArrayTools)
library(dplyr)

Prepare.for.GSEA <- function(Data,name1,name2){
  #Data<-dge_Fib
  #name1="Ident"
  #name2="Fibroblast"n 
  #Data<-dge_LE
  #name1<-"Ident"
  #name2<-c("LE_prostate")
  
  #Data: the seurat object we use
  #name1: the variable in the metadata we want to use for the phenotypes
  #Note: Here we should define which dataset and the idents we want to use (e.g Tumor vs. Benign)
  #name2: file name prefix
  
  #normal_LE<-setdiff(colnames(Data),cell_tumor)
  #Data<-SetIdent(object=Data, cells = cell_tumor, value = "Tumor")
  #Data<-SetIdent(object=Data, cells = normal_LE, value = "Benign")
  
  #Comment or uncomment the following line, depends on whether we want to downsample the dataset or not. 
  #Data.subset <- subset(x = Data,downsample=257)
  #Data<-RenameCells(object=Data,new.names = as.character(1:ncol(Data)))
  Data.subset<-Data
  
  #Process subsetted data
  DefaultAssay(Data.subset) <- "RNA"
  Data.subset <- DietSeurat(Data.subset, counts = TRUE, data = TRUE, scale.data = FALSE)
  Data.subset <- NormalizeData(Data.subset)
  Data.subset <- FindVariableFeatures(Data.subset,selection.method = 'vst',nfeatures = 10000)
  Data.subset <- ScaleData(Data.subset, verbose = FALSE,vars.to.regress = c('nCount_RNA'))
  #Data.subset <- RunPCA(Data.subset, npcs = 50, verbose = FALSE)
  Data.subset <- RunPCA(Data.subset, npcs = min((ncol(Data.subset)-1),50), verbose = FALSE)
  
  #Prepare metadata for GSEA
  Data.subset.meta <- Data.subset@meta.data
  #Create a column for our phenotype (i.e Ident)
  Data.subset.meta$Ident<-Data.subset@active.ident
  output.cls.name <- paste(name2,'phenotypes',sep='_')
  #name1 should be the column you want for the phenotype file of the GSEA analysis
  #.cls file generated 
  output.cls(Data.subset.meta,name1,filename = output.cls.name)
  
  #Prepare scaled data for GSEA
  Data.subset.scaled <- as.data.frame(as.matrix(Data.subset[['RNA']]@scale.data))
  #Convert format
  require(Biobase)
  data.ssGSEA <-new("ExpressionSet", exprs=as.matrix(Data.subset.scaled))
  output.gct.name <- paste(name2,'ssGSEA',sep='_')
  #gct file would be the expression file we use for the GSEA analysis
  output.gct(normal = data.ssGSEA,filename = output.gct.name)
}


###################################################