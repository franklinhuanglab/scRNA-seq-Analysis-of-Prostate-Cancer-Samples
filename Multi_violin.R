setwd("~/Desktop/Seqwell_combined/Revision/")
####Multi_violin
dge<-SetIdent(club_PCA,value = as.vector(club_PCA$target))
for (i in 1:length(unique(dge@active.ident))) {
selected_cells <- names(dge@active.ident[dge@active.ident == unique(dge@active.ident)[i]])
data <- FetchData(dge,
                  #vars = c("PIGR","LTF","MMP7","SCGB1A1","NKX3-1"),
                  #vars = c("KLK3","KLK2","ACPP","NKX3-1","AR","PIGR","LTF","MMP7","SCGB1A1"),
                  vars = c("KLK3","KLK2","ACPP","NKX3-1","AR","KRT5","KRT15","KRT17","DST","TP63"),
                  cells = selected_cells ,
                  slot = "data")

head(data)

long_data <- melt(data)
head(long_data)
ggplot(long_data,
       aes(x = variable, y = value)) +
  geom_boxplot() 
  #geom_jitter(size = 0.1)
ggsave(file=paste0("Multi_violinplot_",unique(dge@active.ident)[i],".pdf"),width = 30,height = 15,units = "cm")
}
DotPlot(club_PCA,features = c("KLK3","KLK2","ACPP","NKX3-1","AR","PIGR","LTF","MMP7","SCGB1A1"),cols = c("red","white"))+RotatedAxis()

ggsave(file=paste0("Multi_violinplot_Club.pdf"),width = 30,height = 15,units = "cm")
ggsave(file=paste0("Multi_violinplot_LE.pdf"),width = 30,height = 15,units = "cm")
ggsave(file=paste0("Multi_violinplot_BE.pdf"),width = 30,height = 15,units = "cm")

####Raw signature
dge<-dge_E
DimPlot(dge)
overlap_signature<-rownames(dge)
sig=list()
for (i in 1:length(unique(dge@active.ident))) {
  dge_temp<-SubsetData(dge,ident.use = unique(dge@active.ident)[i])
  dge_temp<-SetIdent(dge_temp,value=as.vector(dge_temp@active.ident))
  av.exp <- AverageExpression(dge_temp)$RNA
  colnames(av.exp)="Expression"
  av.exp$gene<-rownames(av.exp)
  av.exp <- av.exp[order(av.exp$Expression,decreasing = T),]
  mito.genes <- grep(pattern = "^MT-", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
  RPS.genes <- grep(pattern = "^RPS", x = rownames(dge_temp@assays$RNA@data), value = TRUE)
  remove_gene<-c(mito.genes,RPS.genes,overlap_signature)
  av.exp<-av.exp[!(av.exp$gene %in% remove_gene),]
  av.exp<-av.exp[av.exp$Expression>=10,]
  write.table(av.exp,paste0(unique(dge@active.ident)[i],"_RawSignature.txt"),sep="\t",col.names = T,row.names = T)
  sig[[i]]=av.exp
  #signature_temp<-av.exp$gene
  #overlap_signature<-intersect(signature_temp,overlap_signature)
}
names=unique(dge@active.ident)








