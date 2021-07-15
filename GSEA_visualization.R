#####Visualize GSEA results
library(forcats)

####read in upregulated downregulated GSEA results in txt
dot_df_1<-read.table("./GSEA_ERGneg_Tumor_UP.txt",sep = "\t",header = T)
dot_df_2<-read.table("./GSEA_ERGneg_Tumor_DOWN.txt",sep = "\t",header = T)


##### take only the top20 and filter by the FDR < 0.1
dot_df<-rbind(dot_df_1[1:20,],dot_df_2[1:20,])
dot_df<-dot_df[dot_df$FDR.q.val<0.1,]
#### add in the type variable
dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"

#### below is to compute the gene ratio as it is now part of the LEADING.EDGE column by ","
info_column<-as.vector(dot_df$LEADING.EDGE)
DF_nes<-strsplit(info_column, ", ")

dot_df$tags=0
for (i in 1:length(DF_nes)) {
  dot_df$tags[i]<-as.numeric(as.vector(gsub("%", "", as.vector(gsub("tags=", "", DF_nes[[i]][1])))))
}
dot_df$GeneRatio<-dot_df$SIZE*dot_df$tags/100
dot_df$Counts<-dot_df$SIZE*dot_df$tags/100
dot_df$GeneRatio<-dot_df$tags/100

##### plot
p <- ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(NAME, GeneRatio))) + 
  geom_point(aes(size = Counts, color = FDR.q.val)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("C2CP GSEA")

p + facet_grid(.~type)
ggsave(file=paste0("./Common_upregulated_ERGneg.eps"),width = 30,height = 20,units = "cm")





