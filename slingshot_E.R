# By HS 01072020
library(DESeq2)
library(edgeR)
library(DEFormats)
library(Seurat)
library(cowplot)
library(dplyr)
library(GEOquery)
library(ggplot2)
library(ggsignif)
library(data.table)
library(ggrepel)
library(cluster)
library(reshape2)
library(qusage)
library(ggplot2)
library(ggsignif)
library(gplots)
#BiocManager::install("RTCGAToolbox")
#BiocManager::install("Homo.sapiens")
library(RTCGAToolbox)
library(Homo.sapiens)
#install.packages("ggrepel") #don't use the developmental version cause it's not compatible with ggplot2.
theme_set(theme_cowplot())
library(slingshot)
library(SingleCellExperiment)
library(destiny)
library(RColorBrewer)
library(mclust, quietly = TRUE)
library(gam)
setwd("/Users/hsong/Desktop/Seqwell_combined/")
mkdir("./slingshot")
setwd("./slingshot/")

dge<-dge_E
cellname="E_PCA"
nametype=as.vector(unique(dge@active.ident))
sim<-SingleCellExperiment(assays=List(counts=dge@assays$RNA@counts),colData=as.vector(dge@active.ident))
geneFilter <- apply(assays(sim)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sim <- sim[geneFilter, ]

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

assays(sim)$norm <- FQnorm(assays(sim)$counts)
pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)


dm <- DiffusionMap(t(log1p(assays(sim)$norm)))
rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)


cl1 <- Mclust(rd1)$classification
cl1<-as.vector(sim@colData@listData$X)
cl1[cl1=="BE"]=1
cl1[cl1=="Club"]=2
cl1[cl1=="LE"]=3
cl1[cl1=="ERGpos_Tumor"]=4
cl1[cl1=="ERGneg_Tumor"]=5
cl1<-as.numeric(cl1)
colData(sim)$GMM <- cl1

pdf(paste0(cellname,"_PCA.pdf"))
plot(rd2, col = brewer.pal(9,"Set1")[cl1], pch=16, asp=1)
legend("bottomleft", legend=nametype[unique(sim$GMM)],
       col=unique(brewer.pal(9,"Set1")[cl1]), pch=16,box.lty=0,bg="transparent")
dev.off()

cl2 <- kmeans(rd1, centers = 5)$cluster
colData(sim)$kmeans <- cl2
pdf(paste0(cellname,"_PCA_kmeans.pdf"))
plot(rd2, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)
legend("bottomleft", legend=nametype[unique(sim$GMM)],
       col=unique(brewer.pal(9,"Set1")[cl2]), pch=16,box.lty=0,bg="transparent")
dev.off()

sim <- slingshot(sim, clusterLabels = 'X', reducedDim = 'PCA')
summary(sim$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]

pdf(paste0(cellname,"_pseudotime_1.pdf"))
plot(reducedDims(sim)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col='black')
legend("bottomleft", legend=unique(sim$X),
       col=plotcol, pch=16,box.lty=1,bg="transparent")
dev.off()



pdf(paste0(cellname,"_pseudotime_2.pdf"))
plot(reducedDims(sim)$PCA, col = brewer.pal(9,'Set1')[sim$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')
legend("bottomleft", legend=unique(sim$X),
       col=unique(brewer.pal(9,"Set1")[cl1]), pch=16,box.lty=0,bg="transparent")
dev.off()

lin1 <- getLineages(rd2, cl1, start.clus = '1')
pdf(paste0(cellname,"_pseudotime_2.pdf"))
plot(rd2, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin1, lwd = 3, col = 'black')
legend("bottomleft", legend=unique(sim$X),
       col=unique(brewer.pal(9,"Set1")[cl1]), pch=16,box.lty=0,bg="transparent")
dev.off()


crv1 <- getCurves(lin1)
crv1
pdf(paste0(cellname,"_pseudotime_3.pdf"))
plot(rd2, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(crv1, lwd = 3, col = 'black')
legend("bottomleft", legend=unique(sim$X),
       col=unique(brewer.pal(9,"Set1")[cl1]), pch=16,box.lty=0,bg="transparent")
dev.off()


sim5 <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA',
                  approx_points = 5)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sim5$slingPseudotime_1, breaks=100)]

plot(reducedDims(sim5)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim5), lwd=2, col='black')
legend("bottomleft", legend=unique(sim5$X),
       col=unique(brewer.pal(9,"Set1")[cl1]), pch=16,box.lty=0,bg="transparent")

# t <- sim$slingPseudotime_1
# 
# # for time, only look at the 100 most variable genes
# Y <- log1p(assays(sim)$norm)
# var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
# Y <- Y[var100,]
# 
# # fit a GAM with a loess term for pseudotime
# gam.pval <- apply(Y,1,function(z){
#   d <- data.frame(z=z, t=t)
#   suppressWarnings({
#     tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
#   })
#   p <- summary(tmp)[3][[1]][2,3]
#   p
# })
# 
# topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
# mito.genes <- grep(pattern = "^MT-", x = topgenes, value = TRUE)
# RPS.genes <- grep(pattern = "^RPS", x = topgenes, value = TRUE)
# RPL.genes <- grep(pattern = "^RPL", x = topgenes, value = TRUE)
# remove_gene<-c(mito.genes,RPS.genes)#,RPL.genes)
# topgenes<-topgenes[!( topgenes %in% remove_gene)]
# heatdata <- assays(sim)$norm[topgenes, order(t, na.last = NA)]
# heatclus <- sim$GMM[order(t, na.last = NA)]
# 
# pdf(paste0(cellname,"_marker_transition.pdf"))
# heatmap.2(log1p(heatdata), ColSideColors = brewer.pal(9,"Set1")[heatclus],Colv=F,Rowv = T,trace="none",scale="row",col=colorRampPalette(c("blue","white","red"))(256))
# legend(xpd =T ,x=0.06,y=21,  legend=nametype[unique(heatclus)],
#        col=unique(brewer.pal(9,"Set1")[heatclus]), pch=15,box.lty=0,bg="transparent")
# dev.off()

setwd("../")










               