
library(SingleR)
library(DESeq2)
library(edgeR)
library(DEFormats)
library(Seurat)
library(cowplot)
library(dplyr)
library(jjb)
library(GEOquery)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggpubr)
library(ggsignif)
theme_set(theme_cowplot())
library(DoubletFinder)
hpca.se <- HumanPrimaryCellAtlasData()
bp.se <- BlueprintEncodeData()
hpca.se
library(scRNAseq)
dge_temp<-SubsetData(dge,ident.use = "15")
dge_temp<-club_PCA
hESCs <- as.matrix(HBO_final@assays$RNA@data)

# SingleR() expects log-counts, but the function will also happily take raw
# counts for the test dataset. The reference, however, must have log-values.

pred.combined_HPCA <- SingleR(test = hESCs, 
                              ref = hpca.se, 
                              labels = hpca.se$label.main)

pred.combined_BP <- SingleR(test = hESCs, 
                            ref = bp.se, 
                            labels = bp.se$label.main)

table(pred.combined_HPCA$labels)
plotScoreHeatmap(pred.combined_HPCA)

table(pred.combined_BP$labels)
plotScoreHeatmap(pred.combined_BP)

to.remove <- pruneScores(pred.combined_HPCA)
summary(to.remove)
plotScoreDistribution(pred.combined_HPCA, show = "delta.med", ncol = 3, show.nmads = 3)

all.markers <- metadata(pred.combined_BP)$de.genes
hESCs$labels <- pred.combined_BP$labels

plotHeatmap(hESCs, order_columns_by="labels",
            features=unique(unlist(all.markers$beta))) 

new.pruned <- pred.grun$labels
new.pruned[pruneScores(pred.grun, nmads=5)] <- NA
table(new.pruned, useNA="always")

dge$ID_BP<-pred.combined_BP$labels
DimPlot(dge,group.by = "ID_BP",label = T)

dge$ID_HPCA<-pred.combined_HPCA$labels
DimPlot(dge,group.by = "ID_HPCA",label = T)

library(DoubletFinder)
sweep.res.list_subs <- paramSweep_v3(dge, PCs = 1:100)
sweep.stats_subs <- summarizeSweep(sweep.res.list_subs, GT = F)
bcmvn_subs <- find.pK(sweep.stats_subs)


## Homotypic Doublet Proportion Estimate

# TODO: Adjust this number. By super Poisson distribution (for 10x and seqwell 10000 cells -> ~10% doublet rate)
doublet_formation_rate <- 0.08 ####https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php
doublet_formation_rate <- 0.06
doublet_formation_rate <- 0.1
homotypic.prop <- modelHomotypic(dge@active.ident)           ## subs@active.ident is the cell type annotations
nExp_poi <- round(doublet_formation_rate*length(Cells(dge)))  ## Total num doublets--Assuming x% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies

subs <- doubletFinder_v3(dge, PCs = 1:100, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE)
# reuse.paNN might need to change this, look in metadata
subs <- doubletFinder_v3(dge, PCs = 1:100, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE)


## Plot results

# TODO: change DF.classifications_0.25_0.09_143 to the correct column in subs meta data
subs@meta.data[,"DF_hi.lo"] <- subs@meta.data$DF.classifications_0.25_0.09_1192
subs@meta.data$DF_hi.lo[which(subs@meta.data$DF_hi.lo == "Doublet" & subs@meta.data$DF.classifications_0.25_0.09_1192 == "Singlet")] <- "Doublet_lo"
subs@meta.data$DF_hi.lo[which(subs@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
UMAPPlot(subs, group.by="DF_hi.lo")

