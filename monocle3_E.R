setwd("/Users/hsong/Desktop/Seqwell_combined/")
mkdir("./monocle")
setwd("./monocle/")
library(monocle)
library(DESeq2)
library(edgeR)
library(DEFormats)
library(Seurat)
library(cowplot)
library(dplyr)
library(icesTAF)
library(GEOquery)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ArrayTools)
library(dplyr)
# This script is designed to pipeline Seurat object output into Monocle3
# This script also contains all the necessary functionalities in Monocle3 as 
# This could also convert a 2D Seurat object and visualize/analyze it in 3D
# It will start from reading in a Seurat object with the count sparse matrix, UMAP coordinates, and cluster information
# This script is originally written for local machines but adaptations have also been included in annotations


### Require:: 'BiocManager', 'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'reticulate', 'htmlwidgets'

# This is a required python package
#reticulate::py_install("louvain")

# This is installing the actual monocle3
#devtools::install_github('cole-trapnell-lab/monocle3')


### Installing the packages

library(Seurat)
library(monocle)
library(htmlwidgets)
library(monocle3)
seurat_all<-readRDS("../E_analysis/E_refined.rds")
seurat_all<-dge_combine
samplelist<-unique(seurat_all$orig.ident)
samplelist_t=c("AUG_PB1","PR5249_T","PR5251_T","PR5254_T","PR5261_T","PR5269","MAY_PB1","MAY_PB2","PR5186","PR5196","PR5199")
samplelist_n<-c("PR5249_N","PR5251_N","PR5254_N","PR5261_N")
j=2
for (i in j:j) {
i=15
mkdir(paste0("./",samplelist[i]))
setwd(paste0("./",samplelist[i],"/"))

mkdir(paste0("./","Henry"))
setwd(paste0("./","Henry","/"))
#seurat<-SubsetData(seurat_all,cells = colnames(seurat_all)[seurat_all$orig.ident==samplelist[i]])
seurat<-SubsetData(seurat_all,cells = colnames(seurat_all)[seurat_all$orig.ident %in% samplelist_n],ident.remove = c("ERGpos_Tumor","ERGneg_Tumor"))
seurat<-Club_Henry
seurat$refinedID<-as.vector(seurat@active.ident)

#Extract data, phenotype data, and feature data from the SeuratObject
gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
gene_annotation <- as.data.frame(rownames(seurat), row.names = rownames(seurat))

colnames(gene_annotation) <- "gene_short_name"

# part two, cell information

cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
#colnames(cell_metadata) <- "barcode"
cell_metadata<-cbind(cell_metadata,as.vector(seurat$refinedID))
colnames(cell_metadata)<-c("barcode","refinedID")
# part three, counts sparse matrix

New_matrix <- seurat@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix


### Construct the basic cds object

HSMM <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

#View data
#pData(HSMM)
#fData(HSMM)

HSMM <- estimate_size_factors(HSMM)

#print(dim(exprs(HSMM)))

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM<-preprocess_cds(HSMM,method = "PCA")
HSMM <- reduce_dimension(HSMM)

# First decide what you want to color your cells by
#print(head(pData(HSMM)))

## order cells change colors and theta to match your plot
HSMM<-cluster_cells(HSMM,reduction_method = "UMAP")
HSMM<-learn_graph(HSMM)
fig1<-plot_cells(HSMM,reduction_method = "UMAP",color_cells_by = "refinedID",label_cell_groups = T,label_roots = T,label_groups_by_cluster = T,cell_size = 1)
pdf(paste0(cellname,"_cellplot.pdf"), width = 10, height = 10)
fig1
dev.off()

fig2<-plot_cells(HSMM,1,2,reduction_method = "UMAP",color_cells_by = "refinedID",label_cell_groups = T,label_roots = T,cell_size = 1)+facet_wrap(~refinedID)
pdf(paste0(cellname,"_cellplot_separate.pdf"), width = 10, height = 10)
fig2
dev.off()

plot_cells(HSMM,reduction_method = "UMAP",color_cells_by = "refinedID",label_cell_groups = T,label_roots = T,label_groups_by_cluster = T,cell_size = 1)


HSMM <- order_cells(HSMM)
fig3<-plot_cells(HSMM,reduction_method = "UMAP",show_trajectory_graph = T,color_cells_by = "refinedID",label_cell_groups = T,label_roots = T,label_groups_by_cluster = T,cell_size = 1)
pdf(paste0(cellname,"_trajectory.pdf"), width = 10, height = 10)
fig3
dev.off()

fig4<-plot_cells(HSMM,reduction_method = "UMAP",show_trajectory_graph = T,color_cells_by = "pseudotime",label_cell_groups = T,label_roots = T,label_groups_by_cluster = T,cell_size = 1)
pdf(paste0(samplelist[i],"_pseudotime.pdf"), width = 10, height = 10)
fig4
dev.off()

ciliated_cds_pr_test_res <- graph_test(HSMM, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))


gene_module_df <- find_gene_modules(HSMM[pr_deg_ids,], resolution=c(0.0001,10^seq(-3,-1)))
cell_group_df <- tibble::tibble(cell=row.names(colData(HSMM)), 
                                cell_group=colData(HSMM)$refinedID)
agg_mat <- aggregate_gene_expression(HSMM, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
fig5<-pheatmap::pheatmap(agg_mat,scale="column", clustering_method="ward.D2")
pdf(paste0(samplelist[i],"_DE.pdf"), width = 10, height = 10)
fig5
dev.off()
write.table(gene_module_df[order(gene_module_df$module,decreasing = F),],"DE_module.txt",col.names = T,row.names = T,quote = F,sep = "\t")
 plot_cells(HSMM,
            genes=gene_module_df %>% filter(module %in% c(14,2,6)),
            label_cell_groups=FALSE,
            show_trajectory_graph=FALSE)

genes <- c("KLK3","AR","FOS","JUN","SCGB3A1","LCN2","PSCA","KRT13","HSPA8","SERPINH1")
lineage_HSMM <- HSMM[rowData(HSMM)$gene_short_name %in% genes,]
fig6<-plot_genes_in_pseudotime(lineage_HSMM)
pdf(paste0(cellname,"_marker.pdf"), width = 5, height = 20)
fig6
dev.off()
 setwd("../")
}



