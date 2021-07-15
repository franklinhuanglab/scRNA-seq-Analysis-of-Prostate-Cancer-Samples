s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps

dge <- NormalizeData(dge)
dge <- FindVariableFeatures(dge, selection.method = "vst")
dge <- ScaleData(dge, features = rownames(dge))
dge <- RunPCA(dge, features = VariableFeatures(dge), ndims.print = 6:10, nfeatures.print = 10)
DimHeatmap(dge, dims = c(8, 10))
dge <- CellCycleScoring(dge, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(dge[[]])
RidgePlot(dge, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
RidgePlot(dge_E, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2,group.by = "ID")
ggsave(file="Differentiation_markers_ridgeplot.eps",width = 20,height = 20,units = "cm")
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
dge <- RunPCA(dge, features = c(s.genes, g2m.genes))
DimPlot(dge)
ggsave(file="G2M_UMAP_Organoid.eps",width = 20,height = 20,units = "cm")
tab<-table(dge@active.ident,dge$old.ident)
write.table(tab,"G2M_Organoid.txt",sep = "\t",col.names = T,row.names = T )


dge <- ScaleData(dge, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(dge))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
dge <- RunPCA(dge, features = VariableFeatures(dge), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
dge <- RunPCA(dge, features = c(s.genes, g2m.genes))
DimPlot(dge)
DimPlot(dge,group.by = "old.ident")
ggsave(file="G2M_UMAP_orgsanoid_2.eps",width = 20,height = 20,units = "cm")
tab<-table(dge@active.ident,dge$old.ident)
write.table(tab,"G2M_LE_2.txt",sep = "\t",col.names = T,row.names = T )


