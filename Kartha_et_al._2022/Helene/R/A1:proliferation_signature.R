# Are there cells with a proliferation signature?
# marker genes for monocyte ID


library(Seurat)
library(ggplot2)

pbmc <- readRDS("Helene/data/pmbc_seurat_object.rds")

# predefined cell cycle genes
cc.genes <- Seurat::cc.genes
cc.genes$s.genes

# check proliferation/cell cycle genes
pbmc <- CellCycleScoring(
  object = pbmc,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE  
)
# add proliferation signature
proliferation_genes <- c("MKI67", "PCNA", "TOP2A")  
proliferation_score <- colMeans(GetAssayData(pbmc, slot = "data")[proliferation_genes, , drop = FALSE])
pbmc <- AddMetaData(pbmc, metadata = proliferation_score, col.name = "Proliferation1")
# show 
FeaturePlot(pbmc, features = "Proliferation1")


proliferative_cells <- subset(pbmc, subset = Proliferation1 > 0.5) 

# yes there are cells with a prolifaration signature and it is distributed across the clusters. 


