# assemble the seurat object

library(Seurat)
library(Matrix)


metadata <- read.table( "Helene/og_data/GSE178429_PBMCs_stim_scRNAseq_cellMeta.txt", header = TRUE, row.names = 1, sep = "\t")
counts <- readMM("Helene/og_data/PBMCs_stim_scRNAseq_counts.txt")
gene_names <- read.table("Helene/og_data/GSE178429_PBMCs_stim_scRNAseq_geneNames.txt", header = FALSE)

rownames(counts) <-  gene_names[, 1] 

seurat_obj <- CreateSeuratObject(counts, project = "LPS_PBMC", assay = "RNA",
  min.cells = 0, min.features = 0, names.field = 1,
  names.delim = "_", meta.data = metadata)


#cleaning
pbmc <- FindVariableFeatures(seurat_obj)
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 0.2)
pbmc <- RunUMAP(object = pbmc, dims = 1:10)
umap <- DimPlot(pbmc, reduction = "umap")
ggsave("Helene/plots/umap.png", plot = umap, width = 8, height = 6, dpi = 300)

# save
SaveSeuratRds(pbmc, "Helene/data/pmbc_seurat_object.rds")

