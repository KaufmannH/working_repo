
# env basic packages: regnoise_env.yml

library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
library(ggplot2)



# example analysis


obj_seurat <- readRDS("pipeline/data/so.rds")
obj_seurat[["RNA"]] <- as(object = obj_seurat[["RNA"]], Class = "Assay5")

obj <- readRDS("pipeline/data/so.rds")
obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay5")
dim(obj[['RNA']])



# # Write the counts layer to a directory
 write_matrix_dir(mat = obj[["RNA"]]$counts, dir = 'pipeline/data/so_temp.RDS')
 counts.mat <- open_matrix_dir(dir = 'pipeline/data/so_temp.RDS')
 obj[["RNA"]]$counts <- counts.mat

# plots
plot <- VlnPlot(obj, features = c("Sox10", "Slc17a7", "Aif1"), ncol = 3, layer = "counts", alpha = 0.1)
ggsave("pipeline/plots/vln_before_norm.png", plot)
#normalize
obj <- NormalizeData(obj, normalization.method = "LogNormalize")
VlnPlot(obj, features = c("Sox10", "Slc17a7", "Aif1"), ncol = 3, layer = "counts", alpha = 0.1)

plot <- VlnPlot(obj, features = c("Sox10", "Slc17a7", "Aif1"), ncol = 3, layer = "data", alpha = 0.1)
ggsave("pipeline/plots/vln_after_norm.png", plot)





object <- obj # or obj_seruat

start <- Sys.time()
object <- FindVariableFeatures(object, verbose = F)
object <- ScaleData(object, verbose = F)
object <- RunPCA(object, verbose = F)
object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
object <- FindClusters(object, algorithm = 1)
object <- RunUMAP(object, reduction = "pca", dims = 1:30, return.model = T, verbose = F)
end <- Sys.time()

umap <- DimPlot(object, group.by = "sample", reduction = "umap") +
 theme(
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8)) +
  theme_classic()
ggsave("pipeline/plots/umap.png", umap,  width = 13, height = 8)

head(object@meta.data)



# how long did it take? 
duration <- end - start






