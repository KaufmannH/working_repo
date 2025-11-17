
# env basic packages: regnoise_env.yml

library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
library(ggplot2)
library(pryr) # RAM monitoring


run_pipeline_bpcells <- function(obj_path, object_name, output_dir = "pipeline") {

 
  obj <- readRDS(obj_path)
  obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay5")

  # re save in new format
  matrix_dir = file.path(output_dir, "data", paste0(object_name, "_counts"))
  if (dir.exists(matrix_dir)) {unlink(matrix_dir, recursive = TRUE)} # will throw an error otherwise

  write_matrix_dir(mat = obj[["RNA"]]$counts, dir = matrix_dir)
  counts.mat <- open_matrix_dir(dir = file.path(output_dir, "data", paste0(object_name, "_counts")))
  obj[["RNA"]]$counts <- counts.mat
  rm(counts.mat)
  gc()
  print("new matrix is done")


  # Normalisation
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)

  # seurat pipeline
  start <- Sys.time()
  ram_start <- mem_used()

  obj <- FindVariableFeatures(obj, verbose = FALSE)
  print("var_feat done")
  obj <- ScaleData(obj, verbose = FALSE)
  print("scaling done")
  obj <- RunPCA(obj, verbose = FALSE)
  print("PCA done")
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
  print("find neighbors done")
  obj <- FindClusters(obj, algorithm = 1, verbose = FALSE)
  print("find clusters done")
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, return.model = TRUE, verbose = FALSE)
  print("umap done")
  end <- Sys.time()
  ram_end <- mem_used()
  duration <- end - start

  # umap
  umap_plot <- DimPlot(obj, group.by = "sample", reduction = "umap") +
    theme(
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 8)
    ) +
    theme_classic()
    
  ggsave(file.path(output_dir, "plots", paste0("umap_", object_name, ".png")),
         umap_plot, width = 13, height = 8, dpi = 300)
  print("plot done")

# time info
  log_file <- file.path(output_dir, "logs", paste0("pipeline_log_", object_name, ".txt"))
  writeLines(
  c(
    paste("Object:", object_name),
    paste("Start:", start),
    paste("End:", end),
    paste("Duration (seconds):", as.numeric(duration, units = "secs")),
    paste("RAM at start:", round(as.numeric(ram_start) / 1024^3, 3), "GB"),
    paste("RAM at end:  ", round(as.numeric(ram_end) / 1024^3, 3), "GB"),
    paste("RAM diff:    ", round(as.numeric(ram_end - ram_start) / 1024^3, 3), "GB")),
    con = log_file)
  return(obj)
}



run_pipeline_seurat <- function(obj_path, object_name, output_dir = "pipeline") {
  
  obj <- readRDS(obj_path)
  obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay5")

  # Normalisation
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)

  # seurat pipeline
  start <- Sys.time()
  ram_start <- mem_used()

  obj <- FindVariableFeatures(obj, verbose = FALSE)
  print("var_feat done")
  obj <- ScaleData(obj, verbose = FALSE)
  print("scaling done")
  obj <- RunPCA(obj, verbose = FALSE)
  print("PCA done")
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
  print("find neighbors done")
  obj <- FindClusters(obj, algorithm = 1, verbose = FALSE)
  print("find clusters done")
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, return.model = TRUE, verbose = FALSE)
  print("umap done")
  end <- Sys.time()
  ram_end <- mem_used()
  duration <- end - start

  # umap
  umap_plot <- DimPlot(obj, group.by = "sample", reduction = "umap") +
    theme(
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 8)
    ) +
    theme_classic()
    
  ggsave(file.path(output_dir, "plots", paste0("umap_", object_name, ".png")),
         umap_plot, width = 13, height = 8, dpi = 300)
  print("plot done")

# time info
  log_file <- file.path(output_dir, "logs", paste0("pipeline_log_", object_name, ".txt"))
  writeLines(
  c(
    paste("Object:", object_name),
    paste("Start:", start),
    paste("End:", end),
    paste("Duration (seconds):", as.numeric(duration, units = "secs")),
    paste("RAM at start:", round(as.numeric(ram_start) / 1024^3, 3), "GB"),
    paste("RAM at end:  ", round(as.numeric(ram_end) / 1024^3, 3), "GB"),
    paste("RAM diff:    ", round(as.numeric(ram_end - ram_start) / 1024^3, 3), "GB")),
    con = log_file)
  return(obj)
}



execute <- function() {
  obj <- run_pipeline_bpcells("pipeline/data/myeloid_object.RDS", "obj")
  obj_seurat <- run_pipeline_seurat("pipeline/data/myeloid_object.RDS", "obj_seurat")
  
}
execute()


