# This script was used to convert the PanSci h.5ad objects to .rds


library(Seurat)
library(zellkonverter)


path <- "data/subsets"
files <- list.files(path, pattern = "\\.h5ad$", full.names = TRUE)

for (f in files) {
  message("Processing: ", f)
  sce <- readH5AD(f)
  seurat_obj <- as.Seurat(sce, counts = "X", data = "X")
  out_file <- sub("\\.h5ad$", ".rds", f)
  saveRDS(seurat_obj, out_file)
  message("Saved: ", out_file)
}