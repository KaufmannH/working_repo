# This script was used to rename the metadata cols in the converted .rds PanSci object to the names used in the pipeine scripts


#library('tidyverse')
library('Seurat')
library('dplyr')


# comparison: what is different that i need? 
c <- readRDS('/home/hkaufm49/analyses/RegulatedNoiseHel/data/rawdata/Spleen_rawdata.rds')
head(c@meta.data)
colnames(c@meta.data)

# example data
t <- readRDS('data/subsets/Female.03_months.rds')
dim(t)
head(t@meta.data)
colnames(t@meta.data)



# do all files

path <- "data/subsets"
files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
out_path <-  "data/prepped_for_pipeline"

for (f in files) {
  message("Processing: ", f)
  so <- readRDS(f)
  
  # store original cell number for later
  so@misc$orig_n_cells <- ncol(so@meta.data)


  #bring into usable shape
  so@meta.data <- so@meta.data %>%
    rename(
      sex = Sex, 
      tissue = Organ_name,
      age = Age_group,
      age_sex = Sex.Age_group,
     # nFeature_originalexp = nFeature_RNA,
     # nCount_originalexp = nCount_RNA,
      cell_ontology_class = Main_cell_type)

  if (DefaultAssay(so) != "originalexp") {
   so <- RenameAssays(so, RNA = "originalexp")
  }
  
  # TODO:extract cluster names for after the pipeline
  cell_type_df <- so@meta.data 
  rownames(cell_type_df) <- NULL

  cell_type_df <- cell_type_df |>
    select(tissue, age_sex, ID, cell_ontology_class, Lineage) |>
    unique() |>
    arrange(ID)
  cell_type_df
  # add save cell type annotation


  # save seurat object
  file_name <- paste0(sub("\\.rds$", "", basename(f)), "_prepped.rds")
  saveRDS(so, paste0(out_path, "/" , file_name))
  message("Saved: ", file_name)

}

