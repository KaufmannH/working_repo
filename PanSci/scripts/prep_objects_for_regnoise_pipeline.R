# This script was used to rename the metadata cols in the converted .rds PanSci object to the names used in the pipeine scripts


#library('tidyverse')
library('Seurat')
library('tidyverse')


# comparison: what is different that i need? 
c <- readRDS('/home/hkaufm49/analyses/RegulatedNoiseHel/data/rawdata/Spleen_rawdata.rds')
head(c@meta.data)
colnames(c@meta.data)

# example data
t <- readRDS('PanSci/data/subsets/Female.03_months.rds')
dim(t)
head(t@meta.data)
colnames(t@meta.data)



# do all files
path <- "PanSci/data/subsets"
files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
out_path <-  "PanSci/data/prepped_for_pipeline"

for (f in files) {
 # f <- files[3]
  message("Processing: ", f)
  so <- readRDS(f)
  
  # store original cell number for later
  so@misc$orig_n_cells <- ncol(so@meta.data)

  #bring into usable shape
  so@meta.data <- so@meta.data |>
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

  # filter for just WT mice
    so_filtered <- subset(so, subset = Genotype == "WT")

  # create new file name to make it compatiple with the pipeline and downstream processing
  tissue <- so_filtered@meta.data |>
    distinct(tissue) |>
    pull(tissue) 
  base <- tools::file_path_sans_ext(basename(f))
  sex <- base |>
    str_extract("^[^.]+") |>
    tolower()                                   
  age_num <- base |>
    str_extract("(?<=\\.)[0-9]{2}(?=_months)")   
  age_short <- paste0(as.integer(age_num), "m") 
  file_name <- paste0(tissue, "_", age_short, "_", sex, "_prepped.rds")

    # save seurat object
  saveRDS(so_filtered, paste0(out_path, "/" , file_name))
  message("Saved: ", file_name)

}



