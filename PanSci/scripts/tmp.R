# rename the finished pansic objects and their age_sex col
library(tidyverse)
library(Seurat)

path <- 'RegulatedNoiseHel/data/processed/pansci_liver/seurat_objects/'
files <- list.files(path, full.names = TRUE)


for (f in files){
    f_path <- f
    print(f_path)
    obj <- readRDS(f_path)

    obj@meta.data <- obj@meta.data |>
        mutate(
        age_sex = str_match(age_sex, "(\\d+)_months_([A-Za-z]+)") |>
        (\(x) paste0(as.integer(x[,2]), "m_", tolower(x[,3])))())
    print(unique( obj@meta.data$age_sex))
    saveRDS(obj, f_path)
}

test <- readRDS(files[1])
unique(test@meta.data$age_sex)


