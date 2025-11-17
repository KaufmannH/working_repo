library(Seurat)
library(tidyverse)

source('src/FunctionLibHel.R')

# deifning the cluster map

files <- list.files("data/processed/original_version/seurat_objects", pattern = "\\.rds$", full.names = TRUE)

read_object <- function(file_paths){
  so <- readRDS(file_paths)

  age_sex_group <- str_match(basename(file_paths), "^[^_]+_([^_]+_[^_]+)_processed\\.rds$")[,2]
  tibble(
    cell_id       = Cells(so),
    age_sex = age_sex_group,
    cluster       = as.character(Idents(so)))
}

cluster_info <- map_dfr(files, read_object)
saveRDS(cluster_info, "data/original_clusters/original_cluster_map.rds")



# overwriting clusters of an so object


# for many seurat files as input
cluster_map_raw <- readRDS("data/original_clusters/original_cluster_map.rds")
cluster_map <- cluster_map_raw %>%
  transmute(age_sex, cell_id, cluster = as.character(cluster)) %>%
  split(.$age_sex) %>%
  map(~ setNames(.x$cluster, .x$cell_id))


impose_og_clusters <- function(file_paths){
  so <- readRDS(file_paths)
  age_sex_group <- str_match(basename(file_paths),
                             "^[^_]+_([^_]+_[^_]+)_processed\\.rds$")[,2]

  og_cluster_vec <- cluster_map[[age_sex_group]]
  if (is.null(og_cluster_vec)) {
    warning("No old cluster map for group: ", age_sex_group, " (", basename(file_paths), ")")
    return(so)
  }

  aligned <- og_cluster_vec[Cells(so)]
  so$seurat_clusters_reprod <- Idents(so)
  so <- AddMetaData(so, aligned, col.name = "seurat_clusters")
  lvl <- unique(as.character(cluster_map_raw$cluster[cluster_map_raw$age_sex == age_sex_group]))
  so$seurat_clusters <- factor(so$seurat_clusters, levels = lvl)
  Idents(so) <- "seurat_clusters"
  return(so)
}

so_files <- list.files("data/processed/reproduced_vs1/seurat_objects", pattern = "\\.rds$", full.names = TRUE)
so_og_clusters <- map(so_files, impose_og_clusters)





# version for the script: only for one so

impose_og_clusters <- function(seurat_object, age_sex_group){

cluster_map_raw <- readRDS("data/original_clusters/original_cluster_map.rds")
cluster_map <- cluster_map_raw %>%
  transmute(age_sex, cell_id, cluster = as.character(cluster)) %>%
  split(.$age_sex) %>%
  map(~ setNames(.x$cluster, .x$cell_id))

  og_cluster_vec <- cluster_map[[age_sex_group]]
  aligned <- og_cluster_vec[Cells(seurat_object)]
  aligned <- setNames(as.character(aligned), Cells(seurat_object))
  seurat_object$seurat_clusters_reprod <- Idents(seurat_object)
  seurat_object <- AddMetaData(seurat_object, aligned, col.name = "seurat_clusters")
  lvl <- unique(as.character(cluster_map_raw$cluster[cluster_map_raw$age_sex == age_sex_group]))
  seurat_object$seurat_clusters <- factor(seurat_object$seurat_clusters, levels = lvl)
  Idents(seurat_object) <- "seurat_clusters"
  return(seurat_object)
}

so_og_clusters <- impose_og_clusters(so, age_sex)


age_sex <- '1m_male'
so <- readRDS("data/processed/reproduced_vs2/seurat_objects/Spleen_1m_male_processed.rds")