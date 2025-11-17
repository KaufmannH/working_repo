# required packages
library(tidyverse)
library(Seurat)
library(patchwork)
#library(scDblFinder)
library(scater)
library(philentropy)
library(foreach)
library(doFuture)




MakeEmbsClusters.df <- function(SeuratObj , EmbsType= 'pca')
{
  # It makes a data.frame which contains Embeddings & Clusters
  Clusters <- SeuratObj[['seurat_clusters']]
  Embs     <- Embeddings(object = SeuratObj[[EmbsType]])
  EmbsClusters <- merge(Clusters,Embs,by=0, all=TRUE)[,-1]
  rownames(EmbsClusters) <- rownames(Embs)
  return(EmbsClusters)
}


FindCentroids      <- function(EmbsClusters)
{
  # It compute the centroids in each cluster
  Centroids <- aggregate(as.matrix(EmbsClusters[,-1]) ~ seurat_clusters,EmbsClusters,mean)
  return(Centroids)
}


FindUnstableCells  <- function(EmbsClusters, Centroids)
{
  # It compute the distance between each cell and each cluster, 
  # if the minimal distance is NOT with the belonging cluster, than it is an unstable cell.
  
  Embs_euDis = data.frame( cell_id = c(),
                           cell_cluster = c(),
                           euclidian = c(), 
                           centroid_cluster=c())
  
  rowCounter = 0
  
  for (cell in 1: nrow(EmbsClusters)){
    
    if (cell %% 100 == 0) {
      
      print(paste(cell, "cells processed"))
      
    }
    
    
    
    cell_emb <- EmbsClusters[cell,-1]
    
    df <- data.frame( cell_id = c(),
                      cell_cluster = c(),
                      euclidian = c(), 
                      centroid_cluster=c())
    
    for (c in 1:nrow(Centroids[,-1]))
    {
      
      rowCounter <- rowCounter + 1
      philentropy:: distance(rbind(cell_emb, Centroids[c,-1]), mute.message = TRUE) -> euDis
      
      d <- data.frame( cell_id = rownames(cell_emb),
                       cell_cluster = EmbsClusters[cell,1],
                       euclidian = euDis ,
                       centroid_cluster = Centroids[c,1],
                       row.names = rowCounter)
      
      df <- rbind(df,d)
      
    }
    df <- df[df$euclidian == min(df$euclidian),]
    Embs_euDis <- rbind(Embs_euDis,df)
  }
  UnstableCells <- Embs_euDis[Embs_euDis$cell_cluster != Embs_euDis$centroid_cluster,]
  return(UnstableCells)
}


AssignUnstableCellFlag <- function(so){
  
  EmbsClusters  <- MakeEmbsClusters.df(so , EmbsType= 'pca')
  Centroids     <- FindCentroids(EmbsClusters)
  
  UnstableCells <- FindUnstableCells(EmbsClusters, Centroids)
  
  so@meta.data$UnstableCells <- row.names(so@meta.data) %in% UnstableCells$cell_id
  return(so)
  
}


# This function takes as input a seurat object and a logical vector for defining a subset of cells
# and calculates the residual variance and geometric mean for each gene across the subset of cells.
# It requires an SCT model and a populated vst.out slot in so[["SCT"]]@misc$vst.out
res_var_subset <- function(so, cells_filter, count_slot_name, selected_cluster){
  #meta_col <- 'seurat_clusters'
  #cells_filter <- so@meta.data[[meta_col]] == '1'
  #count_slot_name <- "originalexp"
  #selected_cluster <- "18"


  # Extract the SCT model object
  ref_model <- SCTResults(so, slot = "model")
  # subset
  sub <- subset(so, cells = Cells(so)[cells_filter])


  # calc residual variance per gene across subset
  sub_res_var <- FetchResiduals(
    sub,
    assay = "SCT",
    reference.SCT.model = ref_model,
    features = rownames(sub),
    replace.value = TRUE,
    verbose = FALSE)


  gmean_mat   <- GetAssayData(sub, assay = "SCT", slot = "data")

  # calc the resvar across cells myself
  res_var_vec <- matrixStats::rowVars(as.matrix(sub_res_var), na.rm = TRUE)
  names(res_var_vec) <- rownames(sub_res_var)

  common_genes <- intersect(rownames(sub_res_var), rownames(gmean_mat))
  common_genes <- common_genes[order(common_genes)]

  gmean_vec <- rowMeans(gmean_mat[common_genes, , drop = FALSE])

  res_var_mean <- tibble::tibble(
    gene    = common_genes,
    res_var = unname(res_var_vec[common_genes]),
    gmean   = unname(gmean_vec),
    cluster = selected_cluster,
    hvg     = res_var_vec[common_genes] > 5)


  return(res_var_mean)
}


# This function takes as input a seurat object and a metadata column. It will reestimate
# the variability and gmean of all genes across the cells in each group (defined by metadata column)
# and return results as a tibble. meta_col should be a factor
res_var_by_group <- function(so, meta_col, hvg_cutoff = 5, count_slot_name, n_cores, vst_out = NULL, integrated = FALSE){
  
  #cluster <- makeCluster(n_cores)
  #registerDoParallel(cluster)

  res_var_tbl <- 
    foreach(cl = levels(so@meta.data[[meta_col]]), .combine = 'rbind', .export = ls(globalenv())) %do% {

      cells_filter <- so@meta.data[[meta_col]] == cl
      print(paste("cluster: ", cl))

      n_cells <- sum(cells_filter)
      print(paste('number of cells: ', n_cells ))
      if (n_cells < 2) {
        message("WARNING: Skipping cluster ", cl, " — only ", n_cells, " cell(s)")
        return(return(tibble::tibble(gene = character(), res_var = numeric(),
                            gmean = numeric(), cluster = character())))}

      res_var_mean_subs <- res_var_subset(so, cells_filter, count_slot_name, selected_cluster = cl)
      
      # convert to tibble and add cluster info
        res_var_tbl <- res_var_mean_subs %>%
            mutate(cluster = cl)
    }
  
  #stopCluster(cluster)
  
  res_var_tbl <- res_var_tbl %>%
    mutate(hvg = ifelse(res_var > hvg_cutoff, TRUE, FALSE)) # define hypervariable genes
  
  return(res_var_tbl)
  
}



# This function takes as an input a seurat object, clusters and a number of cycles and will
# then perform that many cycles of bootstrapping by choosing random subsets of cells
# and re-calculating the variability for each gene in those subsets
res_var_bootstrap <- function(so, clusters_slot, n_cycles = 1000, n_cores, hvg_cutoff = 5, res_var_cl, keep_all = FALSE){
  #meta_col <- 'seurat_clusters'
  #cells_filter <- so@meta.data[[meta_col]] == '1'
  #count_slot_name <- "originalexp"
  #selected_cluster <- "5"
 # clusters_slot <- 'seurat_clusters'
 # n_cores <- 1
 # n_cycles = 1000
 # cl_vec <- c("4")

  sct_assay <- 'SCT'

  # Use the exact SCT model that was fit on 'so'
  ref_model <- SCTResults(so, slot = "model")

  # Which counts assay was used when fitting SCT
  umi_assay <- tryCatch(SCTResults(so, slot = "umi.assay"), error = function(e) "originalexp")

  cl_vec <- so@meta.data[[clusters_slot]]
  if (!is.factor(cl_vec)) cl_vec <- factor(cl_vec)

  # parallel setup (optional)
  if (n_cores > 1) {
    doParallel::registerDoParallel(n_cores)
    `%op%` <- `%dopar%`
  } else {
    `%op%` <- `%do%`
  }

  # one bootstrap run on a given cell vector -> same output schema as your function
  .one_iter <- function(cell_names, cl, iter_idx){
    # subset + recompute residuals using the SAME model
    sub <- subset(so, cells = cell_names)
    DefaultAssay(sub) <- umi_assay

    sub_res_var <- FetchResiduals(
      sub,
      assay = sct_assay,
      reference.SCT.model = ref_model,
      features = rownames(sub),
      replace.value = TRUE,
      verbose = FALSE
    )

    # residuals (genes x cells)
    sub_res_var <- GetAssayData(sub, assay = sct_assay, slot = "scale.data")

    # gmean source (same as your non-bootstrap: from SCT@data)
    gmean_mat <- GetAssayData(sub, assay = sct_assay, slot = "data")

    # res_var across cells (same method you used)
    res_var_vec <- matrixStats::rowVars(as.matrix(sub_res_var), na.rm = TRUE)
    names(res_var_vec) <- rownames(sub_res_var)

    # align genes and compute gmean
    common_genes <- intersect(rownames(sub_res_var), rownames(gmean_mat))
    common_genes <- common_genes[order(common_genes)]
    gmean_vec    <- rowMeans(gmean_mat[common_genes, , drop = FALSE])

    tibble::tibble(
      gene    = common_genes,
      res_var = unname(res_var_vec[common_genes]),
      gmean   = unname(gmean_vec),
      cluster = cl,
      boolean_hvg     = res_var_vec[common_genes] > hvg_cutoff,
      iter    = iter_idx
    )
  }

  out_all <- list()
  for (cl in levels(cl_vec)) {
    message("cluster: ", cl)
    cells_filter <- cl_vec == cl
    cell_names   <- colnames(so)[cells_filter]
    n_cells      <- length(cell_names)
    message("number of cells: ", n_cells)

    if (n_cells < 2) {
      warning("Skipping cluster ", cl, " — only ", n_cells, " cell(s)")
      next
    }

    set.seed(42)
    # bootstrap: sample same number of cells with replacement, each iteration
    res_cl <- foreach::foreach(iter = seq_len(n_cycles), .combine = "bind_rows") %op% {
      selected.cells <- sample(cell_names, size = n_cells, replace = TRUE)
      .one_iter(selected.cells, cl, iter)
    }

    out_all[[cl]] <- res_cl
    gc(verbose = FALSE)
  }

  if (length(out_all) == 0) {
    return(tibble::tibble(
      gene = character(), res_var = numeric(), gmean = numeric(),
      cluster = character(), boolean_hvg = logical(), iter = integer()
    ))
  }

  dplyr::bind_rows(out_all)
}



plot_resvar <- function(res_var_cl, bootstrap_res, test_gene, housek_gene, plot_cluster) {
  
  resvar_test_gene <- res_var_cl %>% 
    filter(cluster == plot_cluster, gene == test_gene) %>% pull(res_var)
  
  resvar_plot <- ggplot(subset(bootstrap_res, gene %in% c(test_gene, housek_gene) & cluster == plot_cluster), 
                        aes(x = ResVar, y = gene)) +
    geom_boxplot() +
    geom_point(aes(x = resvar_test_gene, y = test_gene), colour = "red", size = 5, shape = 16)
  
  return(resvar_plot)
}




# An adaptation of res_var_bootstrap to calculate stability for stable genes
# This function takes as an input a seurat object, clusters and a number of cycles and will
# then perform that many cycles of bootstrapping by choosing random subsets of cells
# and re-calculating the variability for each gene in those subsets
res_var_bootstrap_low <- function(so, clusters_slot, n_cycles = 1000, n_cores, lvg_cutoff = 1, res_var_cl, keep_all = FALSE){
  
  n_cl <- 0
  
  for (cl in levels(so@meta.data[[clusters_slot]])){
    
    print(cl)
    
    # Get Number of Hyper Variable Genes
    features_tbl <- rownames(so@assays$originalexp@meta.features) %>% 
      as_tibble() %>% dplyr::rename(features = value)
    
    # 1. Cell filter for cluster
    cells_filter <- so@meta.data[[clusters_slot]] == cl 
    cell_names <- colnames(so)[cells_filter]
    
    # 2. Sample Cells in each Cluster 1000 with replacement
    
    set.seed(42) # Set the seed to be able to reproduce the same set of sampling
    
    # setup a cluster for parallelization 
    #cluster <- makeCluster(n_cores)
    #registerDoParallel(cluster)
    
    #registerDoFuture()
    #plan(multisession, workers = n_cores)
    
    out.matrix <- foreach(iter = 1:n_cycles, .combine = 'rbind', .export = ls(globalenv())) %do% {
      
      selected.cells <- sample(cell_names, size = length(cell_names), replace = TRUE)
      
      umi.subset     <- so@assays$originalexp@counts[,selected.cells]
      vst_out.subset <- so@assays$SCT@misc$vst.out
      vst_out.subset$cell_attr <- vst_out.subset$cell_attr[selected.cells,]
      
      ResVar         <- sctransform::get_residual_var(vst_out = vst_out.subset, umi = umi.subset)
      
      # build output matrix
      features <- names(ResVar)
      boolean_lvg <- features %in% names(ResVar[ResVar <= lvg_cutoff])
      
      cbind(boolean_lvg, ResVar, iter)
      
    }
    
    #stopCluster(cluster)
    
    # lvgs for cluster
    lvgs <- res_var_cl %>% filter(cluster == cl,
                                  res_var <= lvg_cutoff)
    
    # filter results for lvgs
    if (keep_all) {
      tmp.tibble <- as_tibble(out.matrix, rownames = "gene")
    } else {
      tmp.tibble <- as_tibble(out.matrix, rownames = "gene") %>%
        filter(gene %in% lvgs$gene)
    }
    
    rm(out.matrix)
    gc()
    
    tmp.tibble$cluster <- cl
    
    # combine result tibbles
    if (n_cl < 1) {
      out.tibble <- tmp.tibble
    } else {
      out.tibble <- bind_rows(out.tibble, tmp.tibble)
    }
    n_cl <- n_cl + 1
    
  }
  
  return(out.tibble)
}



# function for sampling LVGs to compare with HVGs (still needs weighting adaptation)
sample_LVG <- function(res_var_cl, HVG_cutoff = 5, LVG_cutoff = 1) {
  
  count = 0
  
  for (cl in unique(res_var_cl$cluster)) {
    
    hvg <- res_var_cl %>% filter(cluster == cl, res_var >= HVG_cutoff)
    lvg <- res_var_cl %>% filter(cluster == cl, res_var <= LVG_cutoff)
    
    lvg_sample <- sample(lvg$gene, size = nrow(hvg), replace = FALSE)
    lvg_subset <- lvg %>% filter(gene %in% lvg_sample)
    
    if (count == 0) {
      lvg_subset_tbl <- lvg_subset
    } else {
      lvg_subset_tbl <- rbind(lvg_subset_tbl, lvg_subset)
    }
    
    count = count + 1
    
  }
  
  return(lvg_subset_tbl)
  
}