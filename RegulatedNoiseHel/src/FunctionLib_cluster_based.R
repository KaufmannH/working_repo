# required packages
library(tidyverse)
library(Seurat)
library(patchwork)
#library(scDblFinder)
library(scater)
library(philentropy)
library(foreach)
library(doFuture)
library(doRNG)
library(doParallel)


 

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
  #cells_filter <- so@meta.data[[meta_col]] == '6'
  #count_slot_name <- "originalexp"
  #selected_cluster <- "6"

  # subset
  sub <- subset(so, cells = Cells(so)[cells_filter])

  # rerun SCT
  sub <- SCTransform(sub, assay = "originalexp",  
                          verbose = FALSE)
   

  # residual variance per gene across subset
  sub_res_var <- GetAssayData(sub, assay = "SCT", slot = "scale.data")
  # normalized expression per gene across subset
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
  sct_assay <- "SCT"
  
  # pull residuals and expression
  residuals_mat <- GetAssayData(so, assay = sct_assay, layer = "scale.data")
  if (is.null(residuals_mat) || nrow(residuals_mat) == 0L) stop("No residuals in SCT/scale.data.")
  data_mat <- GetAssayData(so, assay = sct_assay, layer = "data")
  if (is.null(data_mat) || nrow(data_mat) == 0L) stop("No data layer in SCT/data.")

  # align genes
  genes <- intersect(rownames(residuals_mat), rownames(data_mat))
  if (length(genes) == 0L) stop("No overlapping genes between residuals and data.")
  residuals_mat <- residuals_mat[genes, , drop = FALSE]
  data_mat      <- data_mat[genes, , drop = FALSE]

  # get clusters
  cl_vec <- so@meta.data[[clusters_slot]]
  if (!is.factor(cl_vec)) cl_vec <- factor(cl_vec)
  cl_levels <- levels(cl_vec)

  # parallelise across clusters (register once)
  if (n_cores > 1L) {
    doParallel::registerDoParallel(n_cores)
    `%op%` <- `%dopar%`
  } else {
    `%op%` <- `%do%`
  }
  set.seed(42); doRNG::registerDoRNG(42)

  # helper: vectorised bootstrap inside a cluster using counts matrix
  .bootstrap_cluster <- function(residuals_cl, data_cl, B){
    n_cells <- ncol(residuals_cl)

    # counts matrix: each column is a bootstrap draw (sums to n_cells)
    counts_mat <- rmultinom(n = B, size = n_cells, prob = rep(1 / n_cells, n_cells))  # n_cells x B (integer)

    # per-gene sums and sums-of-squares across bootstraps (genes x B)
    sum_residuals      <- residuals_cl %*% counts_mat
    sumsq_residuals    <- (residuals_cl * residuals_cl) %*% counts_mat

    # means and (sample) variances of residuals
    mean_residuals     <- sum_residuals / n_cells
    var_mle_residuals  <- (sumsq_residuals / n_cells) - (mean_residuals * mean_residuals)
    if (n_cells > 1L) {
      var_residuals <- var_mle_residuals * (n_cells / (n_cells - 1))
    } else {
      var_residuals <- var_mle_residuals
    }

    # gmean from data layer (same counts)
    sum_data  <- data_cl %*% counts_mat
    mean_data <- sum_data / n_cells

    list(res_var_mat = var_residuals, gmean_mat = mean_data)
  }

  # run per cluster (in parallel)
  result_list <-
    foreach(cl = cl_levels, .combine = "rbind",
            .packages = c("tibble")) %op% {
      cell_idx <- which(cl_vec == cl)
      n_cells  <- length(cell_idx)
      message("Processing cluster ", cl, " (", n_cells, " cell(s)).")
      if (n_cells < 2L) {
        message("Skip cluster ", cl, " (", n_cells, " cell(s)).")
        return(NULL)
      }

      # slice matrices for cluster
      residuals_cl <- residuals_mat[, cell_idx, drop = FALSE]
      data_cl      <- data_mat[,      cell_idx, drop = FALSE]

      # vectorised bootstrap for this cluster
      boot <- .bootstrap_cluster(residuals_cl, data_cl, n_cycles)

      # assemble long table without per-iteration loops
      # (genes x n_cycles) → long
      res_var_vec <- as.vector(boot$res_var_mat)  # column-major: iter1, iter2, ...
      gmean_vec   <- as.vector(boot$gmean_mat)

      # identifiers (use rep(), not rep.int())
      iter_ids  <- rep(seq_len(n_cycles), each = length(genes))
      gene_ids  <- rep(genes,             times = n_cycles)
      bool_vec  <- res_var_vec > hvg_cutoff

      tibble(
        gene        = gene_ids,
        res_var     = res_var_vec,
        gmean       = gmean_vec,
        cluster     = cl,
        boolean_hvg = bool_vec,
        iter        = iter_ids)
    }

  # bind results
  if (is.null(result_list) || nrow(result_list) == 0L) {
    tibble(
      gene = character(), res_var = numeric(), gmean = numeric(),
      cluster = character(), boolean_hvg = logical(), iter = integer()
    )
  } else {
    result_list
  }

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