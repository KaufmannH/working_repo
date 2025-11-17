# additional functions added by Helene

merge_processed_seurat_objects <- function(path) {

    files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
    objs <- lapply(files, readRDS)
    ids <- sub("^.*_(\\d+m_(male|female)).*$", "\\1", basename(files))
    merged_so <- merge(x = objs[[1]], y = objs[-1], add.cell.ids = ids, project = "Spleen_merged")
    merged_so$age_sex <- sapply(strsplit(colnames(merged_so), "_"), `[`, 1)

    return(merged_so)
}

add_resvar_to_seurat_object <- function(seurat_object, path){
    files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
    ids <- sub("^.*_(\\d+m_(male|female)).*$", "\\1", basename(files))
    #seurat_object@misc$resvar_tables <- setNames(lapply(files, readRDS), ids)
    seurat_object@misc$resvar_table <- dplyr::bind_rows(Map(\(f, id) dplyr::mutate(readRDS(f), age_sex = id), files, ids))

    return(seurat_object)
}

compare_res_var_distributions <- function(df_old, df_new, out_path) {
    # variables for log plotting
    bin_width <- 0.05
    safety_log <- 1e-6
    old_log <- log10(df_old$res_var + safety_log)
    new_log <- log10(df_new$res_var + safety_log)
    log_brks <- seq(floor(min(old_log, new_log, na.rm = TRUE)), ceiling(max(old_log, new_log, na.rm = TRUE)), by = bin_width)

    plot_df_old <- df_old |>
        filter(gmean != 0) |>
        mutate(log_res_var = log10(res_var + 1e-6)) |>
        filter(is.finite(log_res_var)) |> 
        #mutate(bin = cut(res_var, breaks = seq(0, ceiling(max(df_old$res_var, df_old$res_var, na.rm = TRUE) / bin_width) * bin_width, by = bin_width))) |>
         mutate(bin = cut(log_res_var, breaks = log_brks, include.lowest = TRUE)) |>
        group_by(bin) |>
        summarize(n_genes = n())
        plot_df_old

    plot_df_new <- df_new |>
        filter(gmean != 0) |>
        mutate(log_res_var = log10(res_var + 1e-6)) |>
        filter(is.finite(log_res_var)) |> 
        #mutate(bin = cut(res_var, breaks = seq(0, ceiling(max(df_new$res_var, df_new$res_var, na.rm = TRUE) / bin_width) * bin_width, by = bin_width))) |>
         mutate(bin = cut(log_res_var, breaks = log_brks, include.lowest = TRUE)) |>
       # count(bin, name = "n_genes")
        group_by(bin) |>
        summarize(n_genes = n())

    plot_new <- ggplot(plot_df_new, aes(x = bin, y = n_genes)) +
        geom_col(fill = '#5195D3') +
        labs(x = "Residual variance bins (log10)",
            y = "Gene count", title = "New") +
    theme_classic() +
    theme(
        strip.background      = element_blank(),
        axis.text.x           = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y           = element_text(size = 20),
        axis.title.x      = element_text(size = 20),
        axis.title.y      = element_text(size = 20)) +
    scale_x_discrete(breaks = levels(plot_df_new$bin)[seq(1, length(levels(plot_df_new$bin)), by = 10)])


    plot_old <- ggplot(plot_df_old, aes(x = bin, y = n_genes)) +
            geom_col(fill = '#5195D3') +
            labs(x = "Residual variance bins (log10)",
                y = "Gene count", title = "Old") +
        theme_classic() +
        theme(
            strip.background      = element_blank(),
            axis.text.x           = element_text(angle = 45, hjust = 1, size = 15),
            axis.text.y           = element_text(size = 20),
            axis.title.x      = element_text(size = 20),
            axis.title.y      = element_text(size = 20)) +
        scale_x_discrete(breaks = levels(plot_df_old$bin)[seq(1, length(levels(plot_df_old$bin)), by = 10)])

    combined_plot <- plot_old | plot_new 
    ggsave(paste0(out_path, "/res_var_dist.png"),  combined_plot, width = 15, height = 5)

}


correlating_res_var <- function(df_old, df_new, out_path){

    # get gene + cluster unique identifier
    df_old_edit <- df_old |>
        mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
        mutate(res_var_old = res_var) |>
        select(gene_cluster_age_sex, res_var_old)
  
    df_new_edit <- df_new |>
        mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
        mutate(res_var_new = res_var) |>
        select(gene_cluster_age_sex, res_var_new)

    intersect_df <- df_old_edit |>
        inner_join(df_new_edit, by = "gene_cluster_age_sex") |>
        mutate(delta = res_var_new - res_var_old)

    extreme_cases_df <- intersect_df |>
        filter(abs(res_var_old - res_var_new) > 5) |>
        arrange(gene_cluster_age_sex)
    dim(intersect_df)

    plot <- ggplot(intersect_df, aes(x = res_var_old, y = res_var_new, color = delta )) +
        geom_point(alpha = 0.3) +
        geom_abline(slope = 1, intercept = 0, colour = "lightblue") +
        labs(x = "Residual variance (old run)", y = "Residual variance (new run)") +
        scale_colour_gradient2(
        low = "black", mid = "grey70", high = "#1a9850",
        midpoint = 0,
        name = "Δ new–old" ) +
        theme_classic(base_size = 20) +                    
        theme(
            axis.title  = element_text(size = 20),
            axis.text   = element_text(size = 20),
            legend.title= element_text(size = 20),
            legend.text = element_text(size = 18))

    ggsave(paste0(out_path, "/res_var_corr.png"),  plot, width = 10, height = 10)

}


correlating_res_var_heat <- function(df_old, df_new, out_path){

    df_old_edit <- df_old |>
        mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
        mutate(res_var_old = res_var) |>
        select(gene_cluster_age_sex, res_var_old)
  
    df_new_edit <- df_new |>
        mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
        mutate(res_var_new = res_var) |>
        select(gene_cluster_age_sex, res_var_new)

    intersect_df <- df_old_edit |>
        inner_join(df_new_edit, by = "gene_cluster_age_sex") |>
        mutate(delta = res_var_new - res_var_old)

    extreme_cases_df <- intersect_df |>
        filter(abs(res_var_old - res_var_new) > 5) |>
        arrange(gene_cluster_age_sex)
    dim(intersect_df)

    plot <- ggplot(intersect_df, aes(x = res_var_old, y = res_var_new)) +
        ggpointdensity::geom_pointdensity(alpha = 0.4) +
        ggpointdensity::geom_pointdensity(aes(colour = log10(..density..)), size = 1, alpha = 0.4) +
        geom_abline(slope = 1, intercept = 0, colour = "lightblue") +
        labs(x = "Residual variance (old run)", y = "Residual variance (new run)") +
         #scale_colour_distiller(palette = "YlGnBu", direction = -1, name = "Point density") +
           scale_colour_distiller(
            palette = "YlGnBu",
            direction = -1,
            name = "log10(density)"
        ) +
        theme_classic(base_size = 20) +                    
        theme(
            axis.title  = element_text(size = 20),
            axis.text   = element_text(size = 20),
            legend.title= element_text(size = 20),
            legend.text = element_text(size = 18))


    ggsave(paste0(out_path, "/res_var_corr_heat.png"),  plot, width = 10, height = 10)

}



correlating_res_var_log <- function(df_old, df_new, out_path){

    # get gene + cluster unique identifier
    df_old_edit <- df_old |>
        mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
        mutate(res_var_old = res_var) |>
        select(gene_cluster_age_sex, res_var_old)
  
    df_new_edit <- df_new |>
        mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
        mutate(res_var_new = res_var) |>
        select(gene_cluster_age_sex, res_var_new)

    intersect_df <- df_old_edit |>
        inner_join(df_new_edit, by = "gene_cluster_age_sex") |>
        mutate(delta = res_var_new - res_var_old)

    extreme_cases_df <- intersect_df |>
        filter(abs(res_var_old - res_var_new) > 5) |>
        arrange(gene_cluster_age_sex)
    dim(intersect_df)

    plot <- ggplot(intersect_df, aes(x = res_var_old, y = res_var_new, color = delta )) +
        geom_point(alpha = 0.3) +
        geom_abline(slope = 1, intercept = 0, colour = "lightblue") +
        labs(x = "Residual variance log10 (old run)", y = "Residual variance log10 (new run)") +
        scale_colour_gradient2(
        low = "black", mid = "grey70", high = "#1a9850",
        midpoint = 0,
        name = "Δ new–old" ) +
        scale_x_log10(labels = scales::label_log()) +
        scale_y_log10(labels = scales::label_log()) +
        theme_classic(base_size = 20) +                    
        theme(
            axis.title  = element_text(size = 20),
            axis.text   = element_text(size = 20),
            legend.title= element_text(size = 20),
            legend.text = element_text(size = 18))

    ggsave(paste0(out_path, "/res_var_corr_log.png"),  plot, width = 10, height = 10)
}



correlating_res_var_log_heat <- function(df_old, df_new, out_path){
   # individual points are colored by the number of neighboring points

     # get gene + cluster unique identifier
    df_old_edit <- df_old |>
        mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
        mutate(res_var_old = res_var) |>
        select(gene_cluster_age_sex, res_var_old)
  
    df_new_edit <- df_new |>
        mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
        mutate(res_var_new = res_var) |>
        select(gene_cluster_age_sex, res_var_new)

    intersect_df <- df_old_edit |>
        inner_join(df_new_edit, by = "gene_cluster_age_sex") |>
        mutate(delta = res_var_new - res_var_old)

    extreme_cases_df <- intersect_df |>
        filter(abs(res_var_old - res_var_new) > 5) |>
        arrange(gene_cluster_age_sex)
    dim(intersect_df)

    plot <- ggplot(intersect_df, aes(x = res_var_old, y = res_var_new)) +
        #ggpointdensity::geom_pointdensity(size = 1, alpha = 0.4) +
        ggpointdensity::geom_pointdensity(aes(colour = log10(..density..)), size = 1, alpha = 0.4) +
        geom_abline(slope = 1, intercept = 0, colour = "lightblue") +
        labs(x = "Residual variance log10 (old run)", y = "Residual variance log10 (new run)") +
        scale_colour_distiller(palette = "YlGnBu", direction = -1, name = "Point density") +
        scale_x_log10(labels = scales::label_log()) +
        scale_y_log10(labels = scales::label_log()) +
        theme_classic(base_size = 20) +                    
        theme(
            axis.title  = element_text(size = 20),
            axis.text   = element_text(size = 20),
            legend.title= element_text(size = 20),
            legend.text = element_text(size = 18))

    ggsave(paste0(out_path, "/res_var_corr_log_heat.png"),  plot, width = 10, height = 10)
}


umap_cluster_comparison <- function (age_sex_group, old_data_path, new_data_path, out_path){
    # get old object
    old_files <- list.files(
                file.path(old_data_path, "seurat_objects"),
                pattern = "\\.rds$",
                full.names = TRUE)
    old_ids <- sub("^.*_(\\d+m_(male|female))_processed\\.rds$", "\\1", basename(old_files))

    match_idx <- match(age_sex_group, old_ids)
    if (is.na(match_idx)) {
        stop(paste("No Seurat object found for", age_sex_group))}
    seurat_obj_old <- readRDS(old_files[match_idx])

    # get new object
    new_files <- list.files(
                file.path(new_data_path, "seurat_objects"),
                pattern = "\\.rds$",
                full.names = TRUE)
    new_ids <- sub("^.*_(\\d+m_(male|female))_processed\\.rds$", "\\1", basename(new_files))

    match_idx <- match(age_sex_group, new_ids)
    if (is.na(match_idx)) {
        stop(paste("No Seurat object found for", age_sex_group))}
    seurat_obj_new <- readRDS(new_files[match_idx])


    p1 <- UMAPPlot(seurat_obj_old, group.by = "seurat_clusters") 
    p2 <- UMAPPlot(seurat_obj_old, group.by = "cell_ontology_class")
    q1 <- UMAPPlot(seurat_obj_new, group.by = "seurat_clusters")
    q2 <- UMAPPlot(seurat_obj_new, group.by = "cell_ontology_class")

    combined <-(p1 | p2) /
              (q1 | q2) +
            patchwork::plot_annotation(
            title    = stringr::str_replace_all(age_sex_group, "_", " "),
            subtitle = "Top: old  •  Bottom: reproduced",
            theme    = ggplot2::theme(
            plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5)))

    ggsave(paste0(out_path, "/umap_cluster_", age_sex_group, ".png"), combined, width = 12, height = 10)
}



umap_cluster_comparison_grey <- function (age_sex_group, old_data_path, new_data_path, out_path){

    # load old 
    old_files <- list.files(
    file.path(old_data_path, "seurat_objects"),
    pattern = "\\.rds$",
    full.names = TRUE)

    old_ids <- sub("^.*_(\\d+m_(male|female))_processed\\.rds$", "\\1", basename(old_files))
    match_idx <- match(age_sex_group, old_ids)
    if (is.na(match_idx)) stop(paste("No Seurat object found for (OLD)", age_sex_group))
    seurat_obj_old <- readRDS(old_files[match_idx])

    # load new
    new_files <- list.files(
    file.path(new_data_path, "seurat_objects"),
    pattern = "\\.rds$",
    full.names = TRUE)

    new_ids <- sub("^.*_(\\d+m_(male|female))_processed\\.rds$", "\\1", basename(new_files))
    match_idx <- match(age_sex_group, new_ids)
    if (is.na(match_idx)) stop(paste("No Seurat object found for (NEW)", age_sex_group))
    seurat_obj_new <- readRDS(new_files[match_idx])


    # extract dfs
    get_df <- function(so, dataset_label) {
    cl <- so@meta.data[['seurat_clusters']]
    if (!is.factor(cl)) cl <- factor(cl)
    emb <- Embeddings(so, reduction = "umap")
    tibble(
        UMAP_1  = emb[, 1],
        UMAP_2  = emb[, 2],
        cluster = cl,
        dataset = dataset_label )}

    df_old <- get_df(seurat_obj_old, "Old")
    df_new <- get_df(seurat_obj_new, "New")

    all_levels <- sort( unique(c( as.integer((df_old$cluster)), as.integer((df_new$cluster) ))))

    df_old <- df_old %>%
        mutate(cluster = factor(as.integer(cluster), levels = all_levels)) |> arrange(cluster)
    df_new <- df_new %>%
        mutate(cluster = factor(as.integer(cluster), levels = all_levels)) |> arrange(cluster)

    cluster_cols <- setNames(hue_pal(l = 55, c = 100)(length(all_levels)), all_levels)


    plot_df <- bind_rows(df_old, df_new) %>%
    mutate(
        row_id      = row_number(),
        cluster_chr = as.character(cluster) ) %>%
    tidyr::crossing(facet_cluster = all_levels) %>%
    mutate(
        facet_cluster_chr = as.character(facet_cluster),
        is_target         = cluster_chr == facet_cluster_chr,
        colour_hex        = if_else(is_target, cluster_cols[facet_cluster_chr], '#7F7F7F'),
        alpha_value       = if_else(is_target, 0.9, 0.2),
        dataset           = factor(dataset, levels = c("Old", "New")))

        x_rng <- range(plot_df$UMAP_1, na.rm = TRUE)
        y_rng <- range(plot_df$UMAP_2, na.rm = TRUE)


    p_facets <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(colour = colour_hex, alpha = alpha_value), size = 0.7, stroke = 0) +
    scale_color_identity() +
    scale_alpha_identity() +
    facet_grid(rows = vars(dataset), cols = vars(facet_cluster)) +
    coord_cartesian(xlim = x_rng, ylim = y_rng, expand = TRUE) +
    labs(
        x = "UMAP 1",
        y = "UMAP 2"  ) +
    theme_classic(base_size = 20) +
    theme(
        strip.text      = element_text(face = "bold"),
        strip.background= element_blank(),
        legend.position = "none")

    ggsave(
    file.path(out_path, paste0("umap_grey_", gsub(" ", "_", age_sex_group), ".png")),
    p_facets, width = 35, height = 5, dpi = 300)

}




gene_var_group_comparison <- function (df_old, df_new, out_path){
    
    df_old_tagged <- df_old |>
    mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
    mutate(gene_variability = case_when(
            gmean == 0 ~ NA_character_,
            res_var < 1 ~ "LVG",
            res_var > 5 ~ "HVG",
            TRUE        ~ "Intermediate" ))

    df_new_tagged <- df_new |>
    mutate(gene_cluster_age_sex = paste0(gene, "_", cluster, "_", age_sex)) |>
    mutate(gene_variability = case_when(
            gmean == 0 ~ NA_character_,
            res_var < 1 ~ "LVG",
            res_var > 5 ~ "HVG",
            TRUE        ~ "Intermediate" ))

    levels <- c("LVG","Intermediate","HVG")

    comparison <- df_old_tagged |>
    select(gene_cluster_age_sex, gene_variability_old = gene_variability, age_sex) |>
    inner_join(
        df_new_tagged |> select(gene_cluster_age_sex, gene_variability_new = gene_variability, age_sex),
        by = c("gene_cluster_age_sex", "age_sex")) |>
    mutate(
        gene_variability_old = factor(gene_variability_old, levels = levels),
        gene_variability_new = factor(gene_variability_new, levels = levels),
        is_match = case_when(
        is.na(gene_variability_old) | is.na(gene_variability_new) ~ NA,
        gene_variability_old == gene_variability_new ~ TRUE,
        TRUE ~ FALSE ))


        agree_count <- comparison |>
        filter(!is.na(is_match)) |>
        mutate(match_lab = if_else(is_match, "Match", "Mismatch")) |>
        dplyr::count(age_sex, match_lab, name = "n") |>
        group_by(age_sex) |>
        mutate(p = n / sum(n)) |>
        ungroup()


        match_plot <- ggplot(agree_count, aes(x = age_sex, y = p, fill = match_lab)) +
        geom_col() +
        geom_text(aes(label = scales::percent(p, accuracy = 0.1)),
                    position = position_stack(vjust = 0.5), size = 6) +
        scale_fill_manual(values = c("Match" = "#1b9e77", "Mismatch" = "#d95f02" )) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(x = NULL, y = "Percent of matching genes", fill = NULL) +
        theme_classic(base_size = 20) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

        ggsave(paste0(out_path, "/gene_var_overlap.png"),  match_plot, width = 10, height = 10)
}


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
      # add og clusters
      seurat_object <- AddMetaData(seurat_object, aligned, col.name = "seurat_clusters")
      lvl <- unique(as.character(cluster_map_raw$cluster[cluster_map_raw$age_sex == age_sex_group]))
      seurat_object$seurat_clusters <- factor(seurat_object$seurat_clusters, levels = lvl)
      Idents(seurat_object) <- "seurat_clusters"
      # subset if cells in new version dont match 
      seurat_object <- subset(seurat_object, cells = Cells(seurat_object)[!is.na(seurat_object$seurat_clusters)])
      seurat_object$seurat_clusters <- droplevels(seurat_object$seurat_clusters)
     
    

      return(seurat_object)
    }

    resave_with_bp_cells <- function(seurat_object, output_dir) {

        object_name <- "Spleen"

        seurat_object[["originalexp"]] <- as(seurat_object[["originalexp"]], Class = "Assay5")

        # re save in new format
        matrix_dir = file.path(output_dir, "data", paste0(object_name, "_counts"))
        if (dir.exists(matrix_dir)) {unlink(matrix_dir, recursive = TRUE)} # will throw an error otherwise

        write_matrix_dir(mat = seurat_object[["originalexp"]]$counts, dir = matrix_dir)
        counts.mat <- open_matrix_dir(dir = file.path(output_dir, "data", paste0(object_name, "_counts")))
        seurat_object[["originalexp"]]$counts <- counts.mat
        rm(counts.mat)
        gc()
        print("new matrix is done")
        return(seurat_object)
    }

