

# data prep

# used data prep version 2 and now use it for immune cells


# this is the next version of the script where certain steps are improved
## do sum and not mean for promoter motivs that are there multiple times
## age separately




assign_variability_ages <- function(){
  # in different clusters the gene can be differently variable: find uniform way
  marrow_spleen_seurat <- readRDS("data_prep_for_classification/immune_models/data/granulocytes_all_ages.rds")
  hep <-  marrow_spleen_seurat@meta.data
  
# ok this way of processing does not work for seurat, and i think i have too many cells to convert it? 
#try which one, also how do i calc gmean? 

  hep_master_df <- hep %>%                            
      mutate(age_sex = age,
              age = case_when(
              str_detect(age_sex, "18m") ~ "18m",
              str_detect(age_sex, "24m") ~ "24m",
              str_detect(age_sex, "3m")  ~ "3m"),
          cluster_id = paste0(tissue,"_" ,age_sex, "_", cluster)) %>%
      select(gene, cluster_id, age, res_var, gmean, variability) %>%
      group_by(gene) %>%                     
      mutate( 
        var_non_amb = case_when(
      any(variability == "HVG") ~ "HVG",
      mean(variability == "LVG") >= 0.75 ~ "LVG",
      TRUE ~ "Inconclusive")) 
  return(hep_master_df)
}
hep_non_ambiguous_df <- assign_variability_ages()
hep_non_ambiguous_df


merge_hep_promoters <- function(hepatocytes_master_df){
  pred <- read_tsv( "data_prep_for_classification/hepatocyte_models/results/predictors_non_amb.tsv")
  
  hep_promoters_merged_df <- hepatocytes_master_df %>%
    left_join(pred, by = "gene") %>%
    filter(!is.na(TATA_box))     
  return(hep_promoters_merged_df)       
}
hep_and_prom_df <- merge_hep_promoters(hep_non_ambiguous_df)
table(hep_and_prom_df$var_non_amb)




gmean_matching_and_subsample <- function(hep_and_promoters_df){

  df <- hep_and_promoters_df %>%
    filter(var_non_amb %in% c("HVG", "LVG")) %>%
    group_by(gene, age, var_non_amb) %>%
    mutate(median_gmean = median(gmean)) %>%
    ungroup() %>%
    distinct(gene, age, var_non_amb, .keep_all = TRUE)


  set.seed(1)

    df_age <- df 
    hvg <- df_age %>% filter(var_non_amb == "HVG")
    lvg <- df_age %>% filter(var_non_amb == "LVG")

    n_matches <- min(nrow(hvg), nrow(lvg))
    hvg <- hvg %>% slice_sample(n = n_matches) # random sampling without replacement according to setseed
    lvg_pool <- lvg   
    lvg_matches <- vector("list", n_matches)

    for (i in seq_len(n_matches)) {
      val <- hvg$median_gmean[i]
      lvg_candidate <- lvg_pool %>%
        mutate(diff = abs(median_gmean - val)) %>%
        slice_min(order_by = diff, n = 1, with_ties = TRUE) %>%
        sample_n(1)

      lvg_matches[[i]] <- lvg_candidate
      lvg_pool <- lvg_pool %>% filter(gene != lvg_candidate$gene)
    }
    lvg_matched <- bind_rows(lvg_matches)

    # test bias to certian positions in lvg pool
    match_list <- match(lvg_matched$gene, lvg$gene)
    match_clean <- match_list[!is.na(match_list)]
    match_df <- data.frame(position = match_clean)
    plot <- ggplot(match_df, aes(x = position)) +
        geom_histogram(binwidth = 50, fill = "steelblue", colour = "white") +
        labs(x = "Position in pool of possible LVG matches", y = "Count", title = paste("LVG Match Position Distribution for all ages")) +
        theme_classic()
    ggsave(paste0("data_prep_for_classification/hepatocyte_models/results/sampling_distribution_all.png"), plot = plot)

    matched_age <- bind_rows(hvg[1:nrow(lvg_matched), ], lvg_matched) %>%
      mutate(var_non_amb = factor(var_non_amb, levels = c("HVG", "LVG")), age = age) 



  return(matched_age)
}
      
matched_df  <- gmean_matching_and_subsample(hep_and_prom_df)
table(matched_df$age)




gmean_matching_and_subsample_age <- function(hep_and_promoters_df){

  df <- hep_and_promoters_df %>%
    filter(var_non_amb %in% c("HVG", "LVG")) %>%
    group_by(gene, age, var_non_amb) %>%
    mutate(median_gmean = median(gmean)) %>%
    ungroup() %>%
    distinct(gene, age, var_non_amb, .keep_all = TRUE)

  matched_all_ages <- list()
  set.seed(1)

  for (age_group in unique(df$age)) {

    df_age <- df %>% filter(age == age_group)
    hvg <- df_age %>% filter(var_non_amb == "HVG")
    lvg <- df_age %>% filter(var_non_amb == "LVG")

    n_matches <- min(nrow(hvg), nrow(lvg))
    hvg <- hvg %>% slice_sample(n = n_matches) # random sampling without replacement according to setseed
    lvg_pool <- lvg   
    lvg_matches <- vector("list", n_matches)

    for (i in seq_len(n_matches)) {
      val <- hvg$median_gmean[i]
      lvg_candidate <- lvg_pool %>%
        mutate(diff = abs(median_gmean - val)) %>%
        slice_min(order_by = diff, n = 1, with_ties = TRUE) %>%
        sample_n(1)

      lvg_matches[[i]] <- lvg_candidate
      lvg_pool <- lvg_pool %>% filter(gene != lvg_candidate$gene)
    }
    lvg_matched <- bind_rows(lvg_matches)

    # test bias to certian positions in lvg pool
    match_list <- match(lvg_matched$gene, lvg$gene)
    match_clean <- match_list[!is.na(match_list)]
    match_df <- data.frame(position = match_clean)
    plot <- ggplot(match_df, aes(x = position)) +
        geom_histogram(binwidth = 50, fill = "steelblue", colour = "white") +
        labs(x = "Position in pool of possible LVG matches", y = "Count", title = paste("LVG Match Position Distribution for Age", age_group)) +
        theme_classic()
    ggsave(paste0("data_prep_for_classification/hepatocyte_models/results/sampling_distribution_", age_group, ".png"), plot = plot)

    matched_age <- bind_rows(hvg[1:nrow(lvg_matched), ], lvg_matched) %>%
      mutate(var_non_amb = factor(var_non_amb, levels = c("HVG", "LVG")), age = age_group) 
    matched_all_ages[[age_group]] <- matched_age
  }
  final_matched <- bind_rows(matched_all_ages)

  return(final_matched)
}
      
matched_df  <- gmean_matching_and_subsample_age(hep_and_prom_df)
table(matched_df$age)


save_matched_df <- function(matched_df){
    write_tsv(matched_df, file = 'data_prep_for_classification/hepatocyte_models/results/matched_df.tsv')
}   
save_matched_df(matched_df)



# add the TF binding predicitons
add_tf_binding <- function(matched_df){
  # here i check the new data with the TF binding predictions
 
  familial_profiles_list <- readRDS("data_prep_for_classification/familial_binding_profile_predictions/jaspar/tf_profiles_list.rds")
  selected_profile_df <- familial_profiles_list[['_250_50']]

  tf_binding_df <- matched_df |>
    left_join(selected_profile_df, by = 'gene') 
    
  saveRDS(tf_binding_df, 'data_prep_for_classification/hepatocyte_models/data/matched_tf_binding_df.rds')  
  return(tf_binding_df)
}




# QC

check_match_pairs <- function(matched_df){

    for (age_group in unique(matched_df$age)) {

        df_age <- matched_df %>% filter(age == age_group)
        # see how well the median gmeans of HVG and LVG match 
        hvg_df <- df_age %>% filter(var_non_amb == "HVG") 
        lvg_df <- df_age %>% filter(var_non_amb == "LVG") 

        hist <- ggplot(df_age, aes(x = median_gmean, fill = var_non_amb)) +
            geom_histogram(position = "identity", alpha = 0.6, binwidth = 0.2, colour = "white") +
            scale_fill_manual(values = c("HVG" = "darkred", "LVG" = "orange")) +
            labs(
            title = "Median gmeans of matched HVG/LVG pairs",
            x = "Median gmean",
            y = "Number of genes",
            fill = "Variability class" ) +
            theme_classic()

        ggsave(filename = paste0("data_prep_for_classification/hepatocyte_models/results/hist_hvg_lvg_", age_group, ".png"),plot = hist)
    }
}
match <-  check_match_pairs(matched_df)


check_random_pairs <- function(hep_and_prom_df){
  # see how a random sampling makes the median gmeans of HVG and LVG match 
    
    gmean_df <- hep_and_prom_df %>%
        filter(var_non_amb %in% c("HVG", "LVG")) %>%
        group_by(gene, age, var_non_amb) %>%
        mutate(median_gmean = median(gmean)) %>%
        ungroup() %>%
        distinct(gene, age, var_non_amb, .keep_all = TRUE)


    for (age_group in unique(gmean_df$age)) {

        df_age <- gmean_df %>% filter(age == age_group)

        hvg_df <- df_age %>% filter(var_non_amb == "HVG") 
        lvg_df <- df_age %>% filter(var_non_amb == "LVG") 

        set.seed(1)
        n_matches <- min(nrow(hvg_df), nrow(lvg_df))
        hvg_pool <- hvg_df %>% slice_sample(n = n_matches)
        lvg_pool <- lvg_df %>% slice_sample(n = n_matches)
        fused_df <- bind_rows(hvg_pool, lvg_pool)
        dim(fused_df)

        hist <- ggplot(fused_df, aes(x = median_gmean, fill = var_non_amb)) +
            geom_histogram(position = "identity", alpha = 0.6, binwidth = 0.2, colour = "white") +
            scale_fill_manual(values = c("HVG" = "darkred", "LVG" = "orange")) +
            labs(
            title = "Median gmeans of random HVG/LVG pairs",
            x = "Median gmean",
            y = "Number of genes",
            fill = "Variability class" ) +
            theme_classic()

        ggsave(filename = paste0("data_prep_for_classification/hepatocyte_models/results/hist_random", age_group, ".png"), plot = hist)
    }
}
random <- check_random_pairs(hep_and_prom_df)



check_diff <- function(matched_df){

    for (age_group in unique(matched_df$age)) {

        df_age <- matched_df %>% filter(age == age_group)

        hist <- ggplot(df_age, aes(x = diff)) +
            geom_histogram(position = "identity", alpha = 0.6, binwidth = 0.2, colour = "white") +
            scale_fill_manual(values = c("HVG" = "darkred", "LVG" = "orange")) +
            labs(
            title = "Difference of HVG and LVG in median gmean",
            x = "Difference in median gmean",
            y = "Number of genes",
            fill = "Variability class" ) +
            theme_classic()

        ggsave(filename = paste0("data_prep_for_classification/hepatocyte_models/results/hist_diff_", age_group, ".png"), plot = hist)
    }
}
match <-  check_diff(matched_df)
