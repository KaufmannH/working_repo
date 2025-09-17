
# this script makes intersects of the profiles of each TF binding prediction


familial_profiles_list <- list()

for (up_down in c("_30_30", "_50_50", "_250_50")) {
  promoters_bed <- read_tsv(paste0("data_prep_for_classification/familial_binding_profile_predictions/jaspar/promoters", up_down, ".bed"), col_names = F)
  colnames(promoters_bed) <- c("chrom", "start", "end", "promoter", "score", "strand")
  promoters_bed$index <- seq(1:nrow(promoters_bed))
  
  out_intersection <- read_tsv(paste0("data_prep_for_classification/familial_binding_profile_predictions/jaspar/out_intersected", up_down, ".bed"), col_names = F)
  colnames(out_intersection) <- c("chrom", "start", "end", "profile", "score", "strand", "start2", "end2", "rgb", "index")
  
  out_intersection_pivot <- out_intersection %>%
    mutate(profile = str_replace_all(profile, " ", "_")) %>%
    select(index, profile, strand) %>% 
    group_by(index, profile) %>% 
    summarize(n_occur = n()) %>%
    mutate(occur = 1) %>% # set 1 if occurs (don't care how many predicted binding sites)
    select(-n_occur) %>%
    pivot_wider(names_from = profile, values_from = occur, values_fill = 0) # fill non-occurences with 0
  
  familial_profiles <- left_join(promoters_bed, out_intersection_pivot, by = c("index" = "index"))
  
  gene_description <- read_tsv("data_prep_for_classification/general_models/epd_mm/gene_description.txt", col_names = c("promoter", "gene", "description"))
  familial_profiles <- familial_profiles %>% left_join(gene_description, by = c("promoter" = "promoter")) %>%
    relocate(gene, .after = promoter) %>%
    select(-c("chrom", "start", "end", "description", "score"))

  familial_profiles_list[[up_down]] <- familial_profiles %>%
    group_by(gene) %>%
    summarize(across(starts_with("Familial_profile"), mean))
}

saveRDS(familial_profiles_list, "data_prep_for_classification/familial_binding_profile_predictions/jaspar/tf_profiles_list.rds")


# H match the familial profiles with the names
familial_profiles_list <- readRDS("data_prep_for_classification/familial_binding_profile_predictions/jaspar/tf_profiles_list.rds")
familial_profile_metadata <- readr::read_tsv('data_prep_for_classification/familial_binding_profile_predictions/jaspar/clusters.tab.txt')

for (frame_name in names(familial_profiles_list)) {
  frame <- familial_profiles_list[[frame_name]]

  cleaned_df <- frame |>
    pivot_longer(
      cols = -gene,
      names_to = "familial_profile",
      values_to = "value") |>
      mutate(familial_profile_match = str_replace(
          familial_profile,  "Familial_profile_(\\d+)", 
          ~ sprintf("cluster_%03d", as.integer(str_extract(.x, "\\d+"))))) |>
    left_join(familial_profile_metadata, by = c("familial_profile_match" = "cluster")) |>
    rename(tf_names = name) |>
    select(gene, familial_profile, tf_names) 

    saveRDS(cleaned_df, paste0("data_prep_for_classification/familial_binding_profile_predictions/jaspar/tf_names_df", frame_name, ".rds"))
}


