# this is analogous to assemble_TMS_df.R which brings FACS and droplet into the same shape for the stratification pipeline


out_path <- '/home/hkaufm49/analyses/working_repo/PanSci/data/prepped_for_strat/'



# in the TMS the results were saved in a tibble but not for pansci, so i added the code here

quantile_perc <- c(.025, .975) # this was specified in the pipeline_server script
tissue <- 'Liver'
bs_path <- "/home/hkaufm49/analyses/RegulatedNoiseHel/data/processed/pansci_liver/bootstrapping_tables/"
bs_files <- list.files(bs_path, pattern = "\\.rds$", full.names = TRUE)


for (f in bs_files) {

 # f <- bs_files[1]

  message("Reading: ", f)
  age_group <-  str_extract(f, "(?<=_)[0-9]{2}_months_[A-Za-z]+(?=\\.rds)")
  bootstrap_res <- readRDS(f)
  res_var_cl <- readRDS(paste0('/home/hkaufm49/analyses/RegulatedNoiseHel/data/processed/pansci_liver/res_var_tables/res_var_cl_Liver_', age_group, '.rds' ))

  stats <- bootstrap_res %>% 
    group_by(gene, cluster) %>% 
    summarize(mean_resvar_bs = mean(ResVar),
              perc.hvg = sum(boolean_hvg) / n(),
              quant_low = quantile(ResVar, probs = quantile_perc[1]), 
              quant_high = quantile(ResVar, probs = quantile_perc[2]))
  
  results_tmp_tibble <- left_join(res_var_cl, stats, by = c("gene" = "gene", "cluster" = "cluster")) %>%
    mutate(tissue = tissue, age = age_group)

  if (f == bs_files[1]) {
    results_tibble <- results_tmp_tibble
  } else {
    results_tibble <- rbind(results_tibble, results_tmp_tibble)
  }

}

results_tibble[1:10, 4:11]
#saveRDS(results_tibble, file = file.path(out_folder, "results_tibble.rds"))


# filter for hepatocyte clusters


# stopped here: 
#get a df witht he age, sex and cluster and cell type, do it in prepr objects for regnosie pipeline
# the cluster in the resvar table is the seurat cluster
# i think i need to generate an id for the main cell types myself and then relate




assemble_easysci_df <- function(write = FALSE){

  # load data
  print("Loading data ...")
  gene_expression_data <- read.delim("droplet/droplet_raw_data/full_results_10_22.tsv", header = TRUE, sep = "\t")
  annotation_info_raw <- read_excel("droplet/droplet_raw_data/manual_annotation.xlsx")
  metadata_cell_numbers_raw <- read.table('droplet/droplet_raw_data/combined_metadata_for_cell_numbers.txt', header = TRUE, sep = "\t") # there are cluster_ids that are not in the annotaiton info
  print("Data loaded.")


  # prep dfs
  gene_expression_data <- gene_expression_data |>
    mutate(cluster_id = paste0(tissue , "_", age, "_", cluster)) 
  head(gene_expression_data)
  # every biological replicate and tissue have the same number of rows (genes)

  # prep dfs
  annotation_info <- annotation_info_raw  # cluster id exists already

  gene_expression_data <- gene_expression_data |>
    mutate(cluster_id = paste0(tissue , "_", age, "_", cluster))

  metadata_cell_numbers <- metadata_cell_numbers_raw |> # this df contains one row per cell
    mutate(cluster_id = paste0(tissue , "_", age_sex, "_", seurat_clusters)) |>
    group_by(cluster_id) |>
    summarise(num_cells_in_cluster = n(), .groups = 'drop') 


# join dfs
  numbers_joined <- gene_expression_data  |>
    left_join(metadata_cell_numbers, by = "cluster_id") |>
    select(-tissue) 
  colnames(numbers_joined)

  all_joined <- numbers_joined  |>
    left_join(annotation_info, by = c("cluster_id")) 
  colnames(all_joined)

  # tag LVGs
  df_incl_lvgs <- all_joined |>
    mutate(lvg = if_else(res_var < 1, TRUE, FALSE)) |>
    mutate(hvg = if_else(res_var >= 5 & gmean != 0, TRUE, FALSE))
  colnames(df_incl_lvgs)

  # select necessary columns
  df_finished <- df_incl_lvgs |>
    mutate(age = as.integer(str_extract(age, "\\d+(?=m)"))) |>
    rename(perc_hvg = perc.hvg,
           cell_type = manual_final) |>
           mutate(cell_type = str_remove(cell_type, "^(lymphocyte_|leukocyte_)")) |>
    select(gene, gmean, cluster_id, cell_type, tissue, age, res_var, perc_hvg, hvg, lvg)
  print("Df is done.")
 
  if (write) {
      # write combined df to a csv
    write.csv(df_finished, 'droplet/data/combined_data.csv')
    print("Saved df to data.")

    # get cell type names, write to excel for manual selection of innate immune cells
    cell_type_list <- data.frame(unique(df_finished$cell_type))
    write_xlsx(cell_type_list, 'droplet/data/cell_type_list.xlsx')
    print("Saved cell type list to data.")
  }

return(df_finished)
}


