

gene_distribution <- function(df, data_source) {
  # for a selection of genes show their distribution across categories

  gene_set_order <- names(df)[grepl(" gene$", names(df))]

  all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
   "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
   "Not expressed 0")
  # select genes
  #selected_genes_list <- c("Il6", "Tnf", "Il12a") # M1
  # selected_genes_list <- c("Il10", "Tgfb1") # M2
  # selected_genes_list <- c("Cox6b1", "Ndufc1 ", "Ndufb7") # oxphos
  # selected_genes_list <- c("Pkm", "Pfkfb3", "Slc2a1") # glycolysis
  selected_genes_list <- c("Il6", "Tnf", "Il12a", "Il10", "Tgfb2", "Cox6b1", "Ndufc1", "Ndufb7", "Pkm2", "Pfkfb3", "Slc2a1", "Gapdh") # all
  selected_genes_list <- c( "Il1b" , "Gapdh")
  
 # subearch genes
 #df_matches <- df %>%
  #filter(str_detect(gene, "Gapdh")) 

  counts_df <- df |>
    filter(gene %in% selected_genes_list) |>
    pivot_longer(
      cols      = all_of(gene_set_order),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    mutate(
      gene_set    = factor(gene_set, levels = gene_set_order),
      category    = factor(category, levels = all_cats),
      cluster_plot = paste(cell_type, "\n", cluster_name),
      gene        = factor(gene, levels = selected_genes_list)) |>
    distinct(gene, cluster_plot, .keep_all = TRUE) |>
    group_by(gene, category) |>
    count(name = "n_times") |>
    ungroup() |>
    complete(gene, category, fill = list(n_times = 0))


  # color scheme
  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)

  plot <- ggplot(counts_df, aes(category, n_times, fill = category)) +
    geom_col() +
    facet_wrap(~ gene, scales = "free_y") +
    guides(fill = guide_legend(nrow = 6)) +
    scale_y_continuous(limits = c(0, 9), breaks = 1:9 ) +
    labs(
      x = "Expression category",
      y = "Cluster count",
      fill = "Gene group") +
    scale_fill_manual(values = col_vec) +
    theme_classic() +
    guides(colour = guide_legend(ncol = 1, nrow = 6)) +
    theme(
      strip.placement       = "outside",
      strip.text.y.left     = element_text(angle = 0),
      strip.text.y.right    = element_blank(),
      strip.background      = element_blank(),
      axis.text.x           = element_text(size = 20, angle = 45, hjust = 1),
      axis.text.y           = element_text(size = 20),
      axis.title.x           = element_text(size = 20),
      axis.title.y           = element_text(size = 20),
      legend.position       = "top",
      legend.title = element_text(size = 20),   
      legend.text  = element_text(size = 20),  
      strip.text   = element_text(size = 20))

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_12/gene_distribution.png", plot, width = 20, height = 8)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_12/gene_distribution.png", plot, width = 20, height = 8)
    } else {
      stop("Issue when saving.")
    }   

    return(counts_df)
}




dummy_entropy <- function(data_source) {
  # create dummy data to test raos entropy

  categories <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
  "Not expressed 0")
  genes <- c("same_in_all_cats_big","same_in_all_cats_small", "just_in_lvg", "just_in_hvg", "just_in_interm", "just_in_not", "lvg_narrow", "lvg_wide", "sparse_distributed")

  df_dummy <- 
    full_grid <- expand_grid(
      gene = genes,
      category = categories,
      n_times = 0)

  # same amount of genes
  df_dummy$n_times[df_dummy$gene == "same_in_all_cats_small"] <- 1
  df_dummy$n_times[df_dummy$gene == "same_in_all_cats_big"] <- 2
  df_dummy$n_times[df_dummy$gene == "just_in_lvg" & grepl("^LVG", df_dummy$category)] <- 6
  df_dummy$n_times[df_dummy$gene == "just_in_hvg" & grepl("^HVG", df_dummy$category)] <- 6
  df_dummy$n_times[df_dummy$gene == "just_in_interm" & grepl("^Intermediate", df_dummy$category)] <- 6
  df_dummy$n_times[df_dummy$gene == "just_in_not" & grepl("^Not", df_dummy$category)] <- 18

  # different amount of genes
  lvg_values_narrow <- c( "LVG 1" = 4, "LVG 2" = 4, "LVG 3" = 4,"LVG 4" = 0, "LVG 5" = 0, "LVG 6" = 0)
  lvg_values_wide <- c( "LVG 1" = 6, "LVG 2" = 0, "LVG 3" = 6, "LVG 4" = 0, "LVG 5" = 0, "LVG 6" = 6)
  sparse_dist <- c( "LVG 1" = 5, "HVG 1" = 6, "Intermediate 1" = 7)

  for (cat in names(lvg_values_narrow)) {
     df_dummy$n_times[df_dummy$gene == "lvg_narrow" & df_dummy$category == cat] <- lvg_values_narrow[cat]}

  for (cat in names(lvg_values_wide)) {
     df_dummy$n_times[df_dummy$gene == "lvg_wide" & df_dummy$category == cat] <- lvg_values_wide[cat]}
 

  for (cat in names(sparse_dist)) {
     df_dummy$n_times[df_dummy$gene == "sparse_distributed" & df_dummy$category == cat] <- sparse_dist[cat]}

  # color scheme
  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = categories))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)

  plot <- ggplot(df_dummy, aes(category, n_times, fill = category)) +
    geom_col() +
    facet_wrap(~ gene, scales = "free_y") +
    labs(
      x = "Expression category",
      y = "Gene count") +
    scale_fill_manual(values = col_vec) +
    theme_classic() +
    theme(
      strip.placement       = "outside",
      strip.text.y.left     = element_text(angle = 0),
      strip.text.y.right    = element_blank(),
      strip.background      = element_blank(),
      axis.text.x           = element_text(angle = 45, hjust = 1),
      legend.position       = "top")

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_12/distribution_dummy.png", plot, width = 12, height = 5)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_12/distribution_dummy.png", plot, width = 12, height = 5)
    } else {
      stop("Issue when saving.")
    }   
  return(df_dummy)
}



subsample_genes_balanced <- function(df) {
  hvg_df <- strat_df |>
    select(gene, gmean, cluster_id, tissue, res_var, gene_variability, expression_bin, category ) |>
    filter(gene_variability != "Not expressed") |>
    filter(gene_variability == "HVG") |>
    group_by(cluster_id, expression_bin) |>
    summarise(nr_hvgs = n(), .groups = "drop")
  head(hvg_df)
  hvg_df$nr_hvgs


  annotated_df <- strat_df %>%
  select(gene, gmean, cluster_id, tissue, res_var, gene_variability, expression_bin, category ) |>
  filter(gene_variability != "Not expressed") |>
  left_join(hvg_df, by = c("cluster_id", "expression_bin")) |>
    mutate(nr_hvgs = coalesce(nr_hvgs, 0L)) |>
  select(gene, cluster_id, gene_variability, category, nr_hvgs, gene) |>
  arrange(cluster_id, nr_hvgs)
  tail(annotated_df)

 all_categories  <- c(
  "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
  "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6")

  out <- list()
  set.seed(789)
  for (cat in all_categories) {
    dcat <- annotated_df %>% filter(category == cat)
    samp <- dcat |>
      group_by(cluster_id, category)  %>%               
      sample_n(size = min(first(nr_hvgs), n())) %>%  
      ungroup()
  out[[cat]] <- samp 
  }
  sampled_all <- bind_rows(out)

  check <- sampled_all |>
    group_by(gene_variability, category) |>
    summarise(count = n())
  check

sampled_all
}



all_genes_entropy_prep <- function(category_df, data_source) {
  # prepare the staratified data for the entropy calculation

# for doc vs 2
  all_categories  <- c(
  "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
  "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
  "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6")

  n_clusters_df <- sampled_all |>
    filter(category != "Not expressed 0") |> # only the ones from the gene sets are removed not the ones that are were in cluster
    group_by(gene) |>
    summarise(n_clusters = n_distinct(cluster_id), .groups = "drop") |>
    filter(n_clusters > 5) 
  print(n_clusters_df, n = 300)

  gene_category_counts <- strat_df |>
    filter(category != "Not expressed 0") |>
    semi_join(n_clusters_df, by = "gene") |>          
    count(gene, category, name = "n_category") |>      
    complete(gene, category = all_categories,           
            fill = list(n_category = 0)) |>
    left_join(n_clusters_df, by = "gene") |>        
    mutate(n_times = if_else(n_clusters > 0, n_category / n_clusters, 0)) |>
    arrange(gene, category)

  print(gene_category_counts, n = 300)

  return(gene_category_counts)
}




subsample_genes <- function(df) {
  # subsample the genes (for fast entropy calculation)
  all_genes_rao_prep <- gene_category_counts
  all_genes <- unique(all_genes_rao_prep$gene)
  set.seed(4)  
  n   <- length(all_genes)
  k   <- min(3000, n) 
  subsampled_genes <- sample(all_genes, k, replace = FALSE)
  subsampled_df <- all_genes_rao_prep |>
    filter(gene %in% subsampled_genes)
print(subsampled_df, n = 200)

  return(subsampled_df)
}



quadratic_entropy <- function(df, data_source, write = FALSE) {
  # raos q entropy measures diversity while taking the similarities of the categories into account
  df <- gene_category_counts
  data_source <- "droplet"
  gene_list <- unique(df$gene)
  all_categories <- unique(df$category) 

  # distance matrix 
  #D <- readRDS(file = "facs/data/dissimilarity_matrix_rao.rds")
  rownames(D) <-  all_categories
  colnames(D) <- all_categories

  gene_similarity_df <- data.frame(gene = character(), rao_q = numeric(), stringsAsFactors = FALSE)

  for (gene in gene_list) {
    df_gene <- df |> filter(gene == !!gene)

    p <- setNames(rep(0, length(all_categories)), all_categories)
    p[df_gene$category] <- df_gene$n_times 

    rao_q <- sum(outer(p, p) * D)

    gene_similarity_df <- rbind(gene_similarity_df, data.frame(gene = gene, rao_q = rao_q))
  }
  if (write) {
    saveRDS(gene_similarity_df, paste0(data_source, '/data/df_entropy_full.rds'))
    print("Saved df to data.")
  }
  return(gene_similarity_df)
}






plot_raos_gene_selection <- function(df, data_source) {
  # plot the different entropies for a selection of genes in a bar plot

  plot <- ggplot(df, aes(gene, rao_q, fill = gene)) +
    geom_col(fill = 'grey') +
    labs(
      x = "Genes",
      y = "Rao's quadratic entropy") +
    theme_classic() +
    coord_flip() +
    theme(
      strip.placement       = "outside",
      strip.text.y.left     = element_text(angle = 0),
      strip.text.y.right    = element_blank(),
      strip.background      = element_blank(),
      axis.text.x           = element_text(angle = 45, hjust = 1),
      legend.position       = "top")

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_12/distribution_similarity.png", plot, width = 3, height = 6)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_12/distribution_similarity.png", plot, width = 3, height = 6)
    } else {
      stop("Issue when saving.")
    }   
    return(df)
}



plot_raos_all_genes <- function(df, data_source) {
  # plot entropy for a lot of genes in a histogram
  
  df <- readRDS("droplet/data/df_entropy_full.rds")

  plot <- ggplot(df, aes(x = rao_q)) +
    geom_histogram(binwidth = 0.05, fill = "#A2C759") + #A2C759 #A9D18E
   # scale_x_continuous(
    #  name = "Rao’s quadratic entropy",
    #  breaks = seq(0, max(df$rao_q, na.rm = TRUE), by = 0.5),
     # labels = scales::number_format(accuracy = 0.5)) +
    labs(x = "Rao’s quadratic entropy",
         y = "Count of genes") +
  theme_classic() +
  theme(
      strip.background      = element_blank(),
      axis.text.x           = element_text(size = 20),
      axis.text.y           = element_text(size = 20),
      axis.title.x      = element_text(size = 20),
      axis.title.y      = element_text(size = 20))

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_12/hist_entropy.png", plot, width = 12, height = 5)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_12/hist_entropy.png", plot, width = 12, height = 5)
    } else {
      stop("Issue when saving.")
    }   

    return(df)
}



merge_entropy_results <-  function(data_source) {
  # fuse the entropy and the master table
  df_entropy <- readRDS(paste0(data_source, "/data/df_entropy_full.rds")) 
  category_df <- readRDS(paste0(data_source, "/data/strat_df.rds"))

  master_entropy_df <- df_entropy |>
    left_join(category_df, by = 'gene') 
  saveRDS(master_entropy_df, paste0(data_source,"/data/rao_q_metadata.rds"))

  selection_for_checking <- master_entropy_df |>
    select(gene, category, gmean, gene_variability, rao_q, cluster_id) |>
    arrange(rao_q)
  write.csv(selection_for_checking, paste0(data_source, "/data/df_selected.csv"))

  return(master_entropy_df)
}




calc_variability_direction <-  function(entropy_master_df, data_source) {
  
  #entropy_metadata_df <- readRDS("facs/data/rao_q_metadata.rds")

  subtraction_df <- entropy_master_df |>
   # filter(!gene_variability %in% c("Intermediate", "Not expressed")) |>
    select(gene, cluster_id, gene_variability) |>
    unique()  |>
    mutate(var_multiplied = case_when(gene_variability == "HVG" ~ 1,
                                      gene_variability == "LVG" ~ (-1), 
                                      gene_variability == "Intermediate" ~ 0,
                                      gene_variability == "Not expressed" ~ 0)) |>   
    group_by(gene) |>                   
    summarise(variability_direction = sum(var_multiplied))

  # merge back to master df
   merged_subtraction_df <- master_entropy_df |>
    #filter(!gene_variability %in% c("Intermediate", "Not expressed")) |>
   #  select(gene, rao_q) |>
    # unique()  |>
    left_join(subtraction_df, by = "gene") 
  saveRDS(merged_subtraction_df, paste0(data_source,"/data/variability_dir_metadata.rds"))

  plot_subtraction_score_df <- merged_subtraction_df |>
    select(gene, variability_direction, rao_q) |>
    unique()


plot <- ggplot(plot_subtraction_score_df, aes(x = variability_direction, y = rao_q)) +
  geom_jitter(width = 0.2, height = 0, size = 1, alpha = 0.3, color = "darkblue") +
  labs(
    x = "Variability direction",
    y = "Rao's quadratic entropy") +
  theme_classic()


  if (data_source == "facs") {
    print("Saved to: FACS.")
    ggsave("facs/plots/3_m/Test_14/subtraction_score.png", plot, width = 4, height = 4)
  } else if (data_source == "droplet") {
    print("Saved to: droplet")
    ggsave("droplet/plots/3_m/Test_14/subtraction_score.png", plot, width = 4, height = 4)
  } else {
    stop("Issue when saving.")
  } 
  return(merged_subtraction_df)
}




scatter_entropy_hvg <-  function(variability_df, data_source) {
  # are there any patterns with gene set and entropy
  only_hvg_df <- variability_df |>
    filter(gene_variability == "HVG") |>
    filter(perc_hvg > 0.9) |>
    select(gene,  perc_hvg, rao_q, gene_variability) |>
    unique() |>
   # filter(!is.na(perc_hvg))
    arrange((perc_hvg))
  head(only_hvg_df)


plot <- ggplot(only_hvg_df, aes(x = perc_hvg, y = rao_q)) +
  geom_jitter(width = 0.001, height = 0, size = 1, alpha = 0.5, color = "orange") +
  labs(
    x = "HVG bootstrapping percentage",
    y = "Rao's entropy") +
  #  xlim(0, max(only_hvg_df$perc_hvg, na.rm = TRUE)) +
  theme_classic()

  if (data_source == "facs") {
    print("Saved to: FACS.")
    ggsave("facs/plots/3_m/Test_14/scatter_entropy.png", plot, width = 7, height = 7)
  } else if (data_source == "droplet") {
    print("Saved to: droplet")
    ggsave("droplet/plots/3_m/Test_14/scatter_entropy.png", plot, width = 7, height = 7)
  } else {
    stop("Issue when saving.")
  } 
}



# genes with low entropy -> which gene sets? 

low_entropy_gene_sets <- function(data_source){
 data_source <- "droplet"
 rao_df <- readRDS(paste0(data_source, "/data/rao_q_metadata.rds"))
rao_df

 all_cats <- c(
   "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
   "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
   "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
   "Not expressed 0")


  counts_df <- rao_df |>
   # mutate(`Control gene` = `Meiosis gene` | 
                         #  `Sperm DNA condensation gene` | 
                          # `Oocyte maturation gene`) |>                 
   # select(c(-`Meiosis gene`, - `Sperm DNA condensation gene`, - `Oocyte maturation gene`)) |>
    pivot_longer(
      cols      = matches(" gene$"),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    # filter out gene sets i dont need
    filter(!gene_set %in% c("Housekeeping Lin gene", "LPS response gene", "Other gene", "Autoimmunity gene", "Cytokine response gene")) |>
    mutate(
      #gene_set    = factor(gene_set, levels = gene_set_order),
      category    = factor(category, levels = all_cats),
      cluster_plot = paste(cluster_name,  "| ", cell_type ),
      rao_q_low = rao_q > 30) |>
    #count(gene_set, cluster_plot, category, name = "n_genes") |>
    filter(gene_set != "Not expressed") |>
    group_by(gene_set) |>
    summarise(
      total_genes    = n(),
      low_genes      = sum(rao_q_low, na.rm = TRUE),
      fraction_low_rao   = low_genes / total_genes,
      .groups = "drop" )
print(counts_df, n = 30)


  plot <- ggplot(counts_df, aes(fraction_low_rao, gene_set)) +
    geom_col(width = 0.8, fill = "#40CAF0") + #F76C7E #85A5F1 #42D4D5 #40CAF0 #EFDC5F
    guides(fill = guide_legend(nrow = 6)) +
    labs(
      y = "Gene sets",
      x = "Fraction of genes with entropy < 3") +
    theme_classic() +
    theme(
      strip.placement       = "outside",
      strip.text.y.left     = element_text(angle = 0),
      strip.text.y.right    = element_blank(),
      strip.background      = element_blank(),
      axis.text.x           = element_text(angle = 45, hjust = 1, size = 30),
      axis.text.y           = element_text(size = 30),
      axis.title.x      = element_text(size = 30),   
      axis.title.y      = element_blank(),        
      legend.position       = "top") #+
      #geom_signif( comparisons = list(c("TNF response gene", "Cytokine response gene")),
                 #  map_signif_level = TRUE, annotations =  "helene is cool")

#annotations= "***"
 ggsave(paste0(data_source, "/plots/3_m/Test_11/low_entropy_gene_gets.png"), plot, width = 12, height = 7)

}










# there are NA raoq genes, 
# not all genes are in the same number of clusters

rest <- function() {
  t <- strat_df |>
    select(gene, cluster_id) |>
    distinct() |>
    group_by(gene) |>
    summarise(num_clusters = n()) |>
    ungroup() |>
    arrange(desc(num_clusters)) |>
    filter(num_clusters > 11)
  print(t, n = 100)  

  # 18 089 genes'
  # 9 991 = 14
  # 11 400 > 12
  # 12 400 > 11

  D <- matrix(c(
      0,1,2,3,4,5,10,10,10,10,10,10,10,10,10,10,10,10,
      1,0,1,2,3,4,10,10,10,10,10,10,10,10,10,10,10,10,
      2,1,0,1,2,3,10,10,10,10,10,10,10,10,10,10,10,10,
      3,2,1,0,1,2,10,10,10,10,10,10,10,10,10,10,10,10,
      4,3,2,1,0,1,10,10,10,10,10,10,10,10,10,10,10,10,
      5,4,3,2,1,0,10,10,10,10,10,10,10,10,10,10,10,10,
      10,10,10,10,10,10,0,1,2,3,4,5,10,10,10,10,10,10,
      10,10,10,10,10,10,1,0,1,2,3,4,10,10,10,10,10,10,
      10,10,10,10,10,10,2,1,0,1,2,3,10,10,10,10,10,10,
      10,10,10,10,10,10,3,2,1,0,1,2,10,10,10,10,10,10,
      10,10,10,10,10,10,4,3,2,1,0,1,10,10,10,10,10,10,
      10,10,10,10,10,10,5,4,3,2,1,0,10,10,10,10,10,10,
      10,10,10,10,10,10,10,10,10,10,10,10,0,1,2,3,4,5,
      10,10,10,10,10,10,10,10,10,10,10,10,1,0,1,2,3,4,
      10,10,10,10,10,10,10,10,10,10,10,10,2,1,0,1,2,3,
      10,10,10,10,10,10,10,10,10,10,10,10,3,2,1,0,1,2,
      10,10,10,10,10,10,10,10,10,10,10,10,4,3,2,1,0,1,
      10,10,10,10,10,10,10,10,10,10,10,10,5,4,3,2,1,0

    ), nrow = 18, byrow = TRUE)
}