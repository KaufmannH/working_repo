library(readxl)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr) 
library(purrr)
library(tibble)
library(writexl)

source("scripts/src/stratification.R")

# testing for the best way to get the LVG number down
# 1. change threshold
# 2. resample lvg numbers
# 3. deciding if threshold makes sense by enrichment plotting




# 1. check the resvar distibution of lvgs and test other thresholds


# droplet
strat_df <- readRDS("droplet/data/strat_df.rds") 
data_source <- "droplet"
# facs
facs <- readRDS("facs/data/strat_df.rds") 
data_source <- "facs"


lvg_test_df <- strat_df |>
    pivot_longer(
      cols      = ends_with(" gene"),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    select(gene, gmean, gene_set, res_var, category) 



zero_resvar <- lvg_test_df |>
    filter(res_var < 0.05) 
print(zero_resvar, n = 300)
# mainly genes with gmean 0 that were in tissue but not in cluster, house keeping and other genes

# see what the distribution looks like if i only look at expressed genes
no_expression_filtered <- strat_df |>
  filter(!gmean == 0)
no_expression_filtered


plot <- ggplot(strat_df, aes(x = res_var)) +
    geom_histogram(binwidth = 0.05, fill = "#A2C759") +
    labs(x = "Residual variance",
         y = "Count of genes") +
    coord_cartesian(xlim = c(0, 10)) + 
    geom_vline(aes(xintercept = 1), 
             colour = "black", linetype = "dashed", size = 1) +
  theme_classic() +
  theme(
      strip.background      = element_blank(),
      axis.text.x           = element_text(size = 20),
      axis.text.y           = element_text(size = 20),
      axis.title.x      = element_text(size = 20),
      axis.title.y      = element_text(size = 20))

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_16/res_var_distribution.png", plot, width = 12, height = 5)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_16/res_var_distribution.png", plot, width = 12, height = 5)
    } else {
      stop("Issue when saving.")
    }   


lvg_low <- strat_df |>
  filter(res_var < 0.5) |>
   pivot_longer(
      cols      = ends_with(" gene"),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    select(gene, gmean, gene_set, res_var, category) 
print(lvg_low, n = 300)



# 2. lvg bootstrapping

determine_hvg_numbers <- function(df, data_source){

  hvg_per_cluster <- df |>
    filter(gene_variability == "HVG") |>              
    group_by(cluster_id) |>
    summarise(nr_hvgs = n(), .groups = "drop")
  hvg_per_cluster

  return(hvg_per_cluster)
}


sample_once <- function(df, data_source, seed) {

  hvg_per_cluster_df <- determine_hvg_numbers(df, data_source)

  set.seed(seed)

  lvgs_sampled <- df |>
    filter(gene_variability == "LVG") |>
    inner_join(hvg_per_cluster_df, by = "cluster_id") |>
    group_by(cluster_id) |>
    slice(sample.int(n(), size = pmin(first(nr_hvgs), n()))) |>
    ungroup()
  #lvgs_sampled[1:6, 20:23]

  intermediate_sampled <- df |>
    filter(gene_variability == "Intermediate") |>
    inner_join(hvg_per_cluster_df, by = "cluster_id") |>
    group_by(cluster_id) |>
    slice(sample.int(n(), size = pmin(first(nr_hvgs), n()))) |>
    ungroup()

  lvg_and_inter_df <- bind_rows(lvgs_sampled, intermediate_sampled)

  return(lvg_and_inter_df)
}



gene_variability_group_bootstrapping  <- function(df, data_source, R, seed){

  set.seed(seed) 
  seeds <- sample.int(.Machine$integer.max, R)
  all_runs <- data.frame()
      
  for (run in seq_len(R)) {

    # subsample 
    one_run <- sample_once(df, data_source, seeds[run]) |>
      mutate(run_id = run)

    subsampled_df <- df |>
      filter(!gene_variability %in% c("LVG", "Intermediate")) |>
      bind_rows( one_run |> select(names(df)))


    # do stratification
    sampled_strat_df <- stratify_df(df = subsampled_df, data_source = "droplet", cell_type_selection = "Macrophage") |>
                        mutate(run_id = run)
    
    all_runs <- bind_rows(all_runs, sampled_strat_df)
  }

return(all_runs)
}

configure_bootstrap <- function(){
 # for each cluster run the bootstrapping 

  df <- readRDS("droplet/data/gene_set_df.rds")
  data_source <- "droplet"

  R <- 100         
  seed <- 45

  bootstrapped_df <- df   |>
    group_split(cluster_id, .keep = TRUE) |>
    map_dfr(~ gene_variability_group_bootstrapping(.x, data_source, R, seed))
  
  check <- bootstrapped_df |>
    group_by(gene_variability, cluster_id) |>
    count() |>
    arrange(cluster_id)
  print(check)

  saveRDS(bootstrapped_df, file.path(data_source, "data", "strat_subsampled_df.rds"))

  return(bootstrapped_df)
}



check_bootstrapping <- function(){

  per_run <- bootstrapped_df |>
    pivot_longer(
      cols      = ends_with(" gene"),
      names_to  = "gene_set",
      values_to = "in_set") |>
    filter(in_set == 1) |>
    count(run_id, cluster_id, category, name = "n_genes") |>
    group_by(cluster_id) |>
    complete(run_id, category, fill = list(n_genes = 0L)) |>  # treat absent as 0
    ungroup()


  stats <- per_run |>
    group_by(cluster_id, category) |>
    summarise(
      mean_n = mean(n_genes),
      sd_n   = sd(n_genes),
      .groups = "drop")

    all_cats <- c(
    "LVG 1", "LVG 2", "LVG 3", "LVG 4", "LVG 5", "LVG 6",
    "HVG 1", "HVG 2", "HVG 3", "HVG 4", "HVG 5", "HVG 6", 
    "Intermediate 1",  "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5", "Intermediate 6", 
    "Not expressed 0")

  # color scheme
  col_key <- read_excel("color_scheme_categories.xlsx") |>
    mutate(col_category = factor(col_category, levels = all_cats))
  col_vec <- setNames(col_key$hex_code, col_key$col_category)
  stats$category <- factor(stats$category, levels = all_cats)

  plot <- ggplot(stats, aes(category, mean_n, fill = category)) +
    geom_col() +
    geom_errorbar(aes(ymin = pmax(mean_n - sd_n, 0), ymax = mean_n + sd_n),
                width = 0.2, linewidth = 0.5) +
    facet_wrap(~ cluster_id, scales = "free_y") +
    labs(x = "Category", y = "Number of genes") +
     scale_fill_manual(values = col_vec, limits = all_cats, drop = FALSE,
                    guide = guide_legend(nrow = 6)) +
    theme_classic() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 20),
      axis.text.y        = element_text(size = 20),
      axis.title.x        = element_text(size = 20),
      axis.title.y        = element_text(size = 20),
      strip.text.y.left  = element_blank(),
      strip.text.y.right = element_text(angle = 0, size = 20),
      strip.text = element_text(size = 20),
      strip.background.y = element_blank(),
      legend.position    = "top", size = 20,
      legend.text =     element_text(size = 20))
      

  if (data_source == "facs") {
    print("Saved to: FACS.")
    ggsave("facs/plots/3_m/Test_16/lvg_strap.png", plot, width = 16, height = 17)
  } else if (data_source == "droplet") {
    print("Saved to: droplet")
    ggsave("droplet/plots/3_m/Test_16/lvg_strap.png", plot, width = 16, height = 17)
  } else {
    stop("Issue when saving.")
  }    
}






# 3. distribution of resvar and fraction of gene sets to decide what cutoff to use

# if more enriched in HK the further down then we should not cut it off, and same for the immune response genes

data_source <- "droplet"

# bootstrapped each group: then calc the fraction and enrichment for gene sets

# fraction
plot_fraction_gene_set <- function() {
    
  df <- bootstrapped_df

  bins <- seq(0, 10, by = 0.05) 

  per_run_hist <- df %>%
    group_by(run_id) %>%
    summarise(res_var = list(res_var[is.finite(res_var)]), .groups = "drop") %>%
    unnest(res_var) %>%
    mutate(bin = cut(res_var, breaks = bins, include.lowest = TRUE)) %>%
    count(run_id, bin, name = "n_bin") %>%
    group_by(run_id) %>%
    complete(bin, fill = list(n_bin = 0L)) %>% 
    ungroup()


  hist_stats <- per_run_hist %>%
    group_by(bin) %>%
    summarise(
      mean_n = mean(n_bin),
      sd_n   = sd(n_bin),
      .groups = "drop") %>%
    # midpoint for plotting on a numeric x-axis (optional but nicer)
    mutate(bin_mid = (as.numeric(sub("\\((.+),.*", "\\1", bin)) +
                      as.numeric(sub(".*,(.+)]", "\\1", bin))) / 2)


  plot <- ggplot(hist_stats, aes(x = bin_mid, y = mean_n )) +
      geom_col(width = diff(bins)[1]) +
      geom_errorbar(aes(ymin = pmax(mean_n - sd_n, 0), ymax = mean_n + sd_n), width = 0) +
      #scale_fill_manual(values = c("TRUE" = "#5195D3", "FALSE" = "lightgrey")) +
      labs(x = "Residual variance",
          y = "Count of genes") +
      coord_cartesian(xlim = c(0, 10)) + 
      geom_vline(aes(xintercept = 1), 
              colour = "black", linetype = "dashed", size = 1) +
    theme_classic() +
    theme(
        strip.background      = element_blank(),
        axis.text.x           = element_text(size = 20),
        axis.text.y           = element_text(size = 20),
        axis.title.x      = element_text(size = 20),
        axis.title.y      = element_text(size = 20))


    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_16/res_var_HK.png", plot, width = 12, height = 5)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_16/res_var_HK.png", plot, width = 12, height = 5)
    } else {
      stop("Issue when saving.")
    }   

}
colnames(bootstrapped_df)

df <- bootstrapped_df

bins <- seq(0, 10, by = 0.05)
bin_width <- diff(bins)[1]
bin_centers <- bins[-length(bins)] + bin_width/2

per_run_hist <- df %>%
  filter(is.finite(res_var)) %>%
  mutate(Housekeeping = gene %in% `Housekeeping gene`,
         bin = cut(res_var, breaks = bins, include.lowest = TRUE)) %>%
  count(run_id, Housekeeping, bin, name = "n_bin") %>%
  group_by(run_id, Housekeeping) %>%
  complete(bin, levels = levels(bin), fill = list(n_bin = 0L)) %>%
  ungroup()

hist_stats <- per_run_hist %>%
  group_by(Housekeeping, bin) %>%
  summarise(mean_n = mean(n_bin), sd_n = sd(n_bin), .groups = "drop") %>%
  mutate(bin_mid = bin_centers[as.integer(bin)])

plot <- ggplot(hist_stats, aes(x = bin_mid, y = mean_n, fill = hk)) +
  geom_col(position = position_dodge(width = bin_width), width = bin_width) +
  geom_errorbar(aes(ymin = pmax(mean_n - sd_n, 0), ymax = mean_n + sd_n),
                position = position_dodge(width = bin_width), width = 0) +
  scale_fill_manual(values = c(`TRUE` = "#1f78b4", `FALSE` = "grey70"),
                    labels = c(`TRUE` = "Housekeeping", `FALSE` = "Non-HK")) +
  labs(x = "Residual variance", y = "Genes per bin (mean ± SD across runs)", fill = "") +
  coord_cartesian(xlim = c(0, 10)) +
  theme_classic()






hist_df <- lvg_test_df |>
    mutate(Housekeeping = gene_set == "Housekeeping gene") |>
    mutate(`Immune response` = gene_set == "Immune response gene") |>
    mutate(bin = cut(res_var, breaks = seq(0, 10, by = 0.05)))

# enrichment per bin

p_global <- mean(hist_df$Housekeeping)
enrichment_df <- hist_df %>%
  group_by(bin) %>%
  summarise(
    n_bin = n(),
    n_hk = sum(Housekeeping),
    frac_hk = n_hk / n_bin,
    enrichment = frac_hk / p_global
  )

unique(enrichment_df$bin)

plot <- ggplot(enrichment_df, aes(x = bin, y = enrichment)) +
    geom_col(fill = "#5195D3") +
    labs(x = "Residual variance bins",
         y = "Enrichment") +
  theme_classic() +
  theme(
      strip.background      = element_blank(),
      axis.text.x           = element_text(angle = 45, hjust = 1, size = 15),
      axis.text.y           = element_text(size = 20),
      axis.title.x      = element_text(size = 20),
      axis.title.y      = element_text(size = 20)) +
  scale_x_discrete(breaks = levels(enrichment_df$bin)[seq(1, length(levels(enrichment_df$bin)), by = 10)])

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_16/res_var_HK_enrich.png", plot, width = 12, height = 5)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_16/res_var_HK_enrich.png", plot, width = 12, height = 5)
    } else {
      stop("Issue when saving.")
    }   






# enrichment per bin
p_global <- mean(hist_df$`Immune response`)
enrichment_df <- hist_df %>%
  group_by(bin) %>%
  summarise(
    n_imm = n(),
    n_bin = sum(`Immune response`),
    frac_imm = n_imm / n_bin,
    enrichment = frac_imm / p_global
  )

unique(enrichment_df$bin)

plot <- ggplot(enrichment_df, aes(x = bin, y = enrichment)) +
    geom_col(fill = "#EAB448") +
    labs(x = "Residual variance bins",
         y = "Enrichment") +
  theme_classic() +
  theme(
      strip.background      = element_blank(),
      axis.text.x           = element_text(angle = 45, hjust = 1, size = 15),
      axis.text.y           = element_text(size = 20),
      axis.title.x      = element_text(size = 20),
      axis.title.y      = element_text(size = 20)) +
  scale_x_discrete(breaks = levels(enrichment_df$bin)[seq(1, length(levels(enrichment_df$bin)), by = 10)])

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_16/res_var_imm_enrich.png", plot, width = 12, height = 5)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_16/res_var_imm_enrich.png", plot, width = 12, height = 5)
    } else {
      stop("Issue when saving.")
    }   


# fraciton immune
plot <- ggplot(hist_df, aes(x = res_var, fill = `Immune response` )) +
    geom_histogram(binwidth = 0.05) +
    scale_fill_manual(values = c("TRUE" = "#EAB448", "FALSE" = "lightgrey")) +
    labs(x = "Residual variance",
         y = "Count of genes") +
    coord_cartesian(xlim = c(0, 10)) + 
    geom_vline(aes(xintercept = 1), 
             colour = "black", linetype = "dashed", size = 1) +
  theme_classic() +
  theme(
      strip.background      = element_blank(),
      axis.text.x           = element_text(size = 20),
      axis.text.y           = element_text(size = 20),
      axis.title.x      = element_text(size = 20),
      axis.title.y      = element_text(size = 20))

    if (data_source == "facs") {
      print("Saved to: FACS.")
      ggsave("facs/plots/3_m/Test_16/res_var_immune.png", plot, width = 12, height = 5)
    } else if (data_source == "droplet") {
      print("Saved to: droplet")
      ggsave("droplet/plots/3_m/Test_16/res_var_immune.png", plot, width = 12, height = 5)
    } else {
      stop("Issue when saving.")
    }   


