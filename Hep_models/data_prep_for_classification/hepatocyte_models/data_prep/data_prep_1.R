library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(readr)
library(tidyr)

# this is the script for predictor prep from the original script (all ages together), last edited 17th june 2025

##PART I
## took og code and put it into functions, but the main approach stays the same

promoter_mofiv_prep <- function(){
    # if promoter has multiple promoter motives: calc mean of motiv 
    # if gene has mulitple promoters: calc number of  promoters 

    pred<-read_tsv(file = 'data_prep_for_classification/general_models/predictors/promoter_elements.tsv')
    # H: hanged mean to sum of promoters
    pred_master <- pred %>% 
        group_by(gene) %>% 
        summarize(TATA_box = mean(TATA_box), 
                  Inr = mean(Inr), 
                  CCAAT_box = mean(CCAAT_box), 
                  GC_box = mean(GC_box),
                  count = n()) %>% 
        drop_na() 

    write_tsv(pred_master, file =  'data_prep_for_classification/hepatocyte_models/results/predictors_non_amb.tsv')
}
promoter_mofiv_prep()




assign_variability <- function(){
  # in different clusters the gene can be differently variable: find uniform way
  hep <-  read_delim("data_prep_for_classification/hepatocyte_models/data/hepatocytes_master.txt")
  
  hep_master_df <- hep %>%                            
      mutate(age_sex = age,
              age = case_when(
              str_detect(age_sex, "18m") ~ "18m",
              str_detect(age_sex, "24m") ~ "24m",
              str_detect(age_sex, "3m")  ~ "3m"),
          cluster_id = paste0(tissue,"_" ,age_sex, "_", cluster)) %>%
      select(gene, cluster_id, age, res_var, gmean, variability) %>%
      group_by(gene) %>%                     
      mutate( # in og she had summarize here
        var_non_amb = case_when(
      any(variability == "HVG") ~ "HVG",
      mean(variability == "LVG") >= 0.75 ~ "LVG",
      TRUE ~ "Inconclusive")) %>%
      ungroup() %>%
      group_by(gene, var_non_amb) %>%
      summarise(median_gmean = median(gmean), .groups = "drop") 

  return(hep_master_df)
}
hep_non_ambiguous_df <- assign_variability()



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
    filter(var_non_amb %in% c("HVG", "LVG"))

    
  hvg <- df %>% filter(var_non_amb == "HVG")
  lvg <- df %>% filter(var_non_amb == "LVG")

  n_matches <- min(nrow(hvg), nrow(lvg))

  hvg_pool <- hvg %>% slice_sample(n = n_matches) # random sampling without replacement according to setseed
  lvg_sel <- lvg   
  lvg_matches <- vector("list", n_matches)


    for (row in 1:nrow(hvg_pool)) {
      set.seed(1) # always selects the same row if there are ties, gets the first matching gene
      value <- hvg_pool[row,]$median_gmean # set gmean value of hvg
      
      #identify best match in lvgs
      lvg_match <- lvg_sel %>% 
        filter(abs(median_gmean - value) == min(abs(median_gmean - value))) %>%
        sample_n(1)
      # remove match from lvgs
      lvg_sel <- lvg_sel %>% filter(!gene %in% lvg_match$gene)
      
      if(row == 1) {
        lvg_matches <- lvg_match
      } else {
        lvg_matches <- rbind(lvg_matches, lvg_match)
      }
    }
    lvg_matched <- bind_rows(lvg_matches)

    # test bias to certian positions in lvg pool
    match_list <- match(lvg_matched$gene, lvg$gene)
    match_clean <- match_list[!is.na(match_list)]
    match_df <- data.frame(position = match_clean)
    print(ggplot(match_df, aes(x = position)) + geom_histogram(binwidth = 50, fill = "steelblue", colour = "white") + labs(x = "Position in pool of possible LVG matches", y = "Count") + theme_classic())
    ggtitle(paste("LVG Match Position Distribution"))
    ggsave('data_prep_for_classification/hepatocyte_models/results/sampling_distribution_H.png')

    final_matched <- bind_rows(hvg[1:nrow(lvg_matched), ], lvg_matched) %>%
      mutate(var_non_amb = factor(var_non_amb, levels = c("HVG", "LVG")),) 

  return(final_matched)
}
matched_df  <- gmean_matching_and_subsample(hep_and_prom_df)


gmean_matching_and_subsample_improved <- function(hep_and_promoters_df){

  df <- hep_and_promoters_df %>%
    filter(var_non_amb %in% c("HVG", "LVG"))

  hvg <- df %>% filter(var_non_amb == "HVG")
  lvg <- df %>% filter(var_non_amb == "LVG")

  set.seed(1)
  n_matches <- min(nrow(hvg), nrow(lvg))
  hvg_pool <- hvg %>% slice_sample(n = n_matches) # random sampling without replacement according to setseed
  lvg_pool <- lvg   

  lvg_matches <- vector("list", n_matches)

  for (i in seq_len(n_matches)) {
    val <- hvg_pool$median_gmean[i]
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
  print(ggplot(match_df, aes(x = position)) + geom_histogram(binwidth = 50, fill = "steelblue", colour = "white") + labs(x = "Position in pool of possible LVG matches", y = "Count") + theme_classic())
  ggtitle(paste("LVG Match Position Distribution"))
  ggsave('data_prep_for_classification/hepatocyte_models/results/sampling_distribution_H.png')

  final_matched <- bind_rows(hvg[1:nrow(lvg_matched), ], lvg_matched) %>%
    mutate(var_non_amb = factor(var_non_amb, levels = c("HVG", "LVG")),) 

  return(final_matched)
}
      
matched_df  <- gmean_matching_and_subsample_improved(hep_and_prom_df)
table(matched_df$age, matched_df$var_non_amb)



save_matched_df <- function(matched_df){
    write_tsv(matched_df, file = 'data_prep_for_classification/hepatocyte_models/results/matched_df.tsv')
}   
save_matched_df(matched_df)



check_match_pairs <- function(matched_df){
  # see how well the median gmeans of HVG and LVG match 
  hvg_df <- matched_df %>% filter(var_non_amb == "HVG") 
  lvg_df <- matched_df %>% filter(var_non_amb == "LVG") 

  hist <- ggplot(matched_df, aes(x = median_gmean, fill = var_non_amb)) +
    geom_histogram(position = "identity", alpha = 0.6, binwidth = 0.2, colour = "white") +
    scale_fill_manual(values = c("HVG" = "darkred", "LVG" = "orange")) +
    labs(
      title = "Median gmeans of matched HVG/LVG pairs",
      x = "Median gmean",
      y = "Number of genes",
      fill = "Variability class" ) +
    theme_classic()

  ggsave(("data_prep_for_classification/hepatocyte_models/results/hist_hvg_lvg.png"), plot = hist)
}
match <-  check_match_pairs(matched_df)


check_random_pairs <- function(hep_and_prom_df){
  # see how a random sampling makes the median gmeans of HVG and LVG match 
   
  hvg_df <- hep_and_prom_df %>% filter(var_non_amb == "HVG") 
  lvg_df <- hep_and_prom_df %>% filter(var_non_amb == "LVG") 

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

  ggsave(("data_prep_for_classification/hepatocyte_models/results/hist_random.png"), plot = hist)
}
random <- check_random_pairs(hep_and_prom_df)








##PART II
## extracted the relavant bits of the og code (non_amb_HVG_LVG_prep.R)

out_folder <- 'data_prep_for_classification/hepatocyte_models/results'
if (!file.exists(out_folder)) dir.create(out_folder)

pred_master<-read_tsv(file ='data_prep_for_classification/hepatocyte_models/results/predictors_non_amb.tsv')
hep<- read_delim('data_prep_for_classification/hepatocyte_models/data/hepatocytes_master.txt')


hep_master<- hep %>%
    dplyr::select(gene, cluster, tissue, age, tissue_cluster, res_var, gmean, variability)%>%
    mutate(age_sex=age, 
        age = case_when( 
        str_detect(age_sex, "18m") ~ "18m",
        str_detect(age_sex, "24m") ~ "24m",
        str_detect(age_sex, "3m") ~ "3m"))


hep_non_amb<- hep_master %>% 
   group_by(gene) %>%
   summarise(
   var = case_when(
     any(variability == "HVG") ~ "HVG",
     mean(variability == "LVG") >= 0.75 ~ "LVG",
      TRUE ~ "Inconclusive"))


hep_non_amb_pred <- hep_non_amb %>%
    left_join(pred_master, by = c('gene' = 'gene'))

hep_non_amb_pred_filt <- hep_non_amb_pred %>% filter(!is.na(TATA_box)) #H: ??
hep_non_amb_pred_filt %>% group_by(var) %>% summarise(count=n())
table(hep_non_amb_pred_filt$var)


hvg<- hep_non_amb_pred_filt %>% filter(var=='HVG') %>% pull(gene)
lvg<- hep_non_amb_pred_filt %>% filter(var=='LVG') %>% pull(gene)
incon<-hep_non_amb_pred_filt %>% filter(var=='Inconclusive') %>% pull(gene)

hep_master_non_amb <- hep_master %>% 
    mutate(var_non_amb = case_when(
        gene %in% hvg ~ 'HVG',
        gene %in% lvg ~'LVG',
        gene %in% incon ~ 'Intermediate')) %>% 
    filter(!is.na(var_non_amb))
table(hep_master_non_amb$var_non_amb)

#median gmean by var group
med_gmean_df<- hep_master_non_amb %>%
    group_by(gene, var_non_amb)%>%
    summarise(med_gmean=median(gmean))


#filter out intermediate and downsample based on median gmean
hep_non_amb_ds<- hep_non_amb_pred_filt %>%
    filter(var!='Inconclusive')%>%
    left_join(med_gmean_df, by = c('gene'='gene', 'var'='var_non_amb'))

hvg_sel <- hep_non_amb_ds %>% filter(var == "HVG")
lvg_sel <- hep_non_amb_ds %>% filter(var == "LVG")

for (row in 1:nrow(hvg_sel)) {
  set.seed(1) # always selects the same row if there are ties, gets the first matching gene
  value <- hvg_sel[row,]$med_gmean # set gmean value of hvg
  
  #identify best match in lvgs
  lvg_match <- lvg_sel %>% 
    filter(abs(med_gmean - value) == min(abs(med_gmean - value))) %>%
    sample_n(1)
  # remove match from lvgs
  lvg_sel <- lvg_sel %>% filter(!gene %in% lvg_match$gene)
  
  if(row == 1) {
    lvg_matches <- lvg_match
  } else {
    lvg_matches <- rbind(lvg_matches, lvg_match)
  }
}

prediction_data_matched <- rbind(hvg_sel, lvg_matches) %>%
  mutate(var = factor(var, levels = c("HVG", "LVG")))

head(prediction_data_matched)

#write_tsv(prediction_data_matched, file = 'data_prep_for_classification/hepatocyte_models/results/hvg_lvg_predictors_counts.tsv')

# H: checking if the set_seed() biases the selected lvgs in a certain position of the df
match_list <- match(lvg_matches$gene, lvg_sel$gene)
match_clean <- match_list[!is.na(match_list)]
match_df <- data.frame(position = match_clean)

p <- ggplot(match_df, aes(x = position)) +
  geom_histogram(binwidth = 50, fill = "steelblue", colour = "white") +
  labs(x = "Position in pool of possible LVG matches",
       y = "Count") +
  theme_classic()
ggsave('data_prep_for_classification/hepatocyte_models/results/sampling_distribution.png', p)



# comparison of lvg and hvgs selected
e <- read_tsv( file = 'data_prep_for_classification/hepatocyte_models/results/hvg_lvg_predictors_counts.tsv')
h <- read_tsv(file = 'data_prep_for_classification/hepatocyte_models/results/matched_df.tsv')

# same dimensions
# HVGs match
# LVGs dont (only 792): different sampling method
match_list <- match(e$gene[e$var == 'LVG'], h$gene[h$var_non_amb == "LVG"])
sum(!is.na(match_list))




