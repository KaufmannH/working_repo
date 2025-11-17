#define diff cutoffs for LVG and HVG
#attaches the predictors from EPD
#exploratory plots, to see if maybe the promoter elements are visibly different in the diff variability classes
#then aggregates these promoter elements by mean, if multiple are assigned to one genes, 
#also adds the count of how many elements were annotated to that gene

library(tidyverse)
library(tidymodels)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(tibble)
library(circlize)
library(ranger)

out_folder <- 'regulated_noise/06_24_facs_vs_droplet_analysis/promoter_prediction/iteration_plots'
if (!file.exists(out_folder)) dir.create(out_folder)

pred_master<-read_tsv(file ='./regulated_noise/06_24_facs_vs_droplet_analysis/promoter_prediction/promoter_data_explor/predictors_non_amb.tsv')

hep<- read_delim('./regulated_noise/06_24_facs_vs_droplet_analysis/gsea_all_ages/overlap_gsea/hepatocytes_master.txt')

#label a gene as HVG if it's HVG at least once and LVG which are always labelled LVG
#also join with pred_master
#so intermediate is not filtered

hep_master<- hep %>%
    dplyr::select(gene, cluster, tissue, age, tissue_cluster, res_var, gmean, variability)%>%
    mutate(age_sex=age, 
        age = case_when( 
        str_detect(age_sex, "18m") ~ "18m",
        str_detect(age_sex, "24m") ~ "24m",
        str_detect(age_sex, "3m") ~ "3m"))

#in the 1st case
# Assign "LVG" if at least 75% of rows are LVG and no row is HVG
#assign HVG if the gene is ever labelled as hvg
# Otherwise, assign "Inconclusive"

hep_non_amb<- hep_master %>% 
   group_by(gene)%>%
   summarise(
   var = case_when(
     any(variability == "HVG") ~ "HVG",
     mean(variability == "LVG") >= 0.75 ~ "LVG",
      TRUE ~ "Inconclusive"))

#in the 2nd case, LVG is same, HVG is labelled if HVG in at least 8% of clusters. -- which means at least 2 clusters
#566 HVGs in this case, 525 have predictors for them :))

# hep_non_amb<- hep_master %>% 
#     group_by(gene)%>%
#     summarise(
#     var = case_when(
#       any(variability == "HVG") ~ "HVG",
#       mean(variability == "LVG") >= 0.75 ~ "LVG",
#       TRUE ~ "Inconclusive"))

# hep_non_amb %>% group_by(var) %>% summarise(count=n())

# #in the 3rd case, every 3. cluster for the hvgs, rest stays same
# #306 hvgs, 290 after the ones without promoter info
# hep_non_amb<- hep_master %>% 
#     group_by(gene)%>%
#     summarise(
#     var = case_when(
#       mean(variability == "HVG") >=0.1 ~ "HVG",
#       mean(variability == "LVG") >= 0.75 ~ "LVG",
#       TRUE ~ "Inconclusive"))
# hep_non_amb %>% group_by(var) %>% summarise(count=n())

#attach the predictors 

#it also has the number of predictors assigned to that specific gene
hep_non_amb_pred<- hep_non_amb %>%
    left_join(pred_master, by = c('gene' = 'gene'))

hep_non_amb_pred_filt<- hep_non_amb_pred %>% filter(!is.na(TATA_box))

hep_non_amb_pred_filt %>% group_by(var) %>% summarise(count=n())

#exploratory plots
hep_hm <- hep_non_amb_pred_filt %>%
  filter(count<4)%>%
  mutate(descr = paste(gene, var, sep = ".")) %>%  #combine gene and variability class
  select(-gene, -var) %>%  
  column_to_rownames("descr") %>% 
  as.matrix() %>%
  t()  #transpose for heatmap input

#scale the count column so can be shown on the heatmap
hep_hm["count", ] <- (hep_hm["count", ] - min(hep_hm["count", ])) / 
                     (max(hep_hm["count", ]) - min(hep_hm["count", ]))

#prepare column annotations based on variability
col_annotation <- word(colnames(hep_hm), sep = "\\.", 2)  #extract variability
col_annotation[col_annotation == "HVG"] <- "HVG"
col_annotation[col_annotation == "LVG"] <- "LVG"
col_annotation[col_annotation == "Inconclusive"] <- "Inconc."

custom_colors <- colorRamp2(c(0, 0.4, 0.8), c("blue", "white", "red"))

#heatmap
pdf(file = paste0(out_folder,"/counts/Heatmap_predictors_clusters_3.pdf"), width = 5, height = 4)

Heatmap(hep_hm,
        name = " ",
        column_order = colnames(hep_hm),
        col = custom_colors,
        column_split = col_annotation, 
        show_column_names = FALSE,
        row_order = c("TATA_box", "Inr", "GC_box", "CCAAT_box", 'count'))

dev.off()

#boxplots
boxplot_list <- list()
for (pred in c("TATA_box", "Inr", "CCAAT_box", "GC_box", 'count')) {
  # Create a boxplot for the current predictor grouped by 'var'
  boxplot_list[[pred]] <- ggplot(hep_non_amb_pred_filt, 
                                 aes(x = var, y = .data[[pred]])) +
    geom_boxplot(outlier.size = 0.5, fill = 'lightgrey') +
    #geom_jitter(size=0.1, alpha=0.4)+
    theme_classic() +
    labs(x = "", y = '', title = paste(pred, "by variability class")) +
    ylim(1,4)+
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 14))
  
  # Save the plot to the output folder
  ggsave(boxplot_list[[pred]], filename = file.path(out_folder, paste0('counts/',pred, "_boxplot_by_var.png")), 
    width = 4, height = 3)
}

#number of HVG and LVG in hepatocytes
#then balance those numbers by subsampling the LVGs down to the HVG size -- and matching the average expression
#label the original hep_master genes in the clusters with the hvg, lvg labelling

hvg<- hep_non_amb_pred_filt %>% filter(var=='HVG') %>% pull(gene)
lvg<- hep_non_amb_pred_filt %>% filter(var=='LVG') %>% pull(gene)
incon<-hep_non_amb_pred_filt %>% filter(var=='Inconclusive') %>% pull(gene)

hep_master_non_amb<- hep_master %>% 
    mutate(var_non_amb = case_when(
        gene %in% hvg ~ 'HVG',
        gene %in% lvg ~'LVG',
        gene %in% incon ~ 'Intermediate')) %>% 
    filter(!is.na(var_non_amb))

#median gmean by var group
med_gmean_df<- hep_master_non_amb %>%
    group_by(gene, var_non_amb)%>%
    summarise(med_gmean=median(gmean))

p<-ggplot(med_gmean_df, aes(x = var_non_amb, y = med_gmean, fill = var_non_amb)) +
  geom_violin(trim = TRUE, alpha = 0.7) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8)+
  theme_classic() +
  labs(title = "Median gmean by Variability", x = "", y = "med. gmean (log10)") +
  scale_fill_brewer(palette = "Set2")+
  scale_y_log10()+
  theme(legend.position = 'none')

p

ggsave(p, filename = file.path(out_folder, "counts/violin_var_med_gmean.png"), 
    width = 3, height = 3)


#count per var_group barplot
count_df<- hep_master_non_amb %>%
    group_by(var_non_amb)%>%
    summarise(count = n_distinct(gene))

p<-ggplot(count_df, aes(x = var_non_amb, y = count, fill = var_non_amb)) +
        geom_bar(stat = "identity", color = "black") +
        labs(title = 'Gene Counts per Variability',
            x = "",
            y = "Counts") +
        scale_fill_brewer(palette = "Set2") +  
        theme_classic() +
        theme(axis.text.x = element_text(angle = 10, hjust = 1), legend.position = "none")
ggsave(p, filename = paste0(out_folder, '/counts/gene_count_per_var.png'), height = 3, width = 3)

#filter out intermediate and downsample based on median gmean
hep_non_amb_ds<- hep_non_amb_pred_filt %>%
    filter(var!='Inconclusive')%>%
    left_join(med_gmean_df, by = c('gene'='gene', 'var'='var_non_amb'))

hvg_sel <- hep_non_amb_ds %>% filter(var == "HVG")
lvg_sel <- hep_non_amb_ds %>% filter(var == "LVG")

for (row in 1:nrow(hvg_sel)) {
  set.seed(1)
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

# Check results
p<-ggplot(prediction_data_matched, aes(x = var, y = med_gmean, fill = var)) +
  geom_violin(trim = TRUE, alpha = 0.7) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8)+
  theme_classic() +
  labs(title = "Median gmean by Variability", x = "", y = "med. gmean (log10)") +
  scale_y_log10()+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.position = 'none')

ggsave(p, filename = file.path(out_folder, "/counts/violin_var_downsample.png"), 
    width = 3, height = 3)

count_df<- prediction_data_matched %>%
    group_by(var)%>%
    summarise(count = n_distinct(gene))

p<-ggplot(count_df, aes(x = var, y = count, fill = var)) +
        geom_bar(stat = "identity", color = "black") +
        labs(title = 'Gene Counts per Variability',
            x = "",
            y = "Counts") +
        scale_fill_brewer(palette = "Set3") +  
        theme_classic() +
        theme(axis.text.x = element_text(angle = 10, hjust = 1), legend.position = "none")
ggsave(p, filename = paste0(out_folder, '/counts/gene_count_per_var_downsample.png'), height = 3, width = 3)

write_tsv(prediction_data_matched, file = './regulated_noise/06_24_facs_vs_droplet_analysis/promoter_prediction/hvg_lvg_predictors_counts.tsv')

###
#sep for 3m and 18m 
###

hep_3m<- hep_master %>% filter(age == '3m')
hep_18m<-hep_master %>% filter(age == '18m')

hep_3m_non_amb<- hep_3m %>% 
    group_by(gene)%>%
    summarise(
    var = case_when(
      any(variability == "HVG") ~ "HVG",
      mean(variability == "LVG") >= 0.75 ~ "LVG",
      TRUE ~ "Inconclusive"))

hep_3m_non_amb %>% group_by(var) %>% summarise(count=n())

hep_18m_non_amb<- hep_18m %>% 
    group_by(gene)%>%
    summarise(
    var = case_when(
      any(variability == "HVG") ~ "HVG",
      mean(variability == "LVG") >= 0.75 ~ "LVG",
      TRUE ~ "Inconclusive"))

hep_18m_non_amb %>% group_by(var) %>% summarise(count=n())

#attach the predictors 
#196 hvg in 3m
hep_3m_non_amb_pred<- hep_3m_non_amb %>%
    left_join(pred_master, by = c('gene' = 'gene'))

hep_3m_non_amb_pred_filt<- hep_3m_non_amb_pred %>% filter(!is.na(TATA_box))
hep_3m_non_amb_pred_filt %>% group_by(var) %>% summarise(count=n_distinct(gene))

#1053 hvg in 18m
hep_18m_non_amb<- hep_18m_non_amb %>%
    left_join(pred_master, by = c('gene' = 'gene'))

hep_18m_non_amb_pred_filt<- hep_18m_non_amb %>% filter(!is.na(TATA_box))
hep_18m_non_amb_pred_filt %>% group_by(var) %>% summarise(count=n_distinct(gene))

#exploratory plots
hep_hm_3 <- hep_3m_non_amb_pred_filt %>%
  mutate(descr = paste(gene, var, sep = ".")) %>%  #combine gene and variability class
  select(-gene, -var) %>%  
  column_to_rownames("descr") %>% 
  as.matrix() %>%
  t()  #transpose for heatmap input

hep_hm_18 <- hep_18m_non_amb_pred_filt %>%
  mutate(descr = paste(gene, var, sep = ".")) %>%  #combine gene and variability class
  select(-gene, -var) %>%  
  column_to_rownames("descr") %>% 
  as.matrix() %>%
  t()

#prepare column annotations based on variability -- change between 3m and 18m
col_annotation <- word(colnames(hep_hm_18), sep = "\\.", 2)  #extract variability
col_annotation[col_annotation == "HVG"] <- "HVG"
col_annotation[col_annotation == "LVG"] <- "LVG"
col_annotation[col_annotation == "Inconclusive"] <- "Inconc."

custom_colors <- colorRamp2(c(0, 0.4, 0.8), c("blue", "white", "red"))

#heatmap
pdf(file = paste0(out_folder,"/18m/Heatmap_predictors_clusters.pdf"), width = 7, height = 4)

Heatmap(hep_hm_18,
        name = " ",
        column_order = colnames(hep_hm_18),
        col = custom_colors,
        column_split = col_annotation, 
        show_column_names = FALSE,
        row_order = c("TATA_box", "Inr", "GC_box", "CCAAT_box"))

dev.off()

#boxplots -- change between 3m and 18m
boxplot_list <- list()
for (pred in c("TATA_box", "Inr", "CCAAT_box", "GC_box")) {
  # Create a boxplot for the current predictor grouped by 'var'
  boxplot_list[[pred]] <- ggplot(hep_3m_non_amb_pred_filt, 
                                 aes(x = var, y = .data[[pred]])) +
    geom_boxplot(outlier.size = 0.5, fill = 'lightgrey') +
    theme_classic() +
    labs(x = "", y = '', title = paste(pred, "by variability class")) +
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 14))
  
  # Save the plot to the output folder
  ggsave(boxplot_list[[pred]], filename = file.path(out_folder, paste0('3m/',pred, "_boxplot_by_var_2.png")), 
    width = 4, height = 3)
}

#number of HVG and LVG in hepatocytes
#then balance those numbers by subsampling the LVGs down to the HVG size -- and matching the average expression
#label the original hep_master genes in the clusters with the hvg, lvg labelling

hvg<- hep_18m_non_amb_pred_filt %>% filter(var=='HVG') %>% pull(gene)
lvg<- hep_18m_non_amb_pred_filt %>% filter(var=='LVG') %>% pull(gene)
incon<-hep_18m_non_amb_pred_filt %>% filter(var=='Inconclusive') %>% pull(gene)

hep_master_non_amb<- hep_18m %>% 
    mutate(var_non_amb = case_when(
        gene %in% hvg ~ 'HVG',
        gene %in% lvg ~'LVG',
        gene %in% incon ~ 'Intermediate')) %>% 
    filter(!is.na(var_non_amb))

#median gmean by var group
med_gmean_df<- hep_master_non_amb %>%
    group_by(gene, var_non_amb)%>%
    summarise(med_gmean=median(gmean))

p<-ggplot(med_gmean_df, aes(x = var_non_amb, y = med_gmean, fill = var_non_amb)) +
  geom_violin(trim = TRUE, alpha = 0.7) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8)+
  theme_classic() +
  labs(title = "Median gmean by Variability", x = "", y = "med. gmean (log10)") +
  scale_fill_brewer(palette = "Set2")+
  scale_y_log10()+
  theme(legend.position = 'none')

ggsave(p, filename = file.path(out_folder, "18m/violin_var_med_gmean.png"), 
    width = 3, height = 3)


#count per var_group barplot
count_df<- hep_master_non_amb %>%
    group_by(var_non_amb)%>%
    summarise(count = n_distinct(gene))

p<-ggplot(count_df, aes(x = var_non_amb, y = count, fill = var_non_amb)) +
        geom_bar(stat = "identity", color = "black") +
        labs(title = 'Gene Counts per Variability',
            x = "",
            y = "Counts") +
        scale_fill_brewer(palette = "Set2") +  
        theme_classic() +
        theme(axis.text.x = element_text(angle = 10, hjust = 1), legend.position = "none")
ggsave(p, filename = paste0(out_folder, '/18m/gene_count_per_var.png'), height = 3, width = 3)

#filter out intermediate and downsample based on median gmean
hep_non_amb_ds<- hep_18m_non_amb_pred_filt %>%
    filter(var!='Inconclusive')%>%
    left_join(med_gmean_df, by = c('gene'='gene', 'var'='var_non_amb'))

hvg_sel <- hep_non_amb_ds %>% filter(var == "HVG")
lvg_sel <- hep_non_amb_ds %>% filter(var == "LVG")

for (row in 1:nrow(hvg_sel)) {
  set.seed(1)
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

#check results
p<-ggplot(prediction_data_matched, aes(x = var, y = med_gmean, fill = var)) +
  geom_violin(trim = TRUE, alpha = 0.7) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8)+
  theme_classic() +
  labs(title = "Median gmean by Variability", x = "", y = "med. gmean (log10)") +
  scale_y_log10()+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.position = 'none')

ggsave(p, filename = file.path(out_folder, "/18m/violin_var_downsample.png"), 
    width = 3, height = 3)

count_df<- prediction_data_matched %>%
    group_by(var)%>%
    summarise(count = n_distinct(gene))

p<-ggplot(count_df, aes(x = var, y = count, fill = var)) +
        geom_bar(stat = "identity", color = "black") +
        labs(title = 'Gene Counts per Variability',
            x = "",
            y = "Counts") +
        scale_fill_brewer(palette = "Set3") +  
        theme_classic() +
        theme(axis.text.x = element_text(angle = 10, hjust = 1), legend.position = "none")
ggsave(p, filename = paste0(out_folder, '/18m/gene_count_per_var_downsample.png'), height = 3, width = 3)

write_tsv(prediction_data_matched, file = './regulated_noise/06_24_facs_vs_droplet_analysis/promoter_prediction/hvg_lvg_predictors18m.tsv')
