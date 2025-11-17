#preparation of predictors
# aggregate predictor information: if more than one promoter is annotated, create a masterannotation 
#(if any annotated promoter has a given element for a given gene, it's positive, if none has it, it's negative)
#there is 2 promoters with NA in genes so drop na
#plus also add a count of how many elements are assigned to that gene


library(tidyverse)
library(dplyr)
library(stringr)
library(tibble)
library(circlize)
library(ranger)

out_folder <- 'regulated_noise/06_24_facs_vs_droplet_analysis/promoter_prediction/promoter_data_explor'
if (!file.exists(out_folder)) dir.create(out_folder)

pred<-read_tsv(file = './regulated_noise/06_24_facs_vs_droplet_analysis/promoter_prediction/promoter_elements.tsv')

#preparation of predictors
pred_master <- pred %>% 
    group_by(gene) %>% 
    summarize(TATA_box = mean(TATA_box), 
                Inr = mean(Inr), 
                CCAAT_box = mean(CCAAT_box), 
                GC_box = mean(GC_box),
                count = n()) %>% 
    drop_na()

write_tsv(pred_master, file = paste0(out_folder, '/predictors_non_amb.tsv'))

hep<- read_delim('./regulated_noise/06_24_facs_vs_droplet_analysis/gsea_all_ages/overlap_gsea/hepatocytes_master.txt')

hep_master<- hep %>%
    select(gene, cluster, tissue, age, tissue_cluster, res_var, gmean, variability)%>%
    filter(variability!='intermediate')

#ambigous genes (genes sometimes labelled LVG sometimes HVG)

hep_amb <- hep_master %>%
  group_by(gene) %>%
  summarise(unique_classes = n_distinct(variability),
  classes = paste(sort(unique(variability)), collapse = ", ")) %>%
  mutate(status = case_when(
    unique_classes == 1 & classes == "HVG" ~ "Always HVG",
    unique_classes == 1 & classes == "LVG" ~ "Always LVG",
    unique_classes > 1 ~ "Mixed"
  ))

#all the hepatocytes
hep_counts<- hep_amb %>% 
    group_by(status)%>%
    summarise(count = n())

p<-ggplot(hep_counts, aes(x = status, y = count)) +
        geom_bar(stat = "identity", color = "black", fill = 'darkgrey') +
        labs(title = "Distribution of Gene Variability in all Hepatocytes",
            x = "",
            y = "Count of Genes") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")
ggsave(p, filename = paste0(out_folder, '/hep_LVG_HVG_overlaps.png'), height = 3, width = 5)

#separated to age groups
ages<-unique(hep_master$age)
for (age_group in ages){
    hep_master_temp<-hep_master %>%
        filter(age == age_group)

    hep_amb <- hep_master_temp %>%
        group_by(gene) %>%
        summarise(unique_classes = n_distinct(variability),
        classes = paste(sort(unique(variability)), collapse = ", ")) %>%
        mutate(status = case_when(
            unique_classes == 1 & classes == "HVG" ~ "Always HVG",
            unique_classes == 1 & classes == "LVG" ~ "Always LVG",
            unique_classes > 1 ~ "Mixed")) %>%
        group_by(status)%>%
        summarise(count = n())

    hep_amb<- hep_amb %>%
        filter(status!='Always LVG')

    p<-ggplot(hep_amb, aes(x = status, y = count, fill = status)) +
        geom_bar(stat = "identity", color = "black") +
        labs(title = paste0("Distribution of Gene Variability in ", age_group),
            x = "",
            y = "Count of Genes") +
        scale_fill_brewer(palette = "Set2") +  
        theme_classic() +
        theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")
    ggsave(p, filename = paste0(out_folder, '/LVG_HVG_overlaps_',age_group,'_filt.png'), height = 3, width = 5)
}

#further inv into the mixed groups
hep_mixed<-tibble()
ages<-unique(hep_master$age)
for (age_group in ages){
    hep_master_temp<-hep_master %>%
        filter(age == age_group)
    goi<- hep_master_temp %>% group_by(gene, variability) %>% summarise(vount=n()) %>% filter(variability == 'HVG') %>% pull(gene)
    hep_mixed_temp<-hep_master_temp %>% group_by(gene, variability) %>% summarise(count=n()) %>% filter(gene %in% goi) %>% mutate(age = age_group)
    hep_mixed<-rbind(hep_mixed, hep_mixed_temp)
}

ages<-unique(hep_master$tissue_cluster)
for (age_group in ages){
    hep_master_temp<-hep_master %>%
        filter(tissue_cluster == age_group)

    hep_amb <- hep_master_temp %>%
        group_by(gene) %>%
        summarise(unique_classes = n_distinct(variability),
        classes = paste(sort(unique(variability)), collapse = ", ")) %>%
        mutate(status = case_when(
            unique_classes == 1 & classes == "HVG" ~ "Always HVG",
            unique_classes == 1 & classes == "LVG" ~ "Always LVG",
            unique_classes > 1 ~ "Mixed")) %>%
        group_by(status)%>%
        summarise(count = n())

    p<-ggplot(hep_amb, aes(x = status, y = count, fill = status)) +
        geom_bar(stat = "identity", color = "black") +
        labs(title = paste0("Distribution of Gene Variability in ", age_group),
            x = "",
            y = "Count of Genes") +
        scale_fill_brewer(palette = "Set2") +  
        theme_classic() +
        theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")
    ggsave(p, filename = paste0(out_folder, '/LVG_HVG_overlaps_',age_group,'_filt.png'), height = 3, width = 5)
}

#splitting by age only and disregad of sex
hep_master<- hep_master %>% mutate(age_sex=age, 
    age = case_when(  # Extract the age portion
      str_detect(age_sex, "18m") ~ "18m",
      str_detect(age_sex, "24m") ~ "24m",
      str_detect(age_sex, "3m") ~ "3m"))

ages<-unique(hep_master$age)
for (age_group in ages){
    hep_master_temp<-hep_master %>%
        filter(age == age_group)

    hep_amb <- hep_master_temp %>%
        group_by(gene) %>%
        summarise(unique_classes = n_distinct(variability),
        classes = paste(sort(unique(variability)), collapse = ", ")) %>%
        mutate(status = case_when(
            unique_classes == 1 & classes == "HVG" ~ "Always HVG",
            unique_classes == 1 & classes == "LVG" ~ "Always LVG",
            unique_classes > 1 ~ "Mixed")) %>%
        group_by(status)%>%
        summarise(count = n())

    hep_amb <- hep_amb %>% filter(status!='Always LVG')

    p<-ggplot(hep_amb, aes(x = status, y = count, fill = status)) +
        geom_bar(stat = "identity", color = "black") +
        labs(title = paste0("Distribution of Gene Variability in ", age_group),
            x = "",
            y = "Count of Genes") +
        scale_fill_brewer(palette = "Set2") +  
        theme_classic() +
        theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")
    ggsave(p, filename = paste0(out_folder, '/LVG_HVG_overlaps_',age_group,'_filt.png'), height = 3, width = 5)
}

#number of cluster per age group
hep_cl_count<- hep_master %>% group_by(age) %>% summarise(n_cluster=n_distinct(cluster))
p<-ggplot(hep_cl_count, aes(x = age, y = n_cluster, fill = age)) +
        geom_bar(stat = "identity", color = "black") +
        labs(title = 'Number of clusters per Age',
            x = "",
            y = "Number of clusters") +
        scale_fill_brewer(palette = "Set2") +  
        theme_classic() +
        theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")
ggsave(p, filename = paste0(out_folder, '/cluster_counts_per_age.png'), height = 2, width = 3)
