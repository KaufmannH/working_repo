# 1. preparation of predictors
# 2. exploring 
# aggregate predictor information: if more than one promoter is annotated, create a masterannotation 
#(if any annotated promoter has a given element for a given gene, it's positive, if none has it, it's negative)
#there is 2 promoters with NA in genes so drop na
#plus also add a count of how many elements are assigned to that gene



library(dplyr)
library(stringr)
library(tibble)
library(circlize)
library(readr)
library(tidyr)
library(ggplot2)


out_folder <- 'data_prep_for_classification/hepatocyte_models/results'
if (!file.exists(out_folder)) dir.create(out_folder)

pred<-read_tsv(file = 'data_prep_for_classification/general_models/predictors/promoter_elements.tsv')
pred

#preparation of predictors
# proportion of promoter motifs and number of motivs per gene
pred_master <- pred %>% 
    group_by(gene) %>% 
    summarize(TATA_box = mean(TATA_box), 
                Inr = mean(Inr), 
                CCAAT_box = mean(CCAAT_box), 
                GC_box = mean(GC_box),
                count = n()) %>% 
    drop_na() 

pred_master
write_tsv(pred_master, file = paste0(out_folder, '/predictors_non_amb.tsv'))


# hepatocytes
hep<- read_delim('data_prep_for_classification/hepatocyte_models/data/hepatocytes_master.txt')
hep
hep_master<- hep %>%
    select(gene, cluster, tissue, age, tissue_cluster, res_var, gmean, variability)%>%
    filter(variability!='intermediate')
hep_master
#ambigous genes (genes sometimes labelled LVG sometimes HVG)

# tag genes with variable variability classes
hep_amb <- hep_master %>%
  group_by(gene) %>%
  summarise(unique_classes = n_distinct(variability),
  classes = paste(sort(unique(variability)), collapse = ", ")) %>%
  mutate(status = case_when(
    unique_classes == 1 & classes == "HVG" ~ "Always HVG",
    unique_classes == 1 & classes == "LVG" ~ "Always LVG",
    unique_classes > 1 ~ "Mixed"
  ))

# stats of ambiguous variability classes for all heps
hep_counts<- hep_amb %>% 
    group_by(status)%>%
    summarise(count = n())
hep_counts
p<-ggplot(hep_counts, aes(x = status, y = count)) +
        geom_bar(stat = "identity", color = "black", fill = 'darkgrey') +
        labs(title = "Distribution of Gene Variability in all Hepatocytes",
            x = "",
            y = "Count of Genes") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")
ggsave(p, filename = paste0(out_folder, '/var_switching/hep_LVG_HVG_overlaps.png'), height = 3, width = 5)



# stats of ambiguous variability classes for different age+sex
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
    ggsave(p, filename = paste0(out_folder, '/var_switching/LVG_HVG_overlaps_',age_group,'_filt.png'), height = 3, width = 5)
}
hep_master


#further inv into the mixed groups (H: ?)
# get genes that are HVGs but can be lvgs in the same mouse in a different cluster
hep_mixed<-tibble()
ages<-unique(hep_master$age)
for (age_group in ages){

    hep_master_temp<-hep_master %>%
        filter(age == age_group)
    # get list of all hvgs
    goi<- hep_master_temp %>% 
        group_by(gene, variability) %>%
        summarise(vount=n()) %>% 
        filter(variability == 'HVG') %>% 
        pull(gene)

    hep_mixed_temp<-hep_master_temp %>% 
        group_by(gene, variability) %>% 
        summarise(count=n()) %>% 
        filter(gene %in% goi) %>% 
        mutate(age = age_group)

    hep_mixed<-rbind(hep_mixed, hep_mixed_temp)
}
hep_mixed

# stats of ambiguous variability classes for different clusters
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
    ggsave(p, filename = paste0(out_folder, '/var_switching/LVG_HVG_overlaps_',age_group,'_filt.png'), height = 3, width = 5)
}


# stats of ambiguous variability classes for age only (not sex)
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
    ggsave(p, filename = paste0(out_folder, '/var_switching/LVG_HVG_overlaps_',age_group,'_filt.png'), height = 3, width = 5)
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
