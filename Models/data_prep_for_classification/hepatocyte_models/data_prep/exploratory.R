
library(dplyr)
library(stringr)
library(tibble)
library(circlize)
library(readr)
library(tidyr)
library(ggplot2)




# load df
heatmap_motivs <- function(){
  df_matched <-  read_tsv('data_prep_for_classification/hepatocyte_models/results/matched_df.tsv')

  prom_df <- df_matched |>
    pivot_longer(cols =c('TATA_box', 'CCAAT_box', 'GC_box', 'Inr'), names_to = "prom_element", values_to = "value") 

  prom2 <- prom_df |>
    group_by(var_non_amb, prom_element, cluster_id) |>
    mutate(mean_promoter = mean(value)) |>
    arrange(gene) |>
    select(gene, cluster_id, var_non_amb, mean_promoter, age)

  p <- ggplot(prom2, aes(x = var_non_amb, y = mean_promoter, fill = prom_element)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~age) +
    labs(
      x = "Promoter element",
      y = "Mean count",
      colour = "Variability class" ) +
    theme_classic()

  ggsave('data_prep_for_classification/hepatocyte_models/results/predictors/prom.png', p)
}


motiv_expression <- function(){
  df_matched <-  read_tsv('data_prep_for_classification/hepatocyte_models/results/matched_df.tsv')
  df_matched <- write_tsv(prediction_data_matched, file = 'data_prep_for_classification/hepatocyte_models/results/hvg_lvg_predictors_counts.tsv')

  prom_df <- df_matched |>
    pivot_longer(cols =c('TATA_box', 'CCAAT_box', 'GC_box', 'Inr'), names_to = "prom_element", values_to = "value") 



  p <- ggplot(prom_df, aes(x = value , y = gmean, color = var_non_amb)) +
    geom_point(alpha = 0.7) +
    facet_wrap(~prom_element) +
    labs(
      x = "Promoter element",
      y = "Mean count",
      colour = "Variability class" ) +
    theme_classic()

  ggsave('data_prep_for_classification/hepatocyte_models/results/predictors/motif_expression_og.png', p)
}



heatmap_motivs <- function(hep_and_promoter_df){
# does not work
df <- df_matched
    hep_hm <- df %>%
        filter(count<4)%>%                    # H: why 4? 
        mutate(descr = paste(gene, var, sep = "."))  %>%  #combine gene and variability class
        select(-gene, -variability) %>%  
        column_to_rownames("descr") %>% 
        mutate(across(.cols = -any_of("age"), .fns = as.numeric)) %>% # H added  #### TODO: issue
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
    pdf(file = paste0(out_folder,"/predictors/Heatmap_predictors_clusters_all.pdf"), width = 5, height = 4)

    Heatmap(hep_hm,
            name = " ",
            column_order = colnames(hep_hm),
            col = custom_colors,
            column_split = col_annotation, 
            show_column_names = FALSE,
            row_order = c("TATA_box", "Inr", "GC_box", "CCAAT_box", 'count'))
    dev.off()

}
heatmap_motivs(df_matched)
colnames(df_matched)
head(df_matched)
