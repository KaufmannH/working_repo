
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(readr)
library(tidyr)
library(Seurat)

# check what HVG/LVG percentages there are across age in spleen B cells


b_cells <- readRDS("data_prep_for_classification/immune_models/data/B_cells_all_ages.rds")

tagged_cells <- b_cells@meta.data |>
    mutate(variability = case_when(
      SCT_snn_res.1.3 < 1 ~ "LVG",                               
      SCT_snn_res.1.3 > 5 ~ "HVG",           
      TRUE ~ "Other" ))
head(tagged_cells)
# is SCT_snn_res.1.3 res_var? no bootstrapping

# get expression data
mat <- GetAssayData(b_cells, assay = "SCT", slot = "data") # not scaleed bc not real expression
gene_df_wide <- as.data.frame(mat)
gene_df <- gene_df_wide |>
    rownames_to_column("gene") |>
    pivot_longer(cols = -gene, 
                    names_to =  "cell_id",
                    values_to = "sct_gene_expression")

# fuse data frames
meta_df <- tagged_cells |> rownames_to_column("cell_id")
merged_df <- dplyr::left_join(meta_df, gene_df, by = c("cell_id" = "cell_id")) 

final_df <- merged_df |>
  select(gene, age, cell_id, sct_gene_expression, variability) |>
  arrange(gene, cell_id)
head(meta_df)
head(gene_df)
head(final_df)


ratio_cells <- final_df |>
  group_by(age, variability) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(age) |>
  mutate(prop = n / sum(n))

plot <- ggplot(ratio_cells, aes(x = age, y = prop, fill = variability)) +
  geom_col(position = "stack") +
  labs(x = "Age", y = "Proportion of cells", fill = "Variability class") +
  theme_classic()

ggsave("data_prep_for_classification/immune_models/plots/spleen_var_comparison.png", plot)

# this is not the proper data, i am sure. 