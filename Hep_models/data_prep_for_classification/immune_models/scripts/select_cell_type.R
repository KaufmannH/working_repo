
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(readr)
library(tidyr)
library(Seurat)

# here i check how many cells of each cell type there are in 10X for spleen and bone marrow. 

# marrow
marrow_1 <- readRDS('data_prep_for_classification/immune_models/data/Marrow_1m_male_processed.rds')
marrow_3 <- readRDS('data_prep_for_classification/immune_models/data/Marrow_3m_female_processed.rds')
marrow_18 <- readRDS('data_prep_for_classification/immune_models/data/Marrow_18m_female_processed.rds')
marrow_18_m <- readRDS ('data_prep_for_classification/immune_models/data/Marrow_18m_male_processed.rds')
marrow_21 <- readRDS ('data_prep_for_classification/immune_models/data/Marrow_21m_female_processed.rds')
marrow_24 <- readRDS('data_prep_for_classification/immune_models/data/Marrow_24m_male_processed.rds')
marrow_30 <- readRDS('data_prep_for_classification/immune_models/data/Marrow_30m_male_processed.rds')


marrow_merged <- merge(
  x = marrow_1,
  y = list(marrow_3, marrow_18, marrow_18_m, marrow_21, marrow_24, marrow_30),
  project = "Marrow")
head(marrow_merged)
summary_marrow <- marrow_merged@meta.data 
table(summary_marrow$age)

comp_marrow <- marrow_merged@meta.data |>
  count(age, cell_ontology_class) |>
  group_by(age) |>
  mutate(proportion = n / sum(n)) |>
  ungroup()

plot <- ggplot(comp_marrow, aes(x = age, y = proportion, fill = cell_ontology_class)) +
  geom_col() +
  labs(x = "Age", y = "Proportion", fill = "Cell type") +
  theme_classic()

ggsave("data_prep_for_classification/immune_models/plots/marrow_comparison.png", plot)

# -> granulocytes are in marrow all ages at enough numbers -> suitable cell type


# spleen
spleen_1 <- readRDS('data_prep_for_classification/immune_models/data/Spleen_1m_male_processed.rds')
spleen_3 <- readRDS('data_prep_for_classification/immune_models/data/Spleen_3m_male_processed.rds')
spleen_3m <- readRDS('data_prep_for_classification/immune_models/data/Spleen_3m_female_processed.rds')
spleen_18 <- readRDS('data_prep_for_classification/immune_models/data/Spleen_18m_female_processed.rds')
spleen_21 <- readRDS('data_prep_for_classification/immune_models/data/Spleen_21m_female_processed.rds')
spleen_24 <- readRDS('data_prep_for_classification/immune_models/data/Spleen_24m_male_processed.rds')
spleen_30 <- readRDS('data_prep_for_classification/immune_models/data/Spleen_30m_male_processed.rds')


spleen_merged <- merge(
  x = spleen_1,
  y = list(spleen_3, spleen_3m, spleen_18, spleen_21, spleen_24, spleen_30),
  project = "Spleen")
head(spleen_merged)
summary_spleen <- spleen_merged@meta.data 
table(summary_spleen$age)

comp_spleen <- spleen_merged@meta.data |>
  count(age, cell_ontology_class) |>
  group_by(age) |>
  mutate(proportion = n / sum(n)) |>
  ungroup()

plot <- ggplot(comp_spleen, aes(x = age, y = proportion, fill = cell_ontology_class)) +
  geom_col() +
  labs(x = "Age", y = "Proportion", fill = "Cell type") +
  theme_classic()

ggsave("data_prep_for_classification/immune_models/plots/spleen_comparison.png", plot)

# -> B cells are in spleen all ages at enough numbers -> suitable cell type


spleen_filtered <-  subset(spleen_merged, subset = cell_ontology_class == "B cell")
marrow_filtered <-  subset(marrow_merged, subset = cell_ontology_class == "granulocyte")

table(spleen_filtered$age)
table(marrow_filtered$age)

# they have the exact same number of granulocytes??

saveRDS(marrow_filtered, "data_prep_for_classification/immune_models/data/granulocytes_all_ages.rds")
saveRDS(spleen_filtered, "data_prep_for_classification/immune_models/data/B_cells_all_ages.rds")


