# preparing a matrix with predictors to be merged with HVG information (gene-level)

setwd("~/Documents/Regulated_Noise/TMS_downstream/classification/")

library(tidyverse)

gene_description <- read_tsv("epd_mm/gene_description.txt", col_names = c("promoter", "gene", "description"))
promoter_motifs <- read_tsv("epd_mm/promoter_motifs.txt") %>%
  rename(promoter = `#Name`, TATA_box = `TATA-box`, CCAAT_box = `CCAAT-box`, GC_box = `GC-box`)
promoter_expression <- read_tsv("epd_mm/promoter_expression.txt", col_names = c("promoter", "unknown_col1", "unknown_col2", "expression_info"))

predictors <- promoter_motifs %>% left_join(gene_description, by = c("promoter" = "promoter")) %>%
  relocate(gene, .after = promoter)

write_tsv(predictors, file = "predictors/promoter_elements.tsv")
