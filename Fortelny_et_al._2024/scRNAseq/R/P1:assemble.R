
library(readxl)
library(data.table)
library(dplyr)
library(tibble)
library(writexl)



gene_expression_data <- fread("droplet/full_results_10_22.tsv", header = TRUE, sep = "\t")
annotation_info <- read_excel("droplet/manual_annotation.xlsx", sheet = 1)
#last_column <- annotation_info %>% select((ncol(annotation_info)-2):ncol(annotation_info))

print(colnames(gene_expression_data))
print(head(gene_expression_data))



GSE204735_RNA.read.counts.csv