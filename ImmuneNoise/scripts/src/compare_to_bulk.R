# load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr) 
library(purrr)
library(tibble)
library(writexl)


# this is for the doc, i am comparing the bulk gene expression to the one in the TMS

strat_df <- readRDS("droplet/data/strat_df.rds")
# lavin et al. 2019
bulk_df <- read_xlsx("comparison/data/rna_bulk_expression_bins.xlsx") # recalucalte after fix
# Li et al. 2019
bulk_brain_df <- read_xlsx("comparison/data/rna_bulk_expression_bins.xlsx") # fix issue

bulk_annotated <- bulk_df |>
  mutate(tissue = str_to_sentence(tissue)) |>
    mutate(method = "bulk") |>
    select(gene, tissue, expression_bin, method)
unique(bulk_annotated$tissue)

strat_lean <- strat_df |>
    select(gene, tissue, expression_bin) |>
    mutate(method = "single cell") 
unique(strat_lean$tissue)


matched_genes <-  strat_lean |>
    full_join(bulk_annotated, by = c("gene", "tissue", "method", "expression_bin")) |>
    arrange(gene, tissue) |> 
    group_by(gene, tissue) |>
    filter(all(c("bulk", "single cell") %in% method)) |>
    ungroup()
  
 summary <- matched_genes |>
    group_by(gene, tissue, method) |>
    summarise(mean_bin = mean(as.integer(expression_bin))) |>
    mutate(rounded_mean = floor(mean_bin + 0.5))
summary


matched_stingent <- summary |>
  select(gene, tissue, method, rounded_mean) |>
  pivot_wider(names_from = method, values_from = rounded_mean) |>
  filter(!is.na(`bulk`) & !is.na(`single cell`)) |>
  mutate(same_bin = `bulk` == `single cell`) |>
  group_by(tissue) |>
  summarise(
    n_genes = n(),
    n_same_bin = sum(same_bin),
    percent_same = round(100 * n_same_bin / n_genes, 1))
matched_stingent


matched_loose <- summary |>
  select(gene, tissue, method, rounded_mean) |>
  pivot_wider(names_from = method, values_from = rounded_mean) |>
  filter(!is.na(`bulk`) & !is.na(`single cell`)) |>
  mutate(same_bin = abs(`bulk` - `single cell`) <= 1) |>
   group_by(tissue) |>
   summarise(
    n_genes = n(),
    n_same_bin = sum(same_bin),
    percent_same = round(100 * n_same_bin / n_genes, 1))
matched_loose


matched_even_looser <- summary |>
  select(gene, tissue, method, rounded_mean) |>
  pivot_wider(names_from = method, values_from = rounded_mean) |>
  filter(!is.na(`bulk`) & !is.na(`single cell`)) |>
  mutate(same_bin = abs(`bulk` - `single cell`) <= 2) |>
   group_by(tissue) |>
   summarise(
    n_genes = n(),
    n_same_bin = sum(same_bin),
    percent_same = round(100 * n_same_bin / n_genes, 1))
matched_even_looser


# what genes are the ones that dont change? 

matched_loose_genes <- summary |>
  select(gene, tissue, method, rounded_mean) |>
  pivot_wider(names_from = method, values_from = rounded_mean) |>
  filter(!is.na(`bulk`) & !is.na(`single cell`)) |>
  mutate(same_bin = abs(`bulk` - `single cell`) <= 1) |>
  filter(same_bin == TRUE) |>
  arrange(gene, tissue)
matched_loose_genes

