setwd("/Users/stephan/Documents/Regulated_Noise/TMS_downstream/classification/")

library(tidyverse)
library(tidymodels)
library(themis)
library(ranger)
library(patchwork)
library(ComplexHeatmap)
source("pred.fct.R")

pred <- read_tsv(file = "predictors/promoter_elements.tsv")

# aggregate predictor information: if more than one promoter is annotated, create a masterannotation (if any annotated promoter has a given element for a given gene, it's positive, if none has it, it's negative)
pred <- pred %>% group_by(gene) %>% summarize(TATA_box = mean(TATA_box), Inr = mean(Inr), CCAAT_box = mean(CCAAT_box), GC_box = mean(GC_box))

hvg <- read_tsv(file = "../data/full_results_10_22.tsv")

# create redundant master list of genes, joined with predictor information
gene_masterlist <- hvg %>% 
  mutate(hvg_stable = ifelse(hvg, ifelse(perc.hvg >= 0.9, T, F), F)) %>%
  mutate(lvg_soft = ifelse(res_var <= 1, T, F)) %>%
  mutate(var_class = ifelse(hvg_stable, "hvg", ifelse(lvg_soft, "lvg_s", "intermediate"))) %>%
  mutate(var_class = factor(var_class, levels = c("lvg_s", "hvg", "intermediate"))) %>%
  select(gene, cluster, tissue, age, res_var, gmean, var_class) %>%
  left_join(pred, by = c("gene" = "gene"))

# create exploratory figures for the observations (gene list)
source("exploratory_figures.R")

genelist <- gene_masterlist

# probability of HVGs and LVGs to be HVGs/LVGs in other clusters
gene_class_matrix <- genelist %>% 
  mutate(cluster_id = paste(tissue, age, cluster, sep = ".")) %>%
  select(gene, cluster_id, var_class) %>%
  pivot_wider(names_from = cluster_id, values_from = var_class) %>%
  column_to_rownames("gene") %>%
  as.matrix()

saveRDS(gene_class_matrix, "gene_class_matrix.rds")

count = 0
for(gene in rownames(gene_class_matrix)) {
  vec <- gene_class_matrix[gene,]
  
  n_hvg <- length(which(vec == "hvg"))
  n_lvg <- length(which(vec == "lvg_s"))
  n_intermediate <- length(which(vec == "intermediate"))
  n_total <- length(vec)
  if(count == 0) {
    tib <- tibble(gene, n_hvg, n_lvg, n_intermediate, n_total)
  } else {
    tib <- rbind(tib, tibble(gene, n_hvg, n_lvg, n_intermediate, n_total))
  }
  count <- count + 1
  
}

data <- tib %>% mutate(hvg_fraction = n_hvg/n_total, lvg_fraction = n_lvg/n_total, n_lvg_hvg = n_lvg + n_hvg)
saveRDS(data, "hvg_lvg_fractions_across_clusters.rds")

dir.create("probability_cutoff_tests")

a <- ggplot(data, aes(x = lvg_fraction, y = hvg_fraction)) +
  geom_bin2d() +
  scale_y_log10() +
  scale_x_log10() +
  theme_classic() +
  ylab("Fraction of HVG (log10)") +
  xlab("Fraction of stable gene (log10)")

ggsave("probability_cutoff_tests/fraction_hvg_lvg_log.pdf", plot = a, width = 6, height = 5)

# Try out different cutoffs
potentially_hvg <- data %>% filter(hvg_fraction >= 0.01, lvg_fraction <= 0.9)
stable <- data %>% filter(hvg_fraction <= 0.01, lvg_fraction >= 0.5)

b <- ggplot(data, aes(x = lvg_fraction, y = hvg_fraction)) +
  geom_bin2d() +
  scale_y_log10() +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = c(0.5, 0.9)) +
  geom_hline(yintercept = 0.01) +
  ylab("Fraction of HVG (log10)") +
  xlab("Fraction of stable gene (log10)")

ggsave("probability_cutoff_tests/fraction_hvg_lvg_log_cutoffs.pdf", plot = b, width = 6, height = 5)

c <- ggplot(data, aes(x = lvg_fraction, y = hvg_fraction)) +
  geom_bin2d() +
  theme_classic() +
  geom_vline(xintercept = c(0.5, 0.9)) +
  geom_hline(yintercept = 0.01) +
  ylab("Fraction of HVG") +
  xlab("Fraction of stable gene")

ggsave("probability_cutoff_tests/fraction_hvg_lvg_cutoffs.pdf", plot = c, width = 6, height = 5)

data <- data %>% 
  mutate(class = ifelse(gene %in% potentially_hvg$gene, "potentially_hvg", 
                        ifelse(gene %in% stable$gene, "stable", "unclear")))

prediction_data <- data %>% left_join(pred, by = c("gene" = "gene")) %>%
  filter(!is.na(TATA_box), !is.na(Inr), !is.na(CCAAT_box), !is.na(GC_box))

n_class <- ggplot(prediction_data %>% mutate(class = factor(class, levels = c("potentially_hvg", "unclear", "stable"))), aes(x = class, fill = class)) +
  geom_bar(stat = "count") +
  scale_fill_brewer(palette = "Accent") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("probability_cutoff_tests/n_class.pdf", plot = n_class, width = 4, height = 5)

median_gmean <- gene_masterlist %>% 
  filter(gene %in% prediction_data$gene) %>% 
  group_by(gene) %>% 
  summarize(median_gmean = median(gmean))

prediction_data <- left_join(prediction_data, median_gmean, by = c("gene" = "gene"))

median_gmean_plt <- ggplot(prediction_data %>% mutate(class = factor(class, levels = c("potentially_hvg", "unclear", "stable"))), aes(x = class, y = median_gmean, fill = class)) +
  geom_violin(position = position_dodge(0.8), scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.1, outlier.shape = NA) +
  scale_fill_brewer(palette = "Accent") +
  scale_y_log10() +
  theme_classic() +
  theme(legend.position = "none") 

ggsave(filename = "probability_cutoff_tests/median_gmean.pdf", plot = median_gmean_plt, width = 4, height = 5)

ggplot(prediction_data %>% pivot_longer(cols = c(TATA_box, Inr, CCAAT_box, GC_box), names_to = "predictor", values_to = "predictor_value"),
       aes(x = class, y = predictor_value, fill = class)) +
  geom_violin(position = position_dodge(0.8), scale = "width") +
  #geom_boxplot(position = position_dodge(0.8), width = 0.05, outlier.shape = NA) +
  scale_fill_brewer(palette = "Accent") +
  facet_wrap(~predictor) +
  theme_classic()

# prediction_data %>% select(gene, TATA_box, Inr, CCAAT_box, GC_box) %>% 
#   column_to_rownames("gene") %>%
#   Heatmap(name = " ",
#           column_order = colnames(hm), 
# #          column_split = col_annotation, 
#           show_column_names = FALSE)
# #          row_order = c("TATA", "Inr", "GC", "CCAAT"))

hvg_sel <- prediction_data %>% filter(class == "potentially_hvg")
lvg_sel <- prediction_data %>% filter(class == "stable")
for (row in 1:nrow(hvg_sel)) {
  set.seed(1)
  value <- hvg_sel[row,]$median_gmean # set gmean value of hvg
  
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

prediction_data_matched <- rbind(hvg_sel, lvg_matches) %>%
  mutate(class = factor(class, levels = c("stable", "potentially_hvg")))

median_gmean_matched_plt <- ggplot(prediction_data_matched, aes(x = class, y = median_gmean, fill = class)) +
  geom_violin(position = position_dodge(0.8), scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.1, outlier.shape = NA) +
  scale_fill_brewer(palette = "Accent") +
  scale_y_log10() +
  theme_classic() +
  theme(legend.position = "none")

ggsave(filename = "probability_cutoff_tests/median_gmean_matched.pdf", plot = median_gmean_matched_plt, width = 4, height = 5)

n_class_matched <- ggplot(prediction_data_matched %>% mutate(class = factor(class, levels = c("potentially_hvg", "unclear", "stable"))), aes(x = class, fill = class)) +
  geom_bar(stat = "count") +
  scale_fill_brewer(palette = "Accent") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("probability_cutoff_tests/n_class_matched.pdf", plot = n_class_matched, width = 4, height = 5)