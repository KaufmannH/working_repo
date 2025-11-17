library(ComplexHeatmap)

# Exploring the predictors prediction quality
resvar_plots <- list()
gmean_plots <- list()

out_folder <- "exploratory_figures"
if(!dir.exists(out_folder)) {dir.create(out_folder)}

# add titles and fix axis labels
for (pred in c("TATA_box", "Inr", "CCAAT_box", "GC_box")) {
  resvar_plots[[pred]] <- ggplot(gene_masterlist, aes(x = as.factor(round(eval(parse(text = pred)), 2)), y = res_var)) +
    geom_boxplot(outlier.size = 0.1) +
    ylim(0, 6) +
    facet_wrap(~ tissue) +
    theme_classic()
  
  gmean_plots[[pred]] <- ggplot(gene_masterlist, aes(x = as.factor(round(eval(parse(text = pred)), 2)), y = gmean)) +
    geom_boxplot(outlier.size = 0.1) +
    ylim(0, 6) +
    facet_wrap(~ tissue) +
    theme_classic()
  
  ggsave(resvar_plots[[pred]], filename = file.path(out_folder, paste0(pred, "_resvar_pred.png")), width = 15, height = 10)
  ggsave(gmean_plots[[pred]], filename = file.path(out_folder, paste0(pred, "_gmean_pred.png")), width = 15, height = 10)
  
}

# Heatmap of gene_masterlist (percentages of elements in LVG/HVG per cluster)
gml_sum <- gene_masterlist %>%
  filter(!is.na(TATA_box) & !is.na(Inr) & !is.na(CCAAT_box) & !is.na(GC_box)) %>%
  group_by(tissue, age, cluster, var_class) %>%
  summarize(TATA = mean(TATA_box),
            Inr = mean(Inr),
            CCAAT = mean(CCAAT_box),
            GC = mean(GC_box)) %>%
  ungroup()

# overview split by var_class
hm <- gml_sum %>%
  mutate(descr = paste(tissue, age, cluster, var_class, sep = ".")) %>%
  select(-tissue, -age, -cluster, -var_class) %>%
  column_to_rownames("descr") %>%
  as.matrix() %>%
  t()

col_annotation <- word(colnames(hm), sep = "\\.", 4)
col_annotation[col_annotation == "hvg"] <- "HVG"
col_annotation[col_annotation == "lvg_s"] <- "stable"
col_annotation[col_annotation == "intermediate"] <- "intermediate"

pdf(file = file.path(out_folder, "Heatmap_predictors_clusters.pdf"),
     width = 5,
     height = 4)

Heatmap(hm,
        name = " ",
        column_order = colnames(hm), 
        column_split = col_annotation, 
        show_column_names = FALSE,
        row_order = c("TATA", "Inr", "GC", "CCAAT"))

dev.off()

# only hvg split by tissue
hm_hvg <- gml_sum %>%
  filter(var_class %in% c("hvg")) %>%
  mutate(descr = paste(tissue, age, cluster, var_class, sep = ".")) %>%
  select(-tissue, -age, -cluster, -var_class) %>%
  column_to_rownames("descr")

row_annotation <- paste(word(rownames(hm_hvg), sep = "\\.", 1), word(rownames(hm_hvg), sep = "\\.", 4))
Heatmap(hm_hvg, row_split = row_annotation, row_title_rot = 0, show_row_names = FALSE)

# only lvg, split by tissue
hm_lvg <- gml_sum %>%
  filter(var_class %in% c("lvg_s")) %>%
  mutate(descr = paste(tissue, age, cluster, var_class, sep = ".")) %>%
  select(-tissue, -age, -cluster, -var_class) %>%
  column_to_rownames("descr")

row_annotation <- paste(word(rownames(hm_lvg), sep = "\\.", 1), word(rownames(hm_lvg), sep = "\\.", 4))
Heatmap(hm_lvg, row_split = row_annotation, row_title_rot = 0, show_row_names = FALSE)

# plot %age of features for HVGs and LVGs
perc_scatter <- list()
perc_box <- list()
for (feature in c("TATA", "Inr", "CCAAT", "GC")) {
  plot_data <- gml_sum %>% 
    filter(var_class %in% c("lvg_s", "hvg")) %>%
    select(tissue, age, cluster, var_class, feature) %>%
    pivot_wider(names_from = var_class, values_from = feature) %>%
    filter(!is.na(lvg_s) & !is.na(hvg))
  
  perc_scatter[[feature]] <- ggplot(plot_data, aes(x = hvg, y = lvg_s)) +
    geom_point(size = 0.5) +
    geom_abline() +
    xlim(0, 1) +
    ylim(0, 1) +
    xlab("HVG") +
    ylab("stable") +
    ggtitle(feature) +
    theme_classic()
  
  perc_box[[feature]] <- ggplot(gml_sum, aes(x = factor(var_class, levels = c("hvg", "lvg_s", "intermediate")), y = !!rlang::sym(feature))) +
    geom_boxplot() +
    ggtitle(feature) +
    xlab("") +
    ylab("") +
    scale_x_discrete(limits = c("hvg", "intermediate", "lvg_s"), labels = c("hvg" = "HVG", "intermediate" = "intermediate",
                              "lvg_s" = "stable")) +
    theme_classic()
  ggsave(filename = file.path(out_folder, paste0(feature, "box_per_varclass.png")), perc_box[[feature]], width = 3, height = 3)
}

pw_scatter <- wrap_plots(perc_scatter, ncol = 2)
pw_scatter
pw_box <- wrap_plots(perc_box, ncol = 2)
pw_box
ggsave(filename = file.path(out_folder, "box_pw_perc_predictors.png"), pw_box, width = 6, height = 6)

# plot percentage values for all genes, grouped by element and var_class
indiv_genes <- list()
for (feature in c("TATA_box", "Inr", "CCAAT_box", "GC_box")) {
  indiv_genes[[feature]] <- ggplot(gene_masterlist, aes(x = var_class, y = !!rlang::sym(feature))) +
    geom_boxplot() +
    ggtitle(feature) +
    theme_classic()
}

pw_indiv_genes <- wrap_plots(indiv_genes, ncol = 2)
pw_indiv_genes

vln_gml <- list()
features <- c("TATA_box", "Inr", "CCAAT_box", "GC_box")
for (feature in features) {
  
  if (feature == features[length(features)]) {
    vln_gml[[feature]] <- ggplot(gene_masterlist, aes(x = factor(var_class, levels = c("hvg", "lvg_s", "intermediate")), y = !!rlang::sym(feature), fill = var_class)) +
      geom_violin(scale = "width") + 
      theme_classic() + 
      ylab(feature) +
      theme(legend.position = "none", axis.title = element_blank(), 
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"),  
            axis.text.y = element_text(size = 12, colour = "black"),
            axis.title.y = element_text(size = 16, colour = "black", angle = 0, vjust = 0.5)) +
      scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
      scale_fill_brewer(palette = "Dark2") # for custom coloring
  } else {
    vln_gml[[feature]] <- ggplot(gene_masterlist, aes(x = factor(var_class, levels = c("hvg", "lvg_s", "intermediate")), y = !!rlang::sym(feature), fill = var_class)) +
      geom_violin(scale = "width") + 
      theme_classic() + 
      ylab(feature) +
      theme(legend.position = "none", axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 12, colour = "black"),
            axis.title.y = element_text(size = 16, colour = "black", angle = 0, vjust = 0.5),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
      scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
      scale_fill_brewer(palette = "Dark2") # for custom coloring
  }
  
}

pw_vln <- wrap_plots(vln_gml, ncol = 1)
pw_vln
ggsave(filename = file.path(out_folder, "vln_stack_perc_predictors.pdf"), pw_vln, width = 4, height = 5)


# probability of HVGs and LVGs to be HVGs/LVGs in other clusters
gene_class_matrix <- gene_masterlist %>% 
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

a <- ggplot(data, aes(x = lvg_fraction, y = hvg_fraction)) +
  geom_bin2d() +
  scale_y_log10() +
  scale_x_log10() +
  theme_classic() +
  ylab("Fraction of HVG (log10)") +
  xlab("Fraction of stable gene (log10)")


b <- ggplot(data, aes(hvg_fraction)) +
  geom_histogram() +
  theme_classic() +
  ylab("number of genes") +
  xlab("Fraction of clusters in which a given gene is HVG")

c <- ggplot(data, aes(lvg_fraction)) +
  geom_histogram() +
  theme_classic() +
  ylab("number of genes") +
  xlab("Fraction of clusters in which a given gene is stable")

abc <- a / (b | c)

ggsave(filename = file.path(out_folder, "hvg_lvg_fractions_density.pdf"), a, width = 4, height = 3)
ggsave(filename = file.path(out_folder, "hvg_fraction_hist.pdf"), b, width = 4, height = 3)
ggsave(filename = file.path(out_folder, "lvg_fraction_hist.pdf"), c, width = 4, height = 3)
ggsave(filename = file.path(out_folder, "fractions_panel.pdf"), abc, width = 4, height = 3)




