
library(ggplot2)
library(patchwork)


# 1. exploring gene overlap
# 2. make correlation of lvg and hvg percentage across clusters


bound_df <- read_tsv('hep_all_tech_df.tsv')
head(bound_df)

## 1. explore gene overlap between facs, droplet and easysci

# how many clusters per technology? 
clusters_unique <- bound_df |>
  distinct(cluster_id, technology) |>
  group_by(technology) |>
  count()
clusters_unique


# how big is the overlap of detected genes between technologies?
gene_overlap_df <- bound_df |>
   distinct(gene, technology) |>
   add_count(gene, name = "n_tech") |>
    group_by(technology) |>
    summarise(
    total_unique_genes   = n(),                    
    overlapping_genes    = sum(n_tech == 3),        
    tech_specific_genes  = sum(n_tech == 1),      
    .groups = "drop")
gene_overlap_df


# how many are lost in the ensembl conversion? / which genes are in droplet/facs but not in pansci? 
genes_not_in_pansci <- bound_df |>
  distinct(gene, technology) |>
  mutate(present = TRUE) |>
  pivot_wider(
    names_from  = technology,
    values_from = present,
    values_fill = FALSE) |>
  filter(
    facs == TRUE,
    droplet == TRUE,
    easysci == FALSE) |>
    distinct() |>
    pull(gene)
genes_not_in_pansci
# there are 1798 genes that are not in pansci but in both facs and droplet

ensembl_genes = sum(str_detect(genes_not_in_pansci, "^ENS"))
ensembl_genes
# there are 0 genes with retained ensembl id that are in facs and droplet but not pansci

# reverse check: which genes with ensembl ids are shared? 
genes_shared <- bound_df |>
  distinct(gene, technology) |>
  mutate(present = TRUE) |>
  pivot_wider(
    names_from  = technology,
    values_from = present,
    values_fill = FALSE) |>
  filter( facs,  droplet,  easysci ) |>
    distinct() |>
    pull(gene)
genes_shared

shared_ensembl_genes = sum(str_detect(genes_shared, "^ENS"))
shared_ensembl_genes
# none are

# how many of the genes that are in the pansci but not in the droplet or facs have ensembl ids?
genes_pansci <- bound_df |>
  distinct(gene, technology) |>
  mutate(present = TRUE) |>
  pivot_wider(
    names_from  = technology,
    values_from = present,
    values_fill = FALSE) |>
  filter( facs == FALSE,  
          droplet == FALSE,  
          easysci == TRUE ) |>
    distinct() |>
    pull(gene)
genes_pansci

pansci_ensembl_genes = sum(str_detect(genes_pansci, "^ENS"))
pansci_ensembl_genes
# 181 of the genes unique to pansci are ensembl id genes

# are these all ensembl genes? 
n_ens <- sum(str_detect(bound_df, "^ENS")) 
n_ens

unique_bound <- bound_df |>  distinct(gene) |> pull()
n_ens_unique <- sum(str_detect(unique_bound, "^ENS")) 
# yes in the whole df there are also only 181 genes with ensembl ids





## 2. correlation 
#  how many times a gene is hvg/lvg in each technology?

# get genes that are shared between all technologies
genes_shared_list <- bound_df |>
  distinct(gene, technology) |>
  mutate(present = TRUE) |>
  pivot_wider(
    names_from  = technology,
    values_from = present,
    values_fill = FALSE) |>
  filter(
    facs == TRUE,
    droplet == TRUE,
    easysci == TRUE) |>
    pull(gene)

# calc percentages of gene labelled hvg/lvg over all clusters in each technology
gene_cluster_stats <- bound_df |>
  filter(gene %in% genes_shared_list) |>
  mutate(stable_hvg = perc_hvg > 0.9) |>
  distinct(gene, technology, cluster_id, hvg, lvg, stable_hvg) |>
  group_by(gene, technology) |>
  summarise(
    n_clusters_total = n_distinct(cluster_id),
    n_detected = n_distinct(cluster_id),  
    n_hvg      = sum(hvg, na.rm = TRUE),
    n_stable_hvg =  sum(stable_hvg, na.rm = TRUE),
    n_lvg      = sum(lvg, na.rm = TRUE),
    pct_detected = round(100 * n_detected / n_clusters_total, 2),
    pct_hvg      = round(100 * n_hvg / n_clusters_total, 2),
    pct_stable_hvg      = round(100 * n_stable_hvg  / n_clusters_total, 2),
    pct_lvg      = round(100 * n_lvg / n_clusters_total, 2),
    .groups = "drop")

gene_cluster_stats





# correlation of percentage of clusters as HVG

hvg_wide <- gene_cluster_stats |>
  select(gene, technology, pct_hvg) |>
  pivot_wider(
    names_from  = technology, 
    values_from = pct_hvg)

# spearman correlations
cor_fd <- cor(hvg_wide$facs,    hvg_wide$droplet,
              use = "complete.obs", method = "spearman")
cor_fe <- cor(hvg_wide$facs,    hvg_wide$easysci,
              use = "complete.obs", method = "spearman")
cor_de <- cor(hvg_wide$droplet, hvg_wide$easysci,
              use = "complete.obs", method = "spearman")


p_fd <- ggplot(hvg_wide, aes(x = facs, y = droplet)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#1b9e77") +
  geom_smooth(method = "lm", se = FALSE, colour = "#1b9e77") +
  labs(
    x = "FACS: % clusters where gene is HVG",
    y = "Droplet: % clusters where gene is HVG",
    title = paste0("FACS vs Droplet (Spearman r = ",
                   round(cor_fd, 2), ")")) +
  theme_classic()

p_fe <- ggplot(hvg_wide, aes(x = facs, y = easysci)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#d95f02") +
  geom_smooth(method = "lm", se = FALSE, colour = "#d95f02") +
  labs(
    x = "FACS: % clusters where gene is HVG",
    y = "EasySci: % clusters where gene is HVG",
    title = paste0("FACS vs EasySci (Spearman r = ",
                   round(cor_fe, 2), ")")) +
  theme_classic()

p_de <- ggplot(hvg_wide, aes(x = droplet, y = easysci)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#7570b3") +
  geom_smooth(method = "lm", se = FALSE, colour = "#7570b3") +
  labs(
    x = "Droplet: % clusters where gene is HVG",
    y = "EasySci: % clusters where gene is HVG",
    title = paste0("Droplet vs EasySci (Spearman r = ",
                   round(cor_de, 2), ")")) +
  theme_classic()


combined_plot <- p_fd | p_fe | p_de

ggsave(filename = "TMS/comparison/plots/liver/hvg_correlation.png", plot = combined_plot, width = 12, height = 4, dpi = 300)



# correlation of percentage of clusters as STABLE HVGs

hvg_wide <- gene_cluster_stats |>
  select(gene, technology, pct_stable_hvg) |>
  pivot_wider(
    names_from  = technology, 
    values_from = pct_stable_hvg)

# spearman correlations
cor_fd <- cor(hvg_wide$facs,    hvg_wide$droplet,
              use = "complete.obs", method = "spearman")
cor_fe <- cor(hvg_wide$facs,    hvg_wide$easysci,
              use = "complete.obs", method = "spearman")
cor_de <- cor(hvg_wide$droplet, hvg_wide$easysci,
              use = "complete.obs", method = "spearman")


p_fd <- ggplot(hvg_wide, aes(x = facs, y = droplet)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#1b9e77") +
  geom_smooth(method = "lm", se = FALSE, colour = "#1b9e77") +
  labs(
    x = "FACS: % clusters where gene is HVG",
    y = "Droplet: % clusters where gene is HVG",
    title = paste0("FACS vs Droplet (Spearman r = ",
                   round(cor_fd, 2), ")")) +
  theme_classic()

p_fe <- ggplot(hvg_wide, aes(x = facs, y = easysci)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#d95f02") +
  geom_smooth(method = "lm", se = FALSE, colour = "#d95f02") +
  labs(
    x = "FACS: % clusters where gene is HVG",
    y = "EasySci: % clusters where gene is HVG",
    title = paste0("FACS vs EasySci (Spearman r = ",
                   round(cor_fe, 2), ")")) +
  theme_classic()

p_de <- ggplot(hvg_wide, aes(x = droplet, y = easysci)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#7570b3") +
  geom_smooth(method = "lm", se = FALSE, colour = "#7570b3") +
  labs(
    x = "Droplet: % clusters where gene is HVG",
    y = "EasySci: % clusters where gene is HVG",
    title = paste0("Droplet vs EasySci (Spearman r = ",
                   round(cor_de, 2), ")")) +
  theme_classic()


combined_plot <- p_fd | p_fe | p_de

ggsave(filename = "TMS/comparison/plots/liver/hvg_stable_correlation.png", plot = combined_plot, width = 12, height = 4, dpi = 300)





# correlation of percentage of clusters as LVG

lvg_wide <- gene_cluster_stats |>
  select(gene, technology, pct_lvg) |>
  pivot_wider(
    names_from  = technology, 
    values_from = pct_lvg)

# spearman correlations
cor_fd <- cor(lvg_wide$facs,    lvg_wide$droplet,
              use = "complete.obs", method = "spearman")
cor_fe <- cor(lvg_wide$facs,    lvg_wide$easysci,
              use = "complete.obs", method = "spearman")
cor_de <- cor(lvg_wide$droplet, lvg_wide$easysci,
              use = "complete.obs", method = "spearman")


p_fd <- ggplot(lvg_wide, aes(x = facs, y = droplet)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#1b9e77") +
  geom_smooth(method = "lm", se = FALSE, colour = "#1b9e77") +
  labs(
    x = "FACS: % clusters where gene is LVG",
    y = "Droplet: % clusters where gene is LVG",
    title = paste0("FACS vs Droplet (Spearman r = ",
                   round(cor_fd, 2), ")")) +
  theme_classic()

p_fe <- ggplot(lvg_wide, aes(x = facs, y = easysci)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#d95f02") +
  geom_smooth(method = "lm", se = FALSE, colour = "#d95f02") +
  labs(
    x = "FACS: % clusters where gene is LVG",
    y = "EasySci: % clusters where gene is LVG",
    title = paste0("FACS vs EasySci (Spearman r = ",
                   round(cor_fe, 2), ")")) +
  theme_classic()

p_de <- ggplot(lvg_wide, aes(x = droplet, y = easysci)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#7570b3") +
  geom_smooth(method = "lm", se = FALSE, colour = "#7570b3") +
  labs(
    x = "Droplet: % clusters where gene is LVG",
    y = "EasySci: % clusters where gene is LVG",
    title = paste0("Droplet vs EasySci (Spearman r = ",
                   round(cor_de, 2), ")")) +
  theme_classic()


combined_plot <- p_fd | p_fe | p_de

ggsave(filename = "TMS/comparison/plots/liver/lvg_correlation.png", plot = combined_plot, width = 12, height = 4, dpi = 300)






# density plots for resvar distribution

tech_cols_fill <- c( facs = "#74C4E7", droplet  = "#F7D639", easysci  = "#E77E01")

plot <- ggplot(bound_df, aes(x = res_var, fill = technology, color = technology)) +
  geom_density(
    alpha = 0.5,
    linewidth = 0.5) +
  scale_x_log10() +
   scale_fill_manual(values = tech_cols_fill) +
     scale_colour_manual(values = tech_cols_fill) +
  theme_classic(base_size = 20) +
  guides(colour = "none") +  
  labs(
    x = "Residual variance (log10)",
    y = "Density",
    fill = "Technology" )

  ggsave('ImmuneNoise/comparison/plots/liver/resvar_dist_density.png', plot, width = 10, height = 5)




# wow you actaully used the script! LiGrü 🎃