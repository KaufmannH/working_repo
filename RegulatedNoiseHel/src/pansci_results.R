

library('tidyverse')
library('Seurat')

# 1. load resvar and bootstrapping files and merge

path <- "data/processed/pansci_liver/res_var_tables"

files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
df_3 <- readRDS(files[1])
df_6 <- readRDS(files[2])

df_3$source <- "df_3"
df_6$source <- "df_6"
df_new <- dplyr::bind_rows(df1, df2)

so_3 <- readRDS('data/processed/pansci_liver/seurat_objects/Liver_03_months_Female_processed.rds')

# save table





# 2. check result

# residual variance distribution

bin_width <- 0.05
safety_log <- 1e-6

new_log <- log10(df_new$res_var + safety_log)
log_brks <- seq(floor(min(new_log, na.rm = TRUE)), ceiling(max( new_log, na.rm = TRUE)), by = bin_width)


plot_df_new <- df_new |>
    filter(gmean != 0) |>
    mutate(log_res_var = log10(res_var + 1e-6)) |>
    filter(is.finite(log_res_var)) |> 
    #mutate(bin = cut(res_var, breaks = seq(0, ceiling(max(df_new$res_var, df_new$res_var, na.rm = TRUE) / bin_width) * bin_width, by = bin_width))) |>
        mutate(bin = cut(log_res_var, breaks = log_brks, include.lowest = TRUE)) |>
    # count(bin, name = "n_genes")
    group_by(bin) |>
    summarize(n_genes = n())

plot_new <- ggplot(plot_df_new, aes(x = bin, y = n_genes)) +
    geom_col(fill = '#5195D3') +
    labs(x = "Residual variance bins (log10)",
        y = "Gene count", title = "New") +
theme_classic() +
theme(
    strip.background      = element_blank(),
    axis.text.x           = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y           = element_text(size = 20),
    axis.title.x      = element_text(size = 20),
    axis.title.y      = element_text(size = 20)) +
scale_x_discrete(breaks = levels(plot_df_new$bin)[seq(1, length(levels(plot_df_new$bin)), by = 10)])


ggsave(paste0("data/comparison/pansci/res_var_dist.png"),  plot_new, width = 15, height = 5)





