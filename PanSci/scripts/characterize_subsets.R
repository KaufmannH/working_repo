
library('Seurat')
library('dplyr')
library('ggplot2')


# check og objects

path <- "data/subsets"
files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)

df <- data.frame(
  file = character(),
  n_cells = numeric(),
  stringsAsFactors = FALSE)

for (f in files) {
  message("Reading: ", f)
  so <- readRDS(f)
  df <- rbind(df, data.frame(
    file = basename(f),
    n_cells = ncol(so)))
  rm(so)
  gc()
}
df$file <- factor(df$file, levels = df$file[order(df$n_cells)])

options(scipen = 999)
plot <- ggplot(df, aes(x = file, y = n_cells)) +
  geom_col(fill = "lightblue") +

  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold")) +
  labs(
    x = "File",
    y = "Number of cells",
    title = "Cells per Seurat object")

ggsave("plots/cell_counts_per_object.png", plot, width = 10, height = 6, dpi = 150)




# check only WT mice

path <- "data/prepped_for_pipeline"
files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)

df <- data.frame(
  file = character(),
  n_cells = numeric(),
  stringsAsFactors = FALSE)

for (f in files) {
  message("Reading: ", f)
  so <- readRDS(f)
  df <- rbind(df, data.frame(
    file = basename(f),
    n_cells = ncol(so)))
  rm(so)
  gc()
}
df$file <- factor(df$file, levels = df$file[order(df$n_cells)])

options(scipen = 999)
plot <- ggplot(df, aes(x = file, y = n_cells)) +
  geom_col(fill = "lightblue") +

  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold")) +
  labs(
    x = "File",
    y = "Number of cells",
    title = "Cells per Seurat object (WT only)")

ggsave("plots/cell_counts_per_object_WT.png", plot, width = 10, height = 6, dpi = 150)


# plot the resulting cells after the regnoise pipeline



# check only WT mice

path <- 'RegulatedNoiseHel/data/processed/pansci_liver/seurat_objects/'
files <- list.files(path, full.names = TRUE)


df <- data.frame(
  file = character(),
  n_cells = numeric(),
  stringsAsFactors = FALSE)

for (f in files) {
  message("Reading: ", f)
  so <- readRDS(f)
  df <- rbind(df, data.frame(
    file = basename(f),
    n_cells = ncol(so)))
  rm(so)
  gc()
}
df$file <- factor(df$file, levels = df$file[order(df$n_cells)])

options(scipen = 999)
plot <- ggplot(df, aes(x = file, y = n_cells)) +
  geom_col(fill = "lightblue") +

  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold")) +
  labs(
    x = "File",
    y = "Number of cells",
    title = "Cells per Seurat object (post pipeline)")

ggsave("PanSci/plots/cell_counts_per_object_post_pipeline.png", plot, width = 10, height = 6, dpi = 150)

