# Tag macrophages/monocytes - using marker genes

library(Seurat)
library(SeuratData)
library(ggplot2)


palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")


pbmc <- readRDS("Helene/data/pmbc_seurat_object.rds")


#  average expression of key markers
# conversions: CD16 = FCGR3A, HLA-DR is heterodimer: HLA-DRB1 + HLA-DRA, CD11c = ITGAX
macrophage_markers <- c("CD68", "CD14", "CSF1R")
classical_monocyte_markers <- c("CD14", "CCR2", "CCR5", "CD62L")
intermediate_monocyte_markers <- c("CD14", "FCGR3A", "HLA-DRA","HLA-DRB1", "CD68", "ITGAX")
nonclassical_monocyte_markers <- c("CD14", "FCGR3A", "HLA-DRA","HLA-DRB1", "CX3CR1")
all_markers <- c("CD14", "CCR2", "CCR5", "CD62L", "FCGR3A", "HLA-DRA","HLA-DRB1", "CD68", "ITGAX",  "CX3CR1")


# macros
macro_plot <- FeaturePlot(pbmc, features = macrophage_markers)
ggsave("Helene/plots/macrophage_markers.png", plot = macro_plot, width = 8, height = 6, dpi = 300)

# classical monos

classical_mono_plot <- DotPlot(pbmc, features = classical_monocyte_markers, group.by = "seurat_clusters") +
  scale_size(range = c(0.5, 6)) + # Adjust dot size range
  scale_color_gradient(low = "#D6EFC3", high = "#123230") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs( color = "Expression", size = "Percent Expressed") +
    theme_classic()
ggsave("Helene/plots/classical_monocyte_markers.png", plot = classical_mono_plot, width = 8, height = 6, dpi = 300)


# intermediate monos
intermediate_mono_plot <- DotPlot(pbmc, features = intermediate_monocyte_markers, group.by = "seurat_clusters") +
  scale_size(range = c(0.5, 6)) + # Adjust dot size range
  scale_color_gradient(low = "#D6EFC3", high = "#123230") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs( color = "Expression", size = "Percent Expressed") +
    theme_classic()
ggsave("Helene/plots/interm_monocyte_markers.png", plot = intermediate_mono_plot, width = 8, height = 6, dpi = 300)


# nonclassical monos
noncalssical_mono_plot <- DotPlot(pbmc, features = nonclassical_monocyte_markers, group.by = "seurat_clusters") +
  scale_size(range = c(0.5, 6)) + # Adjust dot size range
  scale_color_gradient(low = "#D6EFC3", high = "#123230") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs( color = "Expression", size = "Percent Expressed") +
  theme_classic()
ggsave("Helene/plots/noncalssical_monocyte_markers.png", plot =noncalssical_mono_plot , width = 8, height = 6, dpi = 300)


# all
dot_plot <- DotPlot(pbmc, features = all_markers, group.by = "seurat_clusters") +
  scale_size(range = c(0.5, 6)) + # Adjust dot size range
  scale_color_gradient(low = "#D6EFC3", high = "#123230") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs( color = "Expression", size = "Percent Expressed") +
  theme_classic()

ggsave("Helene/plots/all_monocyte_markers.png", plot = dot_plot , width = 8, height = 6, dpi = 300)




# tag data set
# average expression is low but there is a clear cluster hwere some express it so this is the monocyte cluster

old_clusters <- c("0",  "1",  "2",  "3",  "4",  "5",  "6" , "7",  "8",  "9",  "10")
new_clusters <- c("0",  "1" , "2" , "3",  "Monocytes",  "5",  "6" , "7" , "8" , "Monocytes" , "Monocytes")

Idents(pbmc) <- factor(Idents(pbmc), levels = old_clusters, labels = new_clusters)

umap <- DimPlot(pbmc, reduction = "umap", label = TRUE) + labs(color = "Cluster") +
  xlab("UMAP 1") +
  ylab("UMAP 2")
ggsave("Helene/plots/umap_monocytes.png", plot = umap, width = 8, height = 6, dpi = 300)


# FeaturePlot split by stimType

umap_cond <- DimPlot(pbmc, reduction = "umap", group.by = "Condition") +
  theme_classic() +
   xlab("UMAP 1") +
  ylab("UMAP 2")
  ggsave("Helene/plots/umap_condition.png", plot = umap_cond, width = 8, height = 6, dpi = 300)


# From the dotplots wit looks like it is mostly cluster 4 but in the stimulations it is vivible that in cluster 9 and 10 there are the rest of the conditions


SaveSeuratRds(pbmc, "Helene/data/pmbc_seurat_object_ann.rds")

