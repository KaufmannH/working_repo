# Are there monocytes with a proliferation signature?
# used marker genes for monocyte ID

library(Seurat)
library(ggplot2)

palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")

pbmc <- readRDS("Helene/data/pmbc_seurat_object_ann.rds")


# predefined cell cycle genes
cc.genes <- Seurat::cc.genes 
cc.genes$s.genes

# check proliferation/cell cycle genes
pbmc <- CellCycleScoring(
  object = pbmc,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE  
)

pbmc$Proliferating <- ifelse(pbmc$Phase %in% c("S", "G2M"), "Proliferating", "Non-Proliferating")

cluster_colors <- c("G1" = "#AE9982", "G2M" = "#A4C089", "S" =  "#471339")
umap_proliferation <- DimPlot(pbmc, reduction = "umap", split.by = "Proliferating") +
  scale_color_manual(values = cluster_colors) +
  theme_classic() +
  xlab("UMAP 1") +
  ylab("UMAP 2")
ggsave("Helene/plots/umap_proliferation.png", plot = umap_proliferation , width = 8, height = 6, dpi = 300)

# genes for proliferation extra
# Ki-67 = MKI67(nuclear division), PCNA (S), topoisomerase II alpha = TOP2A (replication)



# yes there are cells with a prolifaration signature and it is distributed across the clusters. 


