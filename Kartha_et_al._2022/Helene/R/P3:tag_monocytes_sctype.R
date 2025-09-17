# Tag macrophages/monocytes - scType 

library(Seurat)
library(ggplot2)
library(devtools)
library(openxlsx)
library(HGNChelper)
library(tidyverse)
library(writexl)


palette <- c("#70916B", "#A4C089", "#D6EFC3", "#F5F5F5", "#C7EAE5" ,"#5AB4AC", 
"#123230",  "#556967",  "#B8C7C7",  "#E0D5C8",  "#D0BFAD",  "#AE9982", "#BA613E", "#843915", "#471339", "#3B052C")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)


# load object
pbmc <- readRDS("Helene/data/pmbc_seurat_object.rds")
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(pbmc[["RNA"]])));
# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(pbmc[["RNA"]]$scale.data) else as.matrix(pbmc[["RNA"]]@scale.data)
# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
# now a few clusters have the same score: took random one bc i only care about monocytes
sctype_scores <- cL_resutls %>% 
    group_by(cluster) %>% 
    top_n(n = 1, wt = scores)  #%>% # mutate(type = ifelse(type =='Non-classical monocytes', 'Non-classical monocytes', 'other')) %>%
    #distinct()
sctype_scores

# set low-confident (low ScType score) clusters to "unknown"
#sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
#print(sctype_scores[,1:3])

# add just sctype
pbmc[['sctype']] <-  pbmc@meta.data %>%
    left_join(sctype_scores, by= c('seurat_clusters' = 'cluster')) %>%
    pull(type)

# add all: needed for checking scores
#meta <-  pbmc@meta.data %>%
 #   left_join(sctype_scores, by= c('seurat_clusters' = 'cluster')) 
#pbmc@meta.data <- meta
#head(pbmc@meta.data)



Idents(pbmc) <- 'sctype'

umap_sctype <- DimPlot(pbmc, reduction = "umap") +
  scale_color_manual(values = palette) +
  theme_classic() +
  xlab("UMAP 1") +
  ylab("UMAP 2")
ggsave("Helene/plots/umap_sctype.png", plot = umap_sctype , width = 8, height = 6, dpi = 300)

check_score <-  pbmc@meta.data %>%
    arrange(descending(scores.x))

# print csv to check scores
#write_xlsx(check_score, './data/cell_type_sctype.xlsx')

#TODO: check score of second best


SaveSeuratRds(pbmc, "Helene/data/pmbc_seurat_object_sctype.rds")

