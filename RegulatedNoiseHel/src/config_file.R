# Select Tabula Muris tissue
# Select from
#datasets <- unique(sapply(list.files(so_folder), function(x) sub("_[^_]+$", "", x), USE.NAMES = FALSE))

SelectedTissues <- "Spleen"

testrun <- FALSE
n_cells <- 500 # only used if testrun == TRUE

min_cells <- 100 # minimum number of cells per subgroup/condition (in case of TMS age/sex group)
min_n_genes <- 500
min_n_counts <- 2500
max_n_counts <- 40000
max_percent_mt <- 5
npc_pca <- 100
dim_umap <- 30
dim_neighbors <- 30
cluster_res <- 1.3 # clustering resolution
num.de_cutoff <- 10 # number of violating genes in doublet detection
quantile_perc <- c(.025, .975) # bootstrapping, CI for ResVar
cycles_bootstrap <- 1000
hvg_cutoff <- 5
