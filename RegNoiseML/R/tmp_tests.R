
# system.time(data <- prepare_matched_gene_sets(input_file = "../data/full_results_10_22.tsv"))
# data <- data[["prediction_data_matched"]]
# gc()
# system.time(pred_list <- list(jaspar = load_jaspar(dir_path = "../classification/JASPAR_new/results_v2"),
#                   hep_chip = load_hepatocyte_chip(dir_path = "../hepatocyte_models/chip_out"),
#                   promoter_elements = load_promoter_elements(dir_path = "../classification/epd_mm")))
#
# convert_vec <- c("jaspar" = F, "hep_chip" = T, "promoter_elements" = F)
#
# pred_list_subset <- lapply(pred_list, function(df) df[, c("gene", head(names(df)[names(df) != "gene"], 200))])
# system.time(preds_final <- merge_filter_predictors(pred_list_subset,
#                                        goi = matched$gene,
#                                        convert_human = convert_vec,
#                                        ortholog_file = "~/Documents/celltype_annotation/scType/ortho_maps/human.txt",
#                                        n_cores = 8,
#                                        block_size = 50))






# system.time(
#   jas <- load_jaspar(dir_path = "../classification/JASPAR_new/results_v2")
# )
#
# system.time(
#   hep_chip <- load_hepatocyte_chip(dir_path = "../hepatocyte_models/chip_out")
# )
# system.time(
#   jas2 <- load_jaspar_2(dir_path = "../classification/JASPAR_new/results_v2")
# )
# system.time(
#   prom <- load_promoter_elements(dir_path = "../classification/epd_mm")
# )
#
#
# ortholog_file = "~/Documents/celltype_annotation/scType/ortho_maps/human.txt"
#
#
# data <- prepare_matched_gene_sets(input_file = "../data/full_results_10_22.tsv")
# matched <- data[["prediction_data_matched"]]











