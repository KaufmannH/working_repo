#' Prepare balanced HVG/LVG gene sets matched by mean expression
#'
#' Reads the output table of the RegNoise pipeline and prepares two balanced
#' gene sets (highly variable vs stable genes) with matched median gmean values.
#' Returns intermediate data frames and plots as a list object, if specified.
#'
#' @param input_file Path to the full_results_10_22.tsv file.
#' @param hvg_perc_hvg Minimum perc.hvg threshold to classify a gene as HVG.
#' @param lvg_res_var Maximum residual variance to classify a gene as LVG.
#' @param hvg_cutoff Minimum hvg_fraction to classify as potentially HVG.
#' @param lvg_cutoff_low Minimum lvg_fraction to classify as stable.
#' @param lvg_cutoff_high Maximum lvg_fraction to classify as potentially HVG.
#' @param return_plots_and_tables Logical; if TRUE, returns plots and tables in the output list.
#' @param seed Random seed for reproducible matching.
#' @param tissue_filter A character string of the tissue that should be filtered for.
#'
#' @return A list with elements:
#'   - gene_masterlist: tibble of all genes with classification
#'   - fraction_data: per-gene HVG/LVG/unclear fractions
#'   - prediction_data_matched: balanced HVG/LVG dataset
#'   - plots: list of ggplot objects (fraction, cutoff, gmean)
#' @export
#'
prepare_matched_gene_sets <- function(
    input_file,
    hvg_perc_hvg = 0.9,
    lvg_res_var = 1,
    hvg_cutoff = 0.01,
    lvg_cutoff_low = 0.5,
    lvg_cutoff_high = 0.9,
    return_plots_and_tables = TRUE,
    seed = 1,
    tissue_filter = NULL
) {
  if (!requireNamespace("tidyverse", quietly = TRUE)) stop("Package 'tidyverse' is required.")
  library(tidyverse)

  # 1. Read data
  hvg <- read_tsv(input_file, show_col_types = FALSE)

  # 2. Create masterlist with classification
  gene_masterlist <- hvg %>%
    mutate(
      hvg_stable = ifelse(hvg, perc.hvg >= hvg_perc_hvg, FALSE),
      lvg = res_var <= lvg_res_var,
      var_class = case_when(
        hvg_stable ~ "hvg",
        lvg ~ "lvg",
        TRUE ~ "unclear"
      ),
      var_class = factor(var_class, levels = c("lvg", "hvg", "unclear"))
    ) %>%
    select(gene, cluster, tissue, age, res_var, gmean, var_class)

  # 3. Optionally filter for a specific tissue
  if (!is.null(tissue_filter)) {
    message("Filtering to tissue: ", tissue_filter)
    gene_masterlist <- gene_masterlist %>%
      filter(tissue == tissue_filter)
  }

  # 4. Count class occurrences
  n_clusters_total <- gene_masterlist %>%
    distinct(tissue, age, cluster) %>%
    nrow()

  fraction_data <- gene_masterlist %>%
    group_by(gene) %>%
    summarize(
      n_hvg = sum(var_class == "hvg", na.rm = TRUE),
      n_lvg = sum(var_class == "lvg", na.rm = TRUE),
      n_unclear = sum(var_class == "unclear", na.rm = TRUE),
      n_detected = n(),
      .groups = "drop"
    ) %>%
    mutate(
      n_total = n_clusters_total,
      hvg_fraction = n_hvg / n_total,
      lvg_fraction = n_lvg / n_total,
      n_lvg_hvg = n_hvg + n_lvg
    )

  # 5. Apply cutoffs
  fraction_data <- fraction_data %>%
    mutate(
      class = case_when(
        hvg_fraction >= hvg_cutoff & lvg_fraction <= lvg_cutoff_high ~ "potentially_hvg",
        hvg_fraction <= hvg_cutoff & lvg_fraction >= lvg_cutoff_low  ~ "stable",
        TRUE ~ "unclear"
      )
    )

  # 6. Compute median gmean per gene and prepare prediction data
  median_gmean <- gene_masterlist %>%
    filter(gene %in% fraction_data$gene) %>%
    group_by(gene) %>%
    summarize(median_gmean = median(gmean, na.rm = TRUE), .groups = "drop")

  prediction_data <- fraction_data %>%
    left_join(median_gmean, by = "gene") %>%
    mutate(class = factor(class, levels = c("stable", "potentially_hvg", "unclear")))

  # 7. Match HVGs and LVGs by gmean
  set.seed(seed)

  hvg_sel <- prediction_data %>% filter(class == "potentially_hvg")
  lvg_sel <- prediction_data %>% filter(class == "stable")

  if (nrow(hvg_sel) == 0 | nrow(lvg_sel) == 0) {
    stop("No HVG or LVG genes found with the given thresholds or tissue filter.")
  }

  lvg_matches <- vector("list", nrow(hvg_sel))  # preallocate list for efficiency

  for (i in seq_len(nrow(hvg_sel))) {
    value <- hvg_sel$median_gmean[i]

    # Find the closest LVG match by median_gmean
    lvg_match <- lvg_sel %>%
      filter(abs(median_gmean - value) == min(abs(median_gmean - value))) %>%
      slice_sample(n = 1)

    # Remove the selected match from the pool
    lvg_sel <- lvg_sel %>% filter(!gene %in% lvg_match$gene)

    lvg_matches[[i]] <- lvg_match
  }

  lvg_matches <- bind_rows(lvg_matches)

  prediction_data_matched <- bind_rows(hvg_sel, lvg_matches) %>%
    mutate(class = factor(class, levels = c("stable", "potentially_hvg")))

  # 8. Optionally generate plots and full return
  if (return_plots_and_tables) {

    p_fraction <- ggplot(fraction_data, aes(x = lvg_fraction, y = hvg_fraction)) +
      geom_bin2d() +
      scale_y_log10() + scale_x_log10() +
      theme_classic() +
      labs(
        x = "Fraction of stable gene (log10)",
        y = "Fraction of HVG (log10)"
      )

    p_cutoff <- p_fraction +
      geom_vline(xintercept = c(lvg_cutoff_low, lvg_cutoff_high)) +
      geom_hline(yintercept = hvg_cutoff)

    p_gmean <- ggplot(prediction_data_matched,
                      aes(x = class, y = median_gmean, fill = class)) +
      geom_violin(scale = "width") +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      scale_y_log10() +
      scale_fill_brewer(palette = "Accent") +
      theme_classic() +
      theme(legend.position = "none")

    return(list(
      gene_masterlist = gene_masterlist,
      fraction_data = fraction_data,
      prediction_data_matched = prediction_data_matched,
      plots = list(
        fraction = p_fraction,
        cutoff = p_cutoff,
        gmean = p_gmean
      )
    ))

  } else {
    # Lightweight mode: only return matched dataset
    return(prediction_data_matched)
  }
}



#' Merge, filter, and process predictors
#'
#' Wrapper function that runs a full workflow:
#' 1. Convert human genes to mouse (optional, as indicated by convert_human)
#' 2. Filter by genes of interest
#' 3. Filter by standard deviation
#' 4. Merge predictors
#' 5. Filter by correlation
#'
#' @param pred_list Named list of tibbles (predictor tables)
#' @param goi Optional character vector of genes of interest
#' @param sd_cutoff Numeric, threshold for SD filtering
#' @param correlation_threshold Numeric, threshold for correlation filtering
#' @param convert_human Optional named logical vector (names have to match pred_list names) indicating which tables to convert from human to mouse gene symbols
#' @param ortholog_file Path to human-mouse ortholog TSV
#' @param n_cores Integer, cores for SD and correlation filtering
#' @param block_size Integer, block size for correlation computation
#' @return List with merged predictors and optionally correlation-filtered predictors
#' @export
#'
merge_filter_predictors <- function(pred_list,
                                    goi = NULL,
                                    sd_cutoff = 0.01,
                                    correlation_threshold = 0.9,
                                    convert_human = NULL,
                                    ortholog_file = NULL,
                                    n_cores = 4,
                                    block_size = 50) {

  # Step 0: check convert_human names
  if (!is.null(convert_human)) {
    if (!all(names(convert_human) %in% names(pred_list))) {
      stop("All names in convert_human must exist in pred_list.")
    }
    if (is.null(ortholog_file)) {
      warning("convert_human provided but ortholog_file is NULL. No conversion applied.")
    }
  }

  # Step 1: Convert human genes if requested
  if (!is.null(convert_human) && !is.null(ortholog_file)) {
    for (name in names(convert_human)) {
      if (convert_human[[name]]) {
        if (!name %in% names(pred_list)) {
          stop("Dataset '", name, "' not found in pred_list")
        }
        pred_list[[name]] <- convert_human_mouse(pred_list[[name]], ortholog_file)
      }
    }
  }

  # Step 2: Filter by genes of interest
  if (!is.null(goi)) {
    pred_list <- filter_goi(pred_list, goi)
  }

  # Step 3: Filter by standard deviation
  pred_list <- filter_sd(pred_list, sd_cutoff = sd_cutoff, n_cores = n_cores)

  # Step 4: Merge all predictors
  merged_pred <- merge_predictors(pred_list)

  # Step 5: Correlation filtering
  high_corr <- compute_correlation_matrix(
    merged_pred,
    correlation_threshold = correlation_threshold,
    n_cores = n_cores,
    block_size = block_size
  )

  # Filter merged predictors based on high correlation
  filtered_pred <- merged_pred |> dplyr::select(-dplyr::all_of(high_corr))

  return(filtered_pred)
}
