#' Load and combine JASPAR transcription factor binding prediction data
#'
#' This function loads JASPAR transcription factor binding prediction TSV files from a directory,
#' renames columns with a prefix containing the intervals around promoter regions, and combines
#' all files by the 'gene' column. Optionally, a table of intervals can be provided to filter
#' which files are loaded.
#'
#' @param dir_path Character. Path to the directory containing TSV files.
#' @param intervals_table Optional data.frame or tibble with columns `interval1` and `interval2`.
#' Only files matching these intervals will be loaded. Default is NULL, which loads all summary files.
#' @return A single tibble with all summary files combined by the 'gene' column.
#' @examples
#' \dontrun{
#' combined_df <- load_jaspar("/path/to/tsv")
#' filtered_df <- load_jaspar("/path/to/tsv", intervals_table = data.frame(interval1 = 10000, interval2 = 5000))
#' }
#' @importFrom readr read_tsv
#' @importFrom stringr str_match
#' @importFrom dplyr rename_with full_join
#' @importFrom purrr reduce
#' @export
load_jaspar <- function(dir_path, intervals_table = NULL) {
  # List summary TSV files only
  files <- list.files(dir_path, pattern = "_summary\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop("No TSV files found in directory.")

  # Filter files if intervals_table is provided
  if (!is.null(intervals_table)) {
    files <- files[sapply(files, function(f) {
      numbers <- stringr::str_match(basename(f), "familial_profiles_(\\d+)_(\\d+)_summary\\.tsv")[, 2:3]
      any(apply(intervals_table, 1, function(row) all(as.numeric(numbers) == as.numeric(row)))) # check if numbers is matching any row in intervals_table
    })]
  }

  # Load, rename, and combine all files by 'gene'
  combined_df <- purrr::reduce(lapply(files, function(f) {
    numbers <- stringr::str_match(basename(f), "familial_profiles_(\\d+)_(\\d+)_summary\\.tsv")[, 2:3]
    prefix <- paste0("JASPAR_", numbers[1], "_", numbers[2], "_")
    df <- readr::read_tsv(f, show_col_types = FALSE)
    df <- dplyr::select(df, -dplyr::matches("index"))
    dplyr::rename_with(df, ~ paste0(prefix, .), .cols = -gene)
  }), dplyr::full_join, by = "gene")

  return(combined_df)
}


#' Load and combine ChIP transcription factor binding data in Hepatocytes
#'
#' This function loads ChIP transcription factor binding TSV files from a directory,
#' optionally filters files by intervals, and combines all files by the 'gene' column.
#' Assumes interval information is already present in the column names.
#'
#' @param dir_path Character. Path to the directory containing TSV files.
#' @param intervals_table Optional data.frame or tibble with columns `interval1` and `interval2`.
#' Only files containing these intervals (in their file names) will be loaded. Default is NULL (load all files).
#' @return A single `tibble` with all files combined by the 'gene' column.
#' @examples
#' \dontrun{
#' combined_dt <- load_hepatocyte_chip("/path/to/chip")
#' filtered_dt <- load_hepatocyte_chip("/path/to/chip",
#'   intervals_table = data.frame(interval1 = 250, interval2 = 50))
#' }
#' @import data.table
#' @importFrom stringr str_match
#' @importFrom tibble as_tibble
#' @export
load_hepatocyte_chip <- function(dir_path, intervals_table = NULL) {
  # List TSV files (e.g., promoters_250up_50down_wide.tsv)
  files <- list.files(dir_path, pattern = "promoters_\\d+up_\\d+down_wide\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop("No TSV files found in directory.")

  # Filter files by intervals in filename if intervals_table is provided
  if (!is.null(intervals_table)) {
    files <- files[sapply(files, function(f) {
      # Extract numbers from filename
      numbers <- stringr::str_match(basename(f), "promoters_(\\d+)up_(\\d+)down_wide\\.tsv")[, 2:3]
      any(apply(intervals_table, 1, function(row) all(as.numeric(numbers) == as.numeric(row))))
    })]
    if (length(files) == 0) stop("No files matched the provided intervals.")
  }

  # Read all files as data.tables
  dt_list <- lapply(files, data.table::fread)

  # Merge all tables by 'gene'
  combined_dt <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), dt_list)

  # convert to tibble
  combined_dt <- tibble::as_tibble(combined_dt)

  return(combined_dt)
}


#' Load and average promoter element data by gene
#'
#' This function reads a TSV file containing promoter element information from
#' the eukaryotic promoter database (EPD) and averages the numeric columns
#' (TATA_box, Inr, CCAAT_box, GC_box) by gene.
#'
#' @param dir_path Character. Path to the EPD promoter data directory.
#' @return A tibble with one row per gene, containing the averaged promoter element values.
#' @examples
#' \dontrun{
#' promoter_elements <- load_promoter_elements("/path/to/epd_mm")
#' }
#' @importFrom readr read_tsv
#' @importFrom dplyr group_by summarize select left_join
#' @export
load_promoter_elements <- function(dir_path) {

  # Read gene descriptions, keep only promoter and gene
  genes <- readr::read_tsv(file.path(dir_path, "gene_description.txt"),
                           col_names = c("promoter", "gene", "description"),
                           show_col_types = FALSE) |>
    dplyr::select(promoter, gene)

  # Read promoter motifs and rename columns
  promoter_motifs <- readr::read_tsv(file.path(dir_path, "promoter_motifs.txt"),
                                     show_col_types = FALSE) |>
    dplyr::rename(promoter = `#Name`,
                  TATA_box = `TATA-box`,
                  CCAAT_box = `CCAAT-box`,
                  GC_box = `GC-box`)

  # Join promoter motifs with genes
  df <- promoter_motifs |>
    dplyr::left_join(genes, by = "promoter")

  # Average by gene (if more than one promoter per gene)
  df <- df |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      TATA_box = mean(TATA_box),
      Inr = mean(Inr),
      CCAAT_box = mean(CCAAT_box),
      GC_box = mean(GC_box)
    )

  return(df)
}


#' Convert human gene symbols to mouse orthologs
#'
#' This function converts human gene symbols in a data frame to their corresponding
#' mouse orthologs using a provided ortholog mapping file. Only 1:1 ortholog pairs
#' (unique human–mouse mappings) are retained.
#'
#' @param df A predictor data frame containing a column named \code{gene} with human gene symbols.
#' @param ortholog_file Character. Path to a TSV file with two columns: human and mouse gene symbols.
#' The file should have no header; column order is assumed to be \code{human}, \code{mouse}.
#'
#' @return A tibble with the same columns as the input, but with human genes replaced by
#' their corresponding mouse orthologs in the \code{gene} column. Rows without a valid
#' 1:1 mouse ortholog are removed.
#'
#' @examples
#' \dontrun{
#' mouse_df <- convert_human_mouse(human_df, "/path/to/orthologs_human_mouse.tsv")
#' }
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr left_join rename filter select
#' @export
convert_human_mouse <- function(df, ortholog_file) {
  if (!"gene" %in% names(df)) {
    stop("Input data frame must contain a column named 'gene'.")
  }

  # Prepare single match orthologs mouse - human
  h2m <- readr::read_tsv(ortholog_file, col_names = c("human", "mouse"), show_col_types = FALSE)
  h2m_single_match <- h2m |>
    dplyr::filter(duplicated(human) == FALSE, duplicated(mouse) == FALSE)
  message("Ortholog map loaded: ", nrow(h2m_single_match), " 1:1 mappings available.")

  # convert gene symbols
  df <- dplyr::left_join(df, h2m_single_match, by = c("gene" = "human")) |>
    dplyr::rename(Gene_symbol_human = gene) |>
    dplyr::rename(gene = mouse) |>
    dplyr::filter(!is.na(gene)) |>
    dplyr::select(gene, everything(), -Gene_symbol_human)

  return(df)
}


#' Filter a list of predictor data frames by genes of interest
#'
#' This function takes a list of tibbles (each containing a \code{gene} column)
#' and filters each one to include only the specified genes of interest.
#'
#' @param pred_list A list of tibbles or data frames, each containing a column named \code{gene}.
#' @param goi A character vector of genes of interest to retain.
#'
#' @return A list of filtered tibbles with only the genes of interest.
#'
#' @examples
#' \dontrun{
#' filtered_list <- filter_goi(pred_list, c("GeneA", "GeneB", "GeneC"))
#' }
#'
#' @importFrom dplyr filter
filter_goi <- function(pred_list, goi) {

  pred_list <- lapply(pred_list, function(df) {
    if (!"gene" %in% names(df)) {
      stop("Each tibble must contain a 'gene' column.")
    }
    dplyr::filter(df, gene %in% goi)
  })

  return(pred_list)
}


#' Filter predictor tables by standard deviation
#'
#' This function computes the standard deviation of all numeric columns in each
#' tibble of predictors within a list and removes columns whose standard deviation
#' is below a given cutoff.
#'
#' @param pred_list A list of tibbles or data frames containing numeric predictor columns.
#' @param sd_cutoff Numeric. Minimum standard deviation required to keep a column.
#' @param n_cores Integer. Number of cores to use for parallel computation (default = 4).
#'
#' @return A list of filtered tibbles with low-variance columns removed.
#'
#' @examples
#' \dontrun{
#' filtered_list <- filter_sd(pred_list, sd_cutoff = 0.01, n_cores = 4)
#' }
#'
#' @importFrom dplyr select
#' @importFrom tidyselect where all_of
#' @importFrom parallel mclapply
#' @export
filter_sd <- function(pred_list, sd_cutoff = 0.01, n_cores = 4) {

  pred_list <- lapply(pred_list, function(df) {
    # select numeric columns
    num_cols <- df |> dplyr::select(tidyselect::where(is.numeric))
    # compute standard deviation
    sd_vector <- parallel::mclapply(num_cols, sd, na.rm = TRUE, mc.cores = n_cores) |>
      unlist()
    # filter
    low_sd_cols <- names(sd_vector[sd_vector < sd_cutoff])
    if (length(low_sd_cols) > 0) {
      message("Removing ", length(low_sd_cols), " columns with SD < ", sd_cutoff, ". ")
    }
    keep_cols <- setdiff(names(df), low_sd_cols)
    df <- df |> dplyr::select(tidyselect::all_of(keep_cols))

    return(df)
  })

  return(pred_list)
}


#' Merge predictor tables by gene
#'
#' This function merges multiple predictor tables (tibbles) by the column \code{gene}.
#' Missing predictor values for genes that are not present in all tables are filled with \code{0}.
#'
#' @param pred_list A list of tibbles or data frames, each containing a column named \code{gene}.
#'
#' @return A tibble combining all predictors by gene, with missing values replaced by 0.
#'
#' @examples
#' \dontrun{
#' merged <- merge_predictors(pred_list)
#' }
#'
#' @importFrom dplyr full_join mutate across
#' @importFrom tidyselect where
#' @importFrom tidyr replace_na
#' @export
merge_predictors <- function(pred_list) {
  if (length(pred_list) == 0) stop("pred_list is empty.")

  # Ensure all have a 'gene' column
  if (!all(vapply(pred_list, function(df) "gene" %in% names(df), logical(1)))) {
    stop("All predictor tables must contain a 'gene' column.")
  }

  # Merge all tables by 'gene'
  merged_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = "gene"), pred_list)

  # Replace NAs with 0 for all numeric columns
  merged_df <- merged_df |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ tidyr::replace_na(., 0)))

  return(merged_df)
}


#' Compute a Correlation Matrix in Parallel Blocks
#'
#' Efficiently computes a full correlation matrix for all numeric columns
#' in a data frame or tibble, using parallel block-wise computation.
#'
#' This function divides the numeric columns into blocks and computes correlations
#' in parallel to handle wide datasets efficiently. The final output is assembled
#' into a complete correlation matrix.
#'
#' @param predictors A data frame, tibble, or data.table containing predictor variables.
#'   Only numeric columns are included in the correlation computation.
#' @param n_cores Integer. Number of CPU cores to use for parallel computation.
#'   On macOS, using multiple cores may not speed up computation due to BLAS threading.
#' @param block_size Integer. Number of columns to include per processing block.
#' @param correlation_threshold Numeric. Predictors with correlations above this value will be filtered.
#' @param return_corr_matrix Boolean. Returns a list of high_corr + cor_matrix if TRUE.
#'
#' @return Column names of predictors that are highly correlated, for filtering.
#' Optional (return_corr_matrix): a list object with the standard return and
#' a numeric correlation matrix of size \code{p × p}, where \code{p}
#'   is the number of numeric columns in \code{predictors}.
#'
#' @details
#' Internally, this function:
#' \enumerate{
#'   \item Selects numeric columns from the input.
#'   \item Converts them to a matrix for fast computation.
#'   \item Splits columns into blocks and computes correlations in parallel using \code{parallel::mclapply()}.
#'   \item Combines all block results into a full symmetric correlation matrix.
#' }
#'
#' @examples
#' set.seed(42)
#' df <- tibble(
#'   x = rnorm(100),
#'   y = rnorm(100),
#'   z = rnorm(100),
#'   group = sample(letters[1:3], 100, replace = TRUE)
#' )
#'
#' cor_mat <- compute_correlation_matrix(df, n_cores = 2, block_size = 2)
#' print(dim(cor_mat))
#'
#' @seealso [parallel::mclapply()], [stats::cor()]
#' @export
compute_correlation_matrix <- function(predictors, correlation_threshold = 0.9,
                                       n_cores = 4, block_size = 50, return_corr_matrix = FALSE) {
  # Extract numeric columns into matrix
  num_mat <- predictors |>
    dplyr::select(tidyselect::where(is.numeric)) |>
    as.matrix()

  if (ncol(num_mat) < 2) {
    stop("Need at least two numeric columns to compute correlations.")
  }

  # Define blocks for parallel computation
  blocks <- split(seq_len(ncol(num_mat)), ceiling(1:ncol(num_mat)/block_size))

  # Function to compute one block
  cor_block <- function(cols) {
    cor(num_mat[, cols, drop = FALSE], num_mat, use = "pairwise.complete.obs")
  }

  # Compute blocks in parallel
  cor_list <- parallel::mclapply(blocks, FUN = cor_block, mc.cores = n_cores)

  # Combine block-wise output to generate correlation matrix
  cor_matrix <- do.call(rbind, cor_list)

  # identify columns with high correlation
  high_corr <- caret::findCorrelation(cor_matrix, cutoff = correlation_threshold, names = TRUE)

  if (return_corr_matrix) {
    res <- list("high_corr" = high_corr,
                "cor_matrix" = cor_matrix)
    return(res)

  } else {
    return(high_corr)
  }
}



