#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title check_size
#' @description
#' Checks if input Seurat object is small enough to run clustOpt
#' @param input object to check
#' @param verbose print warning about size
#' @return NULL
#'
#' @export
#'
check_size <- function(input, verbose = TRUE) {
  if (!methods::is(input, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (ncol(input) >= 2E5) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' leverage_sketch: Uses leverage score based sampling to reduce the size of the
#' input Seurat object and create a sketch assay using this subsample
#'
#' @param input Seurat object
#' @param sketch_size Number of cells to include in the sketch assay
#' @param verbose print messages
#' @return Seurat object with sketch assay
#'
#' @export
#'
leverage_sketch <- function(input, sketch_size, verbose = TRUE) {
  if (is.null(sketch_size)) {
    if (verbose) {
      message("No sketch_size specified, defaulting to 10% of cells")
    }
    sketch_size <- ncol(input) * 0.1
  }
  
  input <- Seurat::NormalizeData(input)
  input <- Seurat::FindVariableFeatures(input)
  input <- Seurat::SketchData(
    object = input,
    ncells = sketch_size,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )
  Seurat::DefaultAssay(input) <- "sketch"
  # Return only the sketch assay
  input <- Seurat::DietSeurat(input, assays = "sketch")
}

#' @title get_valid_samples_by_group
#' @description
#' Checks if each group in a specified metadata variable contains at least 3 samples
#' with 100 or more cells, returning valid sample names if the condition is met.
#' @param input A Seurat object containing metadata
#' @param subject_ids The name of the metadata column containing sample IDs
#' @param group_variable The name of the metadata variable defining groups to be checked
#' @return A vector of sample names meeting the criteria, or NULL if any group fails to meet the criteria
#' @details
#' This function inspects each group in the specified metadata column to ensure they have
#' at least 3 samples with 100 or more cells. If any group does not satisfy this requirement,
#' the function returns NULL and issues a warning listing those groups.
#' @examples
#' # Assuming `seurat_obj` is a Seurat object with metadata columns "sample_id" and "group_var":
#' get_valid_samples_by_group(seurat_obj, "sample_id", "group_var")
#' @export
#'
get_valid_samples_by_group <- function(input, subject_ids, group_variable) {
  
  # Check if the group_variable exists in the metadata
  if (!group_variable %in% colnames(input@meta.data)) {
    stop(paste("The variable", group_variable, "is not found in the metadata."))  # Error if group_variable is missing
  }
  
  # Summarize the number of cells per sample per group
  sample_summary <- input@meta.data %>%
    group_by(.data[[group_variable]], !!sym(subject_ids)) %>%  # Group by the specified group variable and subject IDs
    summarize(cell_count = n(), .groups = "drop") %>%
    ungroup()
  
  # Filter samples to include only those with at least 100 cells
  sufficient_samples <- sample_summary %>%
    filter(cell_count >= 100) %>%  # Filter for samples meeting the cell count criterion
    group_by(.data[[group_variable]]) %>%
    summarize(num_samples = n(), .groups = "drop")  # Summarize the count of sufficient samples per group
  
  # Find groups that do not meet the sample count requirement
  insufficient_groups <- sufficient_samples %>%
    filter(num_samples < 3) %>%  # Identify groups with fewer than 3 sufficient samples
    pull(.data[[group_variable]])
  
  # Check for insufficient groups
  if (length(insufficient_groups) > 0) {
    warning("The following groups have fewer than 3 samples with >= 100 cells each: ", 
            paste(insufficient_groups, collapse = ", "))  # Warning for groups that fail the criterion
    return(NULL)  # Return NULL if any group fails the sample count condition
  }
  
  # Retrieve sample names for groups that meet the requirements
  valid_samples <- sample_summary %>%
    filter(cell_count >= 100) %>%  # Filter to include samples with sufficient cell counts
    pull(!!sym(subject_ids))  # Extract sample names that meet the criteria
  
  # Return valid sample names with confirmation message
  message("All groups have sufficient samples. Returning the sample names that meet the criteria.")
  return(valid_samples)  # Return the list of valid samples
}

#' @title sil_summary
#' @description 
#' Summarises the silhouette score distribution output by clustOpt
#' @param input output of clustOpt
#' @return A data.frame summarising the silhouette score distribution
#'
#' @export
#'
sil_summary <- function(input) {
  input |>
    dplyr::group_by(resolution) |>
    dplyr::summarize(
      mean_score = base::mean(avg_width, na.rm = TRUE),
      variance_score = stats::var(avg_width, na.rm = TRUE),
      standard_error_score = stats::sd(avg_width,
                                      na.rm = TRUE) / sqrt(length(avg_width)),
      cluster_mean_score = base::mean(cluster_avg_widths, na.rm = TRUE),
      cluster_variance_score = stats::var(cluster_avg_widths, na.rm = TRUE),
      cluster_standard_error_score = stats::sd(cluster_avg_widths, na.rm = TRUE) /
        sqrt(length(cluster_avg_widths))
    )
}
