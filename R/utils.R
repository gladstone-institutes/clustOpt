#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title check_size
#' @description
#' Checks if input Seurat object is small enough to run clustOpt
#' @param input object to check
#' @return NULL
#'
#' @export
check_size <- function(input) {
  if (!inherits(result_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (ncol(input) >= 2E5) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @title leverage_sketch
#' @description
#' Uses leverage score based sampling to reduce the size of the
#' input Seurat object and create a sketch assay using this subsample
#'
#' @param input Seurat object
#' @param sketch_size Number of cells to include in the sketch assay
#' @param verbose print messages
#' @return Seurat object with sketch assay
#'
#' @export
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures
#' @importFrom Seurat SketchData DefaultAssay DietSeurat
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

#' @title get_valid_samples
#' @description
#' Checks if there are at least 3 samples with a minimum number of cells,
#' returning valid sample names if the condition is met.
#' @param input A Seurat object containing metadata
#' @param subject_ids The name of the metadata column containing subject IDs
#' @param min_cells Minimum cells per subject
#' @return A vector of sample names meeting the criteria, or NULL if the
#' requirement is not met
#' @details
#' This function inspects the metadata to ensure there are at least 3 samples
#' with min_cells. If this condition is not satisfied, the function
#' returns NULL and issues a warning.
#' @examples
#' \dontrun{
#' # `seurat_obj` is a Seurat object with a metadata column "sample_id":
#' get_valid_samples(seurat_obj, "sample_id", 50)
#' }
#' @export
#'
#' @importFrom dplyr group_by summarize filter
get_valid_samples <- function(input, subject_ids, min_cells) {
  # Summarize the number of cells per sample
  sample_summary <- input@meta.data |>
    dplyr::group_by(!!sym(subject_ids)) |>
    dplyr::summarize(cell_count = n(), .groups = "drop") |>
    dplyr::ungroup()

  sufficient_samples <- sample_summary |>
    dplyr::filter(cell_count >= min_cells)

  if (nrow(sufficient_samples) < 3) {
    warning(
      "There are fewer than 3 samples with at least ",
      min_cells, " cells."
    )
    return(NULL)
  }
  # Retrieve subject names that meet the requirements
  valid_samples <- sufficient_samples |>
    dplyr::pull(!!sym(subject_ids))
  # Return the list of valid subject names with confirmation message
  message(
    "There are sufficient subjects. Returning the subjects names ",
    "have at least ", min_cells, " cells."
  )
  return(valid_samples)
}

#' @title sil_summary
#' @description
#' Summarises the silhouette score distribution output by clustOpt
#' @param input output of clustOpt
#' @return A data.frame summarising the silhouette score distribution
#'
#' @export
#' @importFrom dplyr group_by summarize
sil_summary <- function(input) {
  input |>
    dplyr::group_by(resolution) |>
    dplyr::summarize(
      median_score = stats::median(avg_width, na.rm = TRUE),
      variance_score = stats::var(avg_width, na.rm = TRUE),
      standard_error_score = stats::sd(avg_width,
        na.rm = TRUE
      ) / sqrt(length(avg_width)),
      cluster_median_score = stats::median(cluster_median_widths, na.rm = TRUE),
      cluster_variance_score = stats::var(cluster_median_widths, na.rm = TRUE),
      cluster_standard_error_score = stats::sd(cluster_median_widths,
        na.rm = TRUE
      ) /
        sqrt(length(cluster_median_widths))
    )
}
