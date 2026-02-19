#' @importFrom dplyr group_by summarize filter n sym pull ungroup
#' @importFrom rlang .data
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Input Validation Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Normalize verbose parameter for backward compatibility
#'
#' Converts logical TRUE/FALSE to integer 1/0 and validates numeric input.
#'
#' @param verbose A logical or non-negative numeric value
#' @return An integer verbosity level
#' @keywords internal
normalize_verbose <- function(verbose) {
  if (length(verbose) != 1 || is.na(verbose)) {
    stop("verbose must be a logical or non-negative integer")
  }
  if (is.logical(verbose)) {
    return(as.integer(verbose))
  }
  if (!is.numeric(verbose)) {
    stop("verbose must be a logical or non-negative integer")
  }
  if (verbose < 0) {
    stop("verbose must be a non-negative integer")
  }
  as.integer(verbose)
}
#' @title check_size
#' @description
#' Checks if input Seurat object is small enough to run clustOpt
#' @param input object to check
#' @return NULL
#'
#' @keywords internal
check_size <- function(input) {
  if (!inherits(input, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (ncol(input) >= 2E5) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @title get_valid_samples
#' @description
#' Checks if there are at least 3 samples with a minimum number of cells,
#' returning valid sample names if the condition is met.
#' @param input A Seurat object containing metadata
#' @param subject_ids The name of the metadata column containing subject IDs
#' @param min_cells Minimum cells per subject
#' @param verbose Integer verbosity level (0 = silent, 1 = milestones,
#' 2 = detailed, 3 = includes Seurat output)
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
#' @importFrom dplyr group_by summarize filter n sym
get_valid_samples <- function(input, subject_ids, min_cells, verbose = 0) {
  # Summarize the number of cells per sample
  sample_summary <- input@meta.data |>
    dplyr::group_by(!!dplyr::sym(subject_ids)) |>
    dplyr::summarize(cell_count = dplyr::n(), .groups = "drop")

  # Split samples into sufficient and insufficient
  sufficient <- sample_summary$cell_count >= min_cells
  sufficient_samples <- sample_summary[sufficient, ]
  insufficient_samples <- sample_summary[!sufficient, ]

  # Show removed subjects if any
  if (nrow(insufficient_samples) > 0) {
    removed_subjects <- insufficient_samples |>
      dplyr::pull(!!sym(subject_ids))

    if (verbose >= 1) {
      message(
        "Removing ", nrow(insufficient_samples), " subject(s) with fewer than ",
        min_cells, " cells:\n",
        paste(paste0(removed_subjects, " (", insufficient_samples$cell_count, " cells)"),
          collapse = "\n"
        )
      )
    }
  }

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
  if (verbose >= 1) {
    message(
      "Using ", nrow(sufficient_samples), " subject(s) that have at least ",
      min_cells, " cells:\n",
      paste(paste0(valid_samples, " (", sufficient_samples$cell_count, " cells)"),
        collapse = "\n"
      )
    )
  }

  return(valid_samples)
}
