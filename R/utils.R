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
      cluster_standard_error_score = stats::sd(cluster_median_widths, na.rm = TRUE) /
        sqrt(length(cluster_median_widths))
    )
}
