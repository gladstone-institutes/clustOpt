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
      median_score = stats::median(avg_width, na.rm = TRUE),
      variance_score = stats::var(avg_width, na.rm = TRUE),
      standard_error_score = stats::sd(avg_width,
        na.rm = TRUE
      ) / sqrt(length(avg_width)),
      cluster_median_score = stats::median(cluster_median_widths, na.rm = TRUE),
      cluster_variance_score = stats::var(cluster_median_widths, na.rm = TRUE),
      cluster_standard_error_score = stats::sd(cluster_median_widths, na.rm = TRUE) /
        sqrt(length(avg_width))
    )
}
#' Convert a Seurat object to BPCells on-disk format
#'
#' Converts the counts matrix of specified assays in a Seurat object to BPCells
#' format, saving the matrices on disk and updating the Seurat object to use
#' these on-disk matrices.
#'
#' @param seurat_obj A Seurat object to be converted.
#' @param output_dir Directory where the BPCells matrices will be saved. Defaults
#' to a subdirectory in the system's temporary directory named after the Seurat object.
#' @param assays Character vector of assays to convert. Defaults to "RNA".
#' @return The updated Seurat object using on-disk matrices.
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- readRDS("/path/to/your/seurat_object.rds")
#' seurat_obj <- convert_seurat_to_bpcells(seurat_obj)
#' }
convert_seurat_to_bpcells <- function(seurat_obj, output_dir = NULL,
                                      assays = "RNA") {
  # Derive the name of the seurat object
  obj_name <- deparse(substitute(seurat_obj))

  # Set default output directory to TMPDIR using the object's name
  if (is.null(output_dir)) {
    output_dir <- file.path(tempdir(), obj_name)
  }

  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Iterate over each specified assay in the Seurat object
  for (assay_name in assays) {
    if (!assay_name %in% names(seurat_obj@assays)) {
      warning(paste(
        "Assay", assay_name,
        "not found in the Seurat object. Skipping."
      ))
      next
    }

    assay <- seurat_obj[[assay_name]]

    # Check if the counts matrix is already in BPCells format to avoid reprocessing
    if (inherits(assay@counts, "BPMatrix")) {
      message(paste(
        "Counts matrix for assay",
        assay_name, "is already in BPCells format. Skipping."
      ))
      next
    }

    # Write counts matrix to BPCells format
    counts_dir <- file.path(output_dir, paste0(assay_name, "_counts"))
    BPCells::write_matrix_dir(mat = assay@counts, dir = counts_dir, overwrite = TRUE)

    # Update the counts matrix to on-disk BPCells matrix
    assay@counts <- BPCells::open_matrix_dir(dir = counts_dir)
    seurat_obj[[assay_name]] <- assay
  }

  # Return the updated Seurat object
  return(seurat_obj)
}
