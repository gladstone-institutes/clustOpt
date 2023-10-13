#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' check_size: Checks if input Seurat object is small enough to run
#' sil_score
#'
#' @param input object to check
#' @param verbose print warning about size
#' @return NULL
#'
#' @export
#'
check_size <- function(input, verbose = TRUE) {
  if (!methods::is(input, "Seurat")) {
    stop("Input must be a seurat object")
  }
  if (ncol(input) >= 2E5 | Seurat::Assays(input) != "sketch") {
    if (verbose) {
      warning(paste0(
        "Over 200,000 cells detected, reccommend running with ",
        "downsample = 'sketch'"
      ))
    }
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

  options(Seurat.object.assay.version = "v5")
  options(future.globals.maxSize = 1e9)
  
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
  Seurat::DietSeurat(input, assays = "sketch")
}


