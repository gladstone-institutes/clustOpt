#' @importFrom methods as
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Large Dataset Handling Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title leverage_sketch
#' @description
#' Uses leverage score-based sampling to reduce the size of large Seurat objects
#' by creating a representative sketch assay. This method preserves the most
#' informative cells while dramatically reducing computational requirements.
#' Supports both single-cell RNA-seq and CyTOF data.
#'
#' @param input A Seurat object to be sketched
#' @param sketch_size Integer. Number of cells to include in the sketch assay.
#'   If NULL, defaults to 10\% of total cells
#' @param dtype Character. Type of data: "scRNA" (default) for single-cell
#'   RNA-seq or "CyTOF" for mass cytometry. CyTOF data should be arcsinh
#'   normalized and stored in the counts slot
#' @param skip_norm Logical. Set to TRUE if scRNA-seq data has already been
#'   normalized with `Seurat::NormalizeData()` (default FALSE). CyTOF data
#'   normalization is always skipped
#' @param on_disk Logical. Whether to use BPCells on-disk count matrices to
#'   speed up sketching for very large datasets (default FALSE)
#' @param output_dir Character. Directory path for storing on-disk count
#'   matrices when `on_disk = TRUE`. If NULL, uses temporary directory
#' @param verbose Logical. Whether to print progress messages (default TRUE)
#'
#' @return A Seurat object containing only the sketch assay, renamed to "RNA"
#'   for compatibility with downstream functions
#'
#' @details
#' Large datasets (>200,000 cells) benefit from `on_disk = TRUE` to reduce
#' memory usage during sketching.
#'
#' @examples
#' \dontrun{
#' # Basic sketching for scRNA-seq data (uses variable features)
#' sketched_obj <- leverage_sketch(seurat_obj,
#'   sketch_size = 5000,
#'   dtype = "scRNA"
#' )
#'
#' # Sketching CyTOF data (uses ALL features from marker panel)
#' cytof_sketch <- leverage_sketch(cytof_obj,
#'   sketch_size = 2000,
#'   dtype = "CyTOF"
#' )
#'
#' # Large dataset with on-disk matrices
#' large_sketch <- leverage_sketch(large_obj,
#'   sketch_size = 10000,
#'   on_disk = TRUE, verbose = TRUE
#' )
#' }
#'
#' @export
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures
#' @importFrom Seurat SketchData DefaultAssay DietSeurat
leverage_sketch <- function(input,
                            sketch_size,
                            dtype = "scRNA",
                            skip_norm = FALSE,
                            on_disk = FALSE,
                            output_dir = NULL,
                            verbose = TRUE) {
  # Validate input parameters
  if (!(dtype %in% c("scRNA", "CyTOF"))) {
    stop("dtype must be either 'scRNA' or 'CyTOF'")
  }

  if (is.null(sketch_size)) {
    if (verbose) {
      message("No sketch_size specified, defaulting to 10% of cells")
    }
    sketch_size <- ncol(input) * 0.1
  }

  if (verbose) {
    message(sprintf(
      "Sketching %s data: %d -> %d cells",
      dtype, ncol(input), as.integer(sketch_size)
    ))
  }
  if (!is.null(input@meta.data[["leverage.score"]])) {
    message("\nRemoving previously calculated leverage scores...")
    input@meta.data[["leverage.score"]] <- NULL
  }

  # Convert to on-disk format if requested
  if (on_disk) {
    if (!requireNamespace("BPCells", quietly = TRUE)) {
      stop("The BPCells package must be installed to use on_disk")
    }
    if (verbose) {
      message("Converting to on-disk format before sketching...")
    }
    # Use the convert_seurat_to_bpcells function to convert to on-disk format
    input <- convert_seurat_to_bpcells(input, output_dir = output_dir)
  }


  # Handle data type-specific normalization
  if (dtype == "scRNA" && !skip_norm) {
    if (verbose) {
      message("Normalizing scRNA-seq data...")
    }
    input <- Seurat::NormalizeData(input)
  } else if (dtype == "CyTOF") {
    if (verbose) {
      message("CyTOF data detected - skipping normalization
              (expected to be arcsinh normalized)")
    }
  } else if (dtype == "scRNA" && skip_norm) {
    if (verbose) {
      message("Skipping normalization for scRNA-seq data as requested")
    }
  }

  # Handle feature selection based on data type
  if (dtype == "scRNA") {
    if (verbose) {
      message("Finding variable features for scRNA-seq data...")
    }
    input <- Seurat::FindVariableFeatures(input)
    features_to_use <- Seurat::VariableFeatures(input)
    if (verbose) {
      message(sprintf(
        "Using %d variable features for sketching",
        length(features_to_use)
      ))
    }
  } else if (dtype == "CyTOF") {
    # For CyTOF, use all features since they represent a curated panel of markers
    features_to_use <- rownames(input)
    if (verbose) {
      message(sprintf(
        "Using all %d features for CyTOF sketching",
        length(features_to_use)
      ))
    }
  }

  input <- Seurat::SketchData(
    object = input,
    ncells = sketch_size,
    method = "LeverageScore",
    sketched.assay = "sketch",
    features = features_to_use
  )
  Seurat::DefaultAssay(input) <- "sketch"
  # Return only the sketch assay, renaming it to "RNA"
  # to avoid issues with functions that expect "RNA"
  input <- Seurat::DietSeurat(input, assays = "sketch")

  # Suppress expected warning about key conflict when renaming sketch -> RNA
  return(suppressWarnings(RenameAssays(object = input, sketch = "RNA")))
}

#' @title convert_seurat_to_bpcells
#' @description
#' Convert a Seurat object to BPCells on-disk format
#'
#' Converts the counts matrix of specified assays in a Seurat object to BPCells
#' format, saving the matrices on disk and updating the Seurat object to use
#' these on-disk matrices. Only works for single count layers.
#'
#' @param seurat_obj A Seurat object to be converted.
#' @param output_dir Directory where the BPCells matrices will be saved.
#' Defaults to a subdirectory in the system's temporary directory named after
#' the Seurat object.
#' @param assays Character vector of assays to convert. Defaults to "RNA".
#' @return The updated Seurat object using on-disk matrices.
#' @export
#' @examples
#' \dontrun{
#' # Example usage
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

    # Convert to v5 assay (suppress expected warning about Assay -> Assay5)
    suppressWarnings(
      seurat_obj[[assay_name]] <- as(
        object = seurat_obj[[assay_name]],
        Class = "Assay5"
      )
    )
    # Check if the counts matrix is already in BPCells format
    if (inherits(seurat_obj[[assay_name]]@layers$counts, "BPMatrix")) {
      message(paste(
        "Counts matrix for assay",
        assay_name, "is already in BPCells format. Skipping."
      ))
      next
    }

    # Write counts matrix to BPCells format
    counts_dir <- file.path(output_dir, paste0(assay_name, "_counts"))
    BPCells::write_matrix_dir(
      mat = seurat_obj[[assay_name]]@layers$counts, dir = counts_dir,
      overwrite = TRUE
    )

    # Update the counts matrix to on-disk BPCells matrix
    seurat_obj[[assay_name]]@layers$counts <- BPCells::open_matrix_dir(dir = counts_dir)
  }

  # Return the updated Seurat object
  return(seurat_obj)
}
