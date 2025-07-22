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
  if (!inherits(input, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (ncol(input) >= 2E5) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Leverage Score-based Sketching for Large Datasets
#' 
#' @description
#' Uses leverage score-based sampling to reduce the size of large Seurat objects
#' by creating a representative sketch assay. This method preserves the most 
#' informative cells while dramatically reducing computational requirements.
#' Supports both single-cell RNA-seq and CyTOF data.
#'
#' @param input A Seurat object to be sketched
#' @param sketch_size Integer. Number of cells to include in the sketch assay.
#'   If NULL, defaults to 10% of total cells
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
#' The function automatically handles data type-specific preprocessing:
#' \itemize{
#'   \item For scRNA-seq: Normalizes data (unless `skip_norm = TRUE`) and uses variable features for sketching
#'   \item For CyTOF: Skips normalization (data should be arcsinh normalized) and uses ALL features for sketching since they represent a curated marker panel
#' }
#' 
#' Large datasets (>200,000 cells) benefit from `on_disk = TRUE` to reduce 
#' memory usage during sketching.
#'
#' @examples
#' \dontrun{
#' # Basic sketching for scRNA-seq data (uses variable features)
#' sketched_obj <- leverage_sketch(seurat_obj, sketch_size = 5000, dtype = "scRNA")
#' 
#' # Sketching CyTOF data (uses ALL features from marker panel)
#' cytof_sketch <- leverage_sketch(cytof_obj, sketch_size = 2000, dtype = "CyTOF")
#' 
#' # Large dataset with on-disk matrices
#' large_sketch <- leverage_sketch(large_obj, sketch_size = 10000, 
#'                                on_disk = TRUE, verbose = TRUE)
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
    message(sprintf("Sketching %s data: %d -> %d cells", 
                   dtype, ncol(input), as.integer(sketch_size)))
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
      message(sprintf("Using %d variable features for sketching",
                      length(features_to_use)))
    }
  } else if (dtype == "CyTOF") {
    # For CyTOF, use all features since they represent a curated panel of markers
    features_to_use <- rownames(input)
    if (verbose) {
      message(sprintf("Using all %d features for CyTOF sketching", 
                     length(features_to_use)))
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

  return(RenameAssays(object = input, sketch = "RNA"))
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

    # Convert to v5 assay
    seurat_obj[[assay_name]] <- as(object = seurat_obj[[assay_name]],
                                   Class = "Assay5")
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
#' @importFrom dplyr group_by summarize filter n sym
get_valid_samples <- function(input, subject_ids, min_cells) {
  # Summarize the number of cells per sample
  sample_summary <- input@meta.data |>
    dplyr::group_by(!!dplyr::sym(subject_ids)) |>
    dplyr::summarize(cell_count = dplyr::n(), .groups = "drop") |>
    dplyr::ungroup()

  # Split samples into sufficient and insufficient
  sufficient_samples <- sample_summary |>
    dplyr::filter(cell_count >= min_cells)

  insufficient_samples <- sample_summary |>
    dplyr::filter(cell_count < min_cells)

  # Show removed subjects if any
  if (nrow(insufficient_samples) > 0) {
    removed_subjects <- insufficient_samples |>
      dplyr::pull(!!sym(subject_ids))

    message(
      "Removing ", nrow(insufficient_samples), " subject(s) with fewer than ",
      min_cells, " cells:\n",
      paste(paste0(removed_subjects, " (", insufficient_samples$cell_count, " cells)"),
        collapse = "\n"
      )
    )
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
  message(
    "Using ", nrow(sufficient_samples), " subject(s) that have at least ",
    min_cells, " cells:\n",
    paste(paste0(valid_samples, " (", sufficient_samples$cell_count, " cells)"),
      collapse = "\n"
    )
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

#' Calculate Adjusted Rand Index
#'
#' Computes the Adjusted Rand Index (ARI) to measure the similarity between
#' two clustering assignments. The ARI is a measure of agreement between two
#' partitions, adjusted for chance. Values range from 0 (random partitioning)
#' to 1 (perfect agreement), with negative values indicating worse than random.
#'
#' @param seurat_obj A Seurat object containing the metadata with clustering assignments
#' @param meta1 Character string specifying the first metadata column name containing
#'   cluster assignments
#' @param meta2 Character string specifying the second metadata column name containing
#'   cluster assignments to compare against meta1
#'
#' @return Numeric value representing the Adjusted Rand Index between the two
#'   clustering assignments
#'
#' @details
#' The Adjusted Rand Index is calculated using the formula:
#' ARI = (RI - Expected_RI) / (max(RI) - Expected_RI)
#'
#' Where:
#' - RI is the Rand Index
#' - Expected_RI is the expected value of RI under random partitioning
#' - max(RI) is the maximum possible value of RI
#'
#' @examples
#' \dontrun{
#' # Compare two clustering results in a Seurat object
#' ari_score <- adjusted_rand_index(seurat_obj, "seurat_clusters", "leiden_clusters")
#' }
#'
#' @export
adjusted_rand_index <- function(seurat_obj, meta1, meta2) {
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  if (!is.character(meta1) || !is.character(meta2)) {
    stop("meta1 and meta2 must be character strings")
  }

  if (length(meta1) != 1 || length(meta2) != 1) {
    stop("meta1 and meta2 must be single character strings")
  }

  # Check if metadata columns exist
  meta <- seurat_obj@meta.data
  if (!meta1 %in% colnames(meta)) {
    stop(paste("Column", meta1, "not found in metadata"))
  }

  if (!meta2 %in% colnames(meta)) {
    stop(paste("Column", meta2, "not found in metadata"))
  }

  # Extract groupings
  group1 <- meta[[meta1]]
  group2 <- meta[[meta2]]

  # Check for missing values
  if (any(is.na(group1)) || any(is.na(group2))) {
    warning("Missing values detected in clustering assignments")
    # Remove cells with missing values in either grouping
    valid_cells <- !is.na(group1) & !is.na(group2)
    group1 <- group1[valid_cells]
    group2 <- group2[valid_cells]
  }

  # Check if we have enough data
  if (length(group1) < 2) {
    stop("Need at least 2 observations to calculate ARI")
  }

  if (length(group1) != length(group2)) {
    stop("group1 and group2 must have the same length")
  }

  # Contents of columns are put into a matrix
  tab <- table(group1, group2)
  # n represents total number of observations being compared
  n <- length(group1)

  sum_comb_ij <- sum(choose(tab, 2)) # Σij (nij 2)
  sum_comb_rows <- sum(choose(rowSums(tab), 2)) # Σi (ai 2)
  sum_comb_columns <- sum(choose(colSums(tab), 2)) # Σj (bj 2)
  total_pairs <- choose(n, 2) # (n/2)

  # Several of the sums are operated on and assigned to new variable
  expected_index <- (sum_comb_rows * sum_comb_columns) / total_pairs
  max_index <- 0.5 * (sum_comb_rows + sum_comb_columns)

  # Perfect agreement when both clusterings are identical singletons
  if (max_index == expected_index) {
    return(0)
  }

  # This is the ARI equation and will output our value
  ari_result <- (sum_comb_ij - expected_index) / (max_index - expected_index)
  return(ari_result)
}
