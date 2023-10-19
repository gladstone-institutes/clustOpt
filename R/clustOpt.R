#' @include utils.R
#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Runs the main resolution optimization algorithm
#'
#' @param input Seurat object
#' @param subject_ids Metadata field that identifies unique samples.
#' @param ndim Number of principal components to use.
#' @param dtype Type of data in the Seurat object "scRNA" or "CyTOF", default
#' is "scRNA". CyTOF data is expected to be arcsinh normalized.
#' @param res_range Range of resolutions to test.
#' @param ncores Number of cores.
#' @param verbose print messages.
#' @param within_batch Batch variable, for a given sample only those with the
#' same value for the batch variable will be used for training.
#' @param plan Type of future::plan for parallelization, "multicore" or
#' "multisession". Use "multisession" if running in RStudio or on a Windows OS.
#' @return A data.frame containing a distribution of silhouette scores for each
#' resolution.
#'
#' @export
#'
clust_opt <- function(input,
                      ndim,
                      plan = "multicore",
                      dtype = "scRNA",
                      tmp_dir = "~/tmp",
                      sketch_size = NULL,
                      subject_ids,
                      seed = 1,
                      res_range = c(
                        0.02, 0.04, 0.06, 0.08, 0.1,
                        0.2, 0.4, 0.6, 0.8, 1, 1.2
                      ),
                      within_batch = NA,
                      ncore = 4,
                      verbose = TRUE) {
  sample_names <- as.vector((unique(input@meta.data[[subject_ids]])))

  if (!(dtype %in% c("CyTOF", "scRNA"))) {
    stop("dtype is not one of 'CyTOF' or 'scRNA'")
  }

  # Clear any previous normalizations
  Seurat::DefaultAssay(input) <- "RNA"
  input <- Seurat::DietSeurat(input, assays = "RNA")

  if (check_size(input) || !is.null(sketch_size)) {
    message("Sketching input data")
    input <- leverage_sketch(input, sketch_size)
  } else {
    message("Input is small enough to run with all cells")
  }
  set.seed(seed)
  
  
  # Get every combination of test sample and resolution
  runs <- expand.grid(sample_names, res_range)
  message(paste0(
    "Found ", nrow(runs),
    " combinations of test sample and resolution"
  ))
  res_list <- vector("list", nrow(runs))
  iter <- 1
  for (sam in unique(runs[, 1])) {
    future::plan("sequential")
    message(paste0("Holdout sample: ", sam))
    if (verbose) {
      message(paste0("Preparing training data.."))
    }
    train <- prep_train(
      input = input,
      dtype = dtype,
      subject_ids = subject_ids,
      within_batch = within_batch,
      test_id = sam
    )
    if (verbose) {
      message(paste0("Preparing test data.."))
    }
    test <- prep_test(
      input = input,
      dtype = dtype,
      subject_ids = subject_ids,
      test_id = sam
    )

    message(paste0(
      "Found ",
      length(intersect(
        rownames(test@assays[["SCT"]]@scale.data),
        Seurat::VariableFeatures(train)
      )),
      " shared genes between testing and training data"
    ))
    
    train <- Seurat::RunPCA(train,
      npcs = ndim,
      verbose = FALSE,
      assay = "SCT"
    )
    train <- Seurat::FindNeighbors(
      object = train,
      dims = 1:ndim,
      verbose = verbose
    )
    # Set up parallelization
    future::plan(plan, workers = ncore)
    
    train <- Seurat::FindClusters(
      object = train,
      resolution = res_range,
      verbose = verbose
    )

    # Iterate through the combinations in parallel
    progressr::handlers("progress")
    progressr::with_progress({
      p <- progressr::progressor(along = res_range)
      result <- future.apply::future_lapply(
        future.packages = c("Seurat"),
        FUN = function(x) {
          p()
          return(train)
        },X = list(res_range), future.seed = TRUE
      )
    })
    res_list[[iter]] <- result
    iter <- iter + 1
  }
}

#' Prepare training data for random forest
#'
#' @param input Seurat object
#' @param subject_ids Metadata field that identifies unique samples.
#' @param dtype Type of data in the Seurat object "scRNA" or "CyTOF", default
#' is "scRNA". CyTOF data is expected to be arcsinh normalized.
#' @param within_batch Batch variable, for a given sample only those with the
#' same value for the batch variable will be used for training.
#' @param test_id subject_id for the test sample
#' @return Training data formatted for sil_score, format depends on dtype
#'
#' @export
#'
prep_train <- function(input,
                       subject_ids,
                       dtype = "scRNA",
                       within_batch,
                       test_id) {
  if (dtype == "scRNA") {
    # If within_batch is provided, then use only training samples from the
    # same batch
    if (!is.na(within_batch)) {
      # Get the batch of test_id
      this_batch <- input@meta.data |>
        dplyr::filter(get(subject_ids) == test_id) |>
        dplyr::pull(get(within_batch)) |>
        unique()

      if (length(this_batch) > 1) {
        stop("More than one batch found for this sample")
      }
      
      train_cells <- input@meta.data |>
        dplyr::filter(get(subject_ids) != test_id &
          get(within_batch) == this_batch) |>
        rownames()
      train_seurat <- subset(input, cells = train_cells)
      
      # Normalize the training samples
      train_seurat <- Seurat::SCTransform(train_seurat,
        assay = DefaultAssay(train_seurat),
        vst.flavor = "v2",
        verbose = FALSE
      )
      train_seurat <- Seurat::DietSeurat(train_seurat, assays = "SCT")
      return(train_seurat)
    } else {
      # Return all other samples
      Idents(input) <- subject_ids
      train_cells <- WhichCells(object = input,idents = test_id, invert = TRUE)
      train_seurat <- subset(input, cells = train_cells)
      # Normalize the training samples
      train_seurat <- Seurat::SCTransform(train_seurat,
        assay = DefaultAssay(train_seurat),
        vst.flavor = "v2",
        verbose = FALSE
      )
      train_seurat <- Seurat::DietSeurat(train_seurat, assays = "SCT")
      return(train_seurat)
    }
  } else {
    return(NULL)
  }
}
#' Prepare test data for random forest
#'
#' @param input Seurat object
#' @param subject_ids Metadata field that identifies unique samples.
#' @param dtype Type of data in the Seurat object "scRNA" or "CyTOF", default
#' is "scRNA". CyTOF data is expected to be arcsinh normalized.
#' @param test_id subject_id for the test sample
#' @return Training data formatted for sil_score, format depends on dtype
#'
#' @export
#'
prep_test <- function(input,
                      subject_ids,
                      dtype = "scRNA",
                      test_id) {
  if (dtype == "scRNA") {
    
    Idents(input) <- subject_ids
    test_cells <- WhichCells(object = input,idents = test_id)
    test_seurat <- subset(input, cells = test_cells)
    
    test_seurat <- Seurat::SCTransform(test_seurat,
      vst.flavor = "v2",
      assay = DefaultAssay(test_seurat),
      verbose = FALSE,
      variable.features.n = length(rownames(test_seurat)),
      return.only.var.genes = FALSE,
      min_cells = 1
    )
    test_seurat <- Seurat::DietSeurat(test_seurat, assays = "SCT")
    
    return(test_seurat)
  } else {
    return(NULL)
  }
}
