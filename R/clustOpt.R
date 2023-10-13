#' @include utils.R
#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' clust_opt: Runs the main resolution optimization algorithm
#'
#' @param input Seurat object
#' @param subject_ids Metadata field that identifies unique samples.
#' @param ndim Number of principal components to use.
#' @param norm_method Method to use for normalizing counts, "SCT" or "log",
#' default is "SCT". Ignored if dtype = "CyTOF".
#' @param dtype Type of data in the Seurat object "scRNA" or "CyTOF", default
#' is "scRNA". CyTOF data is expected to be arcsinh normalized.
#' @param res_range Range of resolutions to test.
#' @param ncore Number of cores
#' @param verbose print messages.
#' @param within_batch Batch variable, for a given sample only those with the 
#' same value for the batch variable will be used for training.
#' @return A data.frame containing a distribution of silhouette scores for each
#' resolution.
#'
#' @export
#'
clust_opt <- function(input,
                      ndim,
                      norm_method = "SCT",
                      sketch_size = NULL,
                      subject_ids,
                      seed = 1,
                      res_range = c(0.02, 0.04, 0.06, 0.08, 0.1,
                                    0.2, 0.4, 0.6, 0.8, 1, 1.2),
                      ncore = 4,
                      verbose = TRUE) {
  n_samples <- length(unique(input@meta.data[[subject_ids]]))
  sample_names <- as.vector((unique(input@meta.data[[subject_ids]])))
  
  if (check_size(input)){
    message("Sketching input data")
    input <- leverage_sketch(input, sketch_size)
  } else {
    message("Input is small enough to run with all cells")
  }
  # Set up parallelization
  doParallel::registerDoParallel(cores = ncores)
  foreach::getDoParWorkers()
  foreach::getDoParName()
  set.seed(seed)
  # Get every combination of test sample and resolution
  runs <- expand.grid(sample_names,res_range) |>
    dplyr::rename(subject_ids = Var1,resolution = Var2)

  if (norm_method == "SCT") {

  } else {

  }
  
}
#' sil_score: Trains the RF model, predicts cluster labels for the test sample
#' and returns the average silhouette score and the cluster average silhouette 
#' score for the predicted clusters.
#'
#' @param input Seurat object, SCTnormalized
#' @param test_id
#' @param ndim Number of principal components to use.
#' @param res_range Range of resolutions to test.
#' @param verbose print messages.
#' @param within_batch Batch variable, for a given sample only those with the 
#' same value for the batch variable will be used for training.
#' @return 
#'
#' @export
#'
sil_score <- function() {
  
}

