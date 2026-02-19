#' @include validation.R sketching.R data_preparation.R
#' @importFrom rlang .data
NULL

# Internal timing helper
.elapsed <- function(t0) {
  sprintf("%.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main Algorithm
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title clust_opt
#' @description
#' Runs the main resolution optimization algorithm
#'
#' @param input Seurat object
#' @param ndim Number of principal components to use.
#' @param dtype Type of data in the Seurat object "scRNA" or "CyTOF", default
#' is "scRNA". CyTOF data is expected to be arcsinh normalized (in the counts
#' slot). Sketching is supported for both data types.
#' @param sketch_size Number of cells to use for sketching.
#' @param skip_sketch Skip sketching, by default any input with more than
#' 200,000 cells is sketched to 10\% of the cells.
#' @param subject_ids Metadata field that identifies unique subjects.
#' @param res_range Range of resolutions to test.
#' @param verbose Integer verbosity level: 0 = silent, 1 = key milestones,
#' 2 = detailed progress, 3 = includes Seurat function output.
#' @param within_batch Batch variable, for a given sample only those with the
#' same value for the batch variable will be used for training.
#' @param num_trees Number of trees to use in the random forest.
#' @param train_with Either "odd" or "even" PCs for clustering and training.
#' Default is "even". It is recommended to keep train_with set to "even" so
#' that the 1st PC is in the set used to calculate silhouette scores.
#' @param min_cells Minimum cells per subject, default is 50
#' @param rf_num_threads Number of threads for ranger random forest. Default is
#' 1 to avoid contention with \code{\link[future.apply]{future_lapply}} workers.
#' When using \code{future::plan(sequential)} (the default), increase to utilize
#' available cores.
#' @return A data.frame containing a distribution of silhouette scores for each
#' resolution.
#'
#' @details
#' The clustOpt algorithm works by:
#' \enumerate{
#'   \item Sketching large datasets using leverage score-based sampling (if needed)
#'   \item Splitting principal components into independent odd/even spaces
#'   \item Performing subject-wise cross-validation
#'   \item Training random forests on cluster assignments
#'   \item Evaluating clustering quality using silhouette scores
#' }
#'
#' Both scRNA-seq and CyTOF data types support sketching for improved performance
#' on large datasets. For CyTOF data, normalization is skipped as data should
#' already be arcsinh transformed.
#'
#' @examples
#' \dontrun{
#' # Basic usage with scRNA-seq data
#' results <- clust_opt(seurat_obj, ndim = 50, subject_ids = "donor_id")
#'
#' # CyTOF data analysis
#' cytof_results <- clust_opt(cytof_obj,
#'   ndim = 30, dtype = "CyTOF",
#'   subject_ids = "sample_id"
#' )
#'
#' # Large dataset with custom sketch size
#' large_results <- clust_opt(large_obj,
#'   ndim = 50, sketch_size = 10000,
#'   subject_ids = "donor_id"
#' )
#' }
#'
#' @export
#' @importFrom Seurat DefaultAssay DietSeurat RunPCA FindNeighbors FindClusters
#' @importFrom Seurat ScaleData FindVariableFeatures
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor handlers
#' @importFrom purrr map_df
clust_opt <- function(input,
                      ndim,
                      dtype = "scRNA",
                      sketch_size = NULL,
                      skip_sketch = FALSE,
                      subject_ids,
                      res_range = c(
                        0.02, 0.04, 0.06, 0.08, 0.1,
                        0.2, 0.4, 0.6, 0.8, 1, 1.2
                      ),
                      within_batch = NA,
                      verbose = 0,
                      num_trees = 1000,
                      train_with = "even",
                      min_cells = 50,
                      rf_num_threads = 1) {
  verbose <- normalize_verbose(verbose)
  t0_total <- Sys.time()

  # Make sure a seed is set, only setting if the user has not set one
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    t <- as.integer(Sys.time())
    if (verbose >= 1) message("Setting seed: ", t)
    set.seed(t)
  }

  if (!(train_with %in% c("odd", "even"))) {
    stop("train_with can only 'odd' or 'even'")
  }
  if (!(dtype %in% c("CyTOF", "scRNA"))) {
    stop("dtype is not one of 'CyTOF' or 'scRNA'")
  }


  if (!skip_sketch) {
    if (check_size(input) || !is.null(sketch_size)) {
      if (verbose >= 1) message("Sketching input data")
      t0_sketch <- Sys.time()
      input <- leverage_sketch(input, sketch_size, dtype, verbose = verbose)
      if (verbose >= 1) message(sprintf("[sketch] %s", .elapsed(t0_sketch)))
    } else {
      if (verbose >= 1) message("Input is small enough to run with all cells")
    }
  }


  sample_names <- get_valid_samples(input, subject_ids, min_cells,
                                    verbose = verbose)
  if (is.null(sample_names)) {
    stop(paste0(
      "Unable to perform cluster resolution optimization for this data. ",
      "There are less than 3 subjects with at least", min_cells, "cells."
    ))
  }

  # Get every combination of test sample and resolution
  runs <- expand.grid(sample_names, res_range)
  if (verbose >= 1) {
    message(paste0(
      "Found ", nrow(runs),
      " combinations of test subject and resolution"
    ))
  }

  # Set up progress logging
  progressr::handlers("progress")
  p <- progressr::progressor(along = unique(runs[, 1]))

  res <- NULL
  for (sam in unique(runs[, 1])) {
    t0_sam <- Sys.time()
    if (verbose >= 1) message(paste0("Holdout subject: ", sam))
    if (verbose >= 2) {
      message(paste0("Preparing training data.."))
    }

    t0_prep_train <- Sys.time()
    train <- prep_train(
      input = input,
      dtype = dtype,
      subject_ids = subject_ids,
      within_batch = within_batch,
      test_id = sam,
      verbose = verbose
    )
    if (verbose >= 1) message(sprintf("[prep_train] %s", .elapsed(t0_prep_train)))

    if (verbose >= 2) {
      message(paste0("Preparing test data.."))
    }

    t0_prep_test <- Sys.time()
    test <- prep_test(
      input = input,
      dtype = dtype,
      subject_ids = subject_ids,
      test_id = sam,
      verbose = verbose
    )
    if (verbose >= 1) message(sprintf("[prep_test] %s", .elapsed(t0_prep_test)))

    t0_pca <- Sys.time()
    if (dtype == "scRNA") {
      train <- Seurat::RunPCA(train,
        npcs = ndim,
        verbose = (verbose >= 3),
        assay = "SCT"
      )

      clust_pcs <- switch(train_with,
        odd = "even_pca",
        even = "odd_pca"
      )
      train_with_pcs <- switch(train_with,
        odd = "odd_pca",
        even = "even_pca"
      )

      # Create 2 separate PCA reductions
      train <- split_pca_dimensions(train, verbose)
      if (verbose >= 1) {
        message(sprintf("[RunPCA + split_pca_dimensions] %s", .elapsed(t0_pca)))
      }

      if (verbose >= 2) {
        message(sprintf("Clustering with %s", clust_pcs))
      }

      t0_clust <- Sys.time()
      train <- Seurat::FindNeighbors(
        object = train,
        dims = seq_len(ncol(train@reductions[[clust_pcs]]@cell.embeddings)),
        verbose = (verbose >= 3),
        reduction = clust_pcs
      )

      train <- Seurat::FindClusters(
        object = train,
        resolution = res_range,
        verbose = (verbose >= 3)
      )
      if (verbose >= 1) {
        message(sprintf("[FindNeighbors + FindClusters] %s", .elapsed(t0_clust)))
      }
    } else {
      train <- Seurat::ScaleData(train, features = NULL, verbose = (verbose >= 3))
      train <- Seurat::FindVariableFeatures(train,
        selection.method = "vst", nfeatures = ndim,
        verbose = (verbose >= 3)
      )
      train <- Seurat::RunPCA(train,
        npcs = ndim,
        approx = FALSE,
        verbose = (verbose >= 3)
      )

      test <- Seurat::ScaleData(test, features = NULL, verbose = (verbose >= 3))

      clust_pcs <- switch(train_with,
        odd = "even_pca",
        even = "odd_pca"
      )
      train_with_pcs <- switch(train_with,
        odd = "odd_pca",
        even = "even_pca"
      )
      # Create 2 separate PCA reductions
      train <- split_pca_dimensions(train, verbose)
      if (verbose >= 1) {
        message(sprintf("[RunPCA + split_pca_dimensions] %s", .elapsed(t0_pca)))
      }

      t0_clust <- Sys.time()
      train <- Seurat::FindNeighbors(
        object = train,
        dims = seq_len(ncol(train@reductions[[clust_pcs]]@cell.embeddings)),
        verbose = (verbose >= 3),
        reduction = clust_pcs
      )

      train <- Seurat::FindClusters(
        object = train,
        resolution = res_range,
        verbose = (verbose >= 3)
      )
      if (verbose >= 1) {
        message(sprintf("[FindNeighbors + FindClusters] %s", .elapsed(t0_clust)))
      }
    }
    if (verbose >= 2) {
      message("Clustering complete..")
    }

    t0_project <- Sys.time()

    # Pre-extract cluster assignments using exact column names
    prefix <- if (dtype == "scRNA") "SCT_snn_res." else "RNA_snn_res."
    cluster_by_res <- lapply(res_range, function(r) {
      train@meta.data[[paste0(prefix, r)]]
    })
    names(cluster_by_res) <- as.character(res_range)

    df_list <- project_pca(
      train_seurat = train,
      test_seurat = test,
      train_with_pcs = train_with_pcs,
      clust_pcs = clust_pcs,
      dtype = dtype,
      verbose = verbose
    )
    rm(train, test)
    if (verbose >= 1) message(sprintf("[project_pca] %s", .elapsed(t0_project)))

    if (verbose >= 1) {
      message(sprintf(
        "[data dimensions] train: %d cells, test: %d cells, %d PCs, %d resolutions",
        nrow(df_list[["train_proj_train_with_pcs"]]),
        nrow(df_list[["test_proj_train_with_pcs"]]),
        ncol(df_list[["test_proj_train_with_pcs"]]),
        length(res_range)
      ))
    }

    # Precompute SNN graph (resolution-invariant)
    t0_snn <- Sys.time()
    coords_mat <- df_list[["test_proj_clust_pcs"]]
    rownames(coords_mat) <- seq_len(nrow(coords_mat))
    precomputed_snn <- Seurat::FindNeighbors(
      coords_mat, k.param = 20, compute.SNN = TRUE,
      prune.SNN = 1 / 15, verbose = (verbose >= 3)
    )[["snn"]]
    if (verbose >= 1) message(sprintf("[precompute SNN] %s", .elapsed(t0_snn)))

    # Precompute distance matrix for silhouette (resolution-invariant)
    t0_dist <- Sys.time()
    precomputed_dist <- dist(df_list[["test_proj_clust_pcs"]])
    if (verbose >= 1) {
      message(sprintf(
        "[precompute dist] %s (%d cells)",
        .elapsed(t0_dist), nrow(df_list[["test_proj_clust_pcs"]])
      ))
    }

    t0_rf <- Sys.time()
    this_result <- future.apply::future_lapply(
      rownames(runs[runs$Var1 == sam, ]),
      function(i) {
        train_random_forest(
          res = runs[i, 2],
          df_list = df_list,
          cluster_by_res = cluster_by_res,
          sam = runs[i, 1],
          num_trees = num_trees,
          verbose = verbose,
          snn_graph = precomputed_snn,
          precomputed_dist = precomputed_dist,
          rf_num_threads = rf_num_threads
        )
      },
      future.seed = TRUE,
      future.packages = c("SeuratObject", "Seurat", "ranger", "cluster")
    )
    if (verbose >= 1) {
      message(sprintf(
        "[future_lapply RF] %s (%d resolutions)",
        .elapsed(t0_rf), length(res_range)
      ))
    }
    res <- c(res, this_result)
    p()
    if (verbose >= 1) message(sprintf("[sample %s total] %s", sam, .elapsed(t0_sam)))
  }

  if (verbose >= 1) message(sprintf("[total pipeline] %s", .elapsed(t0_total)))
  purrr::map_df(res, .f = as.data.frame)
}
