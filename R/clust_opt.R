#' @include validation.R sketching.R data_preparation.R checkpoint.R
#' @importFrom rlang .data
NULL

# Internal timing helpers
# Format a duration (in seconds) in a human-readable way: sub-second as
# milliseconds, seconds with one decimal, and a minute or more split into
# whole-second components (e.g. "142ms", "6.3s", "2m 0s", "1h 3m 7s").
.format_duration <- function(secs) {
  if (!is.finite(secs)) {
    return(sprintf("%.1fs", secs))
  }
  if (secs < 1) {
    return(sprintf("%.0fms", secs * 1000))
  }
  if (secs < 60) {
    return(sprintf("%.1fs", secs))
  }
  total <- round(secs)
  if (total < 3600) {
    return(sprintf("%dm %ds", total %/% 60, total %% 60))
  }
  sprintf("%dh %dm %ds", total %/% 3600, (total %% 3600) %/% 60, total %% 60)
}

.elapsed <- function(t0) {
  .format_duration(as.numeric(difftime(Sys.time(), t0, units = "secs")))
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
#' 2 = detailed progress, 3 = includes Seurat function output, 4 = includes
#' output from other packages such as ranger.
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
#' @param checkpoint_dir Optional directory enabling checkpoint/resume (default
#' \code{NULL}, off). Each holdout subject's result is saved as it completes;
#' rerunning with the same directory skips finished subjects, so a timed-out HPC
#' job resumes instead of restarting. Reusing a directory with a different
#' configuration or input throws an error.
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
                      rf_num_threads = 1,
                      checkpoint_dir = NULL) {
  verbose <- normalize_verbose(verbose)
  t0_total <- Sys.time()
  if (verbose >= 1) cli::cli_rule(left = "{.pkg clustOpt}")

  # Make sure a seed is set, only setting if the user has not set one
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    t <- as.integer(Sys.time())
    if (verbose >= 1) cli::cli_alert_info("Setting seed: {t}")
    set.seed(t)
  }

  if (!(train_with %in% c("odd", "even"))) {
    stop("train_with can only 'odd' or 'even'")
  }
  if (!(dtype %in% c("CyTOF", "scRNA"))) {
    stop("dtype is not one of 'CyTOF' or 'scRNA'")
  }

  # Checkpoint setup. Resolve a concrete integer run seed (deterministic from the
  # caller's RNG state, so an upstream set.seed() still controls it) and persist
  # it plus a config fingerprint. On resume the manifest supplies the original
  # seed and rejects a mismatched config.
  checkpointing <- !is.null(checkpoint_dir)
  run_seed <- NULL
  sketch_path <- NULL
  if (checkpointing) {
    if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
    run_seed <- sample.int(.Machine$integer.max, 1L)
    config <- .checkpoint_config(
      input, ndim, dtype, sketch_size, skip_sketch, subject_ids, res_range,
      within_batch, num_trees, train_with, min_cells
    )
    run_seed <- .checkpoint_init_manifest(checkpoint_dir, config, run_seed)
    set.seed(run_seed)
    sketch_path <- file.path(checkpoint_dir, "sketch.rds")
  }

  if (checkpointing && file.exists(sketch_path)) {
    if (verbose >= 1) {
      cli::cli_alert_info("Resuming: loading sketched input from checkpoint")
    }
    input <- readRDS(sketch_path)
  } else if (!skip_sketch) {
    if (check_size(input) || !is.null(sketch_size)) {
      if (verbose >= 1) cli::cli_alert_info("Sketching input data")
      t0_sketch <- Sys.time()
      input <- leverage_sketch(input, sketch_size, dtype, verbose = verbose)
      if (verbose >= 1) cli::cli_alert_success("[sketch] {(.elapsed(t0_sketch))}")
      # Persist the sketched object so a timeout does not force a re-sketch
      if (checkpointing) .save_rds_atomic(input, sketch_path)
    } else {
      if (verbose >= 1) cli::cli_alert_info("Input is small enough to run with all cells")
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
    cli::cli_alert_info(
      "Found {nrow(runs)} combinations of test subject and resolution"
    )
  }

  # Set up progress logging
  progressr::handlers("progress")
  p <- progressr::progressor(along = unique(runs[, 1]))

  unique_samples <- unique(runs[, 1])
  res <- vector("list", length(unique_samples))
  for (sam_idx in seq_along(unique_samples)) {
    sam <- unique_samples[sam_idx]
    t0_sam <- Sys.time()

    # Resume: reload a completed subject and skip its work entirely. Otherwise
    # seed this subject deterministically from its name so its computation (PCA,
    # clustering, and the future.seed RF fan-out) reproduces regardless of which
    # subjects ran before it in this process.
    ckpt_path <- NULL
    if (checkpointing) {
      ckpt_path <- .subject_checkpoint_path(
        checkpoint_dir, sam_idx, sam, length(unique_samples)
      )
      if (file.exists(ckpt_path)) {
        if (verbose >= 1) {
          cli::cli_alert_info("Resuming: loaded subject {sam} from checkpoint")
        }
        res[[sam_idx]] <- readRDS(ckpt_path)
        p()
        next
      }
      set.seed(.subject_seed(run_seed, sam))
    }

    if (verbose >= 1) cli::cli_h2("Holdout subject: {sam}")
    if (verbose >= 2) {
      cli::cli_alert_info("Preparing training data...")
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
    if (verbose >= 1) cli::cli_alert_success("[prep_train] {(.elapsed(t0_prep_train))}")

    if (verbose >= 2) {
      cli::cli_alert_info("Preparing test data...")
    }

    t0_prep_test <- Sys.time()
    test <- prep_test(
      input = input,
      dtype = dtype,
      subject_ids = subject_ids,
      test_id = sam,
      verbose = verbose,
      residual_features = if (dtype == "scRNA") Seurat::VariableFeatures(train) else NULL
    )
    if (verbose >= 1) cli::cli_alert_success("[prep_test] {(.elapsed(t0_prep_test))}")

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
        cli::cli_alert_success("[RunPCA + split_pca_dimensions] {(.elapsed(t0_pca))}")
      }

      if (verbose >= 2) {
        cli::cli_alert_info("Clustering with {.field {clust_pcs}}")
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
        cli::cli_alert_success("[FindNeighbors + FindClusters] {(.elapsed(t0_clust))}")
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
        cli::cli_alert_success("[RunPCA + split_pca_dimensions] {(.elapsed(t0_pca))}")
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
        cli::cli_alert_success("[FindNeighbors + FindClusters] {(.elapsed(t0_clust))}")
      }
    }
    if (verbose >= 2) {
      cli::cli_alert_success("Clustering complete")
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
    if (verbose >= 1) cli::cli_alert_success("[project_pca] {(.elapsed(t0_project))}")

    if (verbose >= 1) {
      n_train <- nrow(df_list[["train_proj_train_with_pcs"]])
      n_test <- nrow(df_list[["test_proj_train_with_pcs"]])
      n_pcs <- ncol(df_list[["test_proj_train_with_pcs"]])
      n_res <- length(res_range)
      cli::cli_alert_info(
        "train: {n_train} cells, test: {n_test} cells, {n_pcs} PCs, {n_res} resolutions"
      )
    }

    # Precompute SNN graph (resolution-invariant)
    t0_snn <- Sys.time()
    coords_mat <- df_list[["test_proj_clust_pcs"]]
    rownames(coords_mat) <- seq_len(nrow(coords_mat))
    precomputed_snn <- Seurat::FindNeighbors(
      coords_mat, k.param = 20, compute.SNN = TRUE,
      prune.SNN = 1 / 15, verbose = (verbose >= 3)
    )[["snn"]]
    if (verbose >= 1) cli::cli_alert_success("[precompute SNN] {(.elapsed(t0_snn))}")

    # Precompute distance matrix for silhouette (resolution-invariant).
    # Under a socket-cluster plan (multisession / remote) globals are serialized
    # to every worker, so shipping the dense O(n^2) dist is wasteful. Skip it and
    # let workers recompute dist() in parallel from the already-shipped small
    # coords via the precomputed_dist = NULL fallback in calculate_silhouette_score().
    t0_dist <- Sys.time()
    serializing_plan <- inherits(future::plan(), "cluster")
    precomputed_dist <- if (serializing_plan) {
      NULL
    } else {
      dist(df_list[["test_proj_clust_pcs"]])
    }
    if (verbose >= 1) {
      n_cells <- nrow(df_list[["test_proj_clust_pcs"]])
      if (serializing_plan) {
        cli::cli_alert_info(
          "[precompute dist] skipped (parallel plan; workers compute dist) ({n_cells} cells)"
        )
      } else {
        cli::cli_alert_success("[precompute dist] {(.elapsed(t0_dist))} ({n_cells} cells)")
      }
    }

    t0_rf <- Sys.time()
    # Iterate the resolution vector directly and capture sam as a scalar so the
    # full runs data.frame is not pulled into (and serialized with) the closure.
    this_result <- future.apply::future_lapply(
      res_range,
      function(r) {
        train_random_forest(
          res = r,
          df_list = df_list,
          cluster_by_res = cluster_by_res,
          sam = sam,
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
      n_res <- length(res_range)
      cli::cli_alert_success("[future_lapply RF] {(.elapsed(t0_rf))} ({n_res} resolutions)")
    }
    res[[sam_idx]] <- this_result
    # Checkpoint this subject atomically (temp file + rename) in the main process
    if (checkpointing) .save_rds_atomic(this_result, ckpt_path)
    p()
    if (verbose >= 1) cli::cli_rule(right = "{sam} {(.elapsed(t0_sam))}")
  }

  if (verbose >= 1) {
    cli::cli_alert_success("[total pipeline] {(.elapsed(t0_total))}")
    cli::cli_rule()
  }
  res <- unlist(res, recursive = FALSE)
  purrr::list_rbind(purrr::map(res, as.data.frame))
}
