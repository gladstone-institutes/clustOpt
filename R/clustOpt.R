#' @include utils.R pca_split.R
#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title clust_opt
#' @description
#' Runs the main resolution optimization algorithm
#'
#' @param input Seurat object
#' @param ndim Number of principal components to use.
#' @param dtype Type of data in the Seurat object "scRNA" or "CyTOF", default
#' is "scRNA". CyTOF data is expected to be arcsinh normalized  (in the counts
#' slot) and sketching is not implemented for CyTOF.
#' @param sketch_size Number of cells to use for sketching.
#' @param skip_sketch Skip sketching, by default any input with more than
#' 200,000 cells is sketched to 10\% of the cells.
#' @param subject_ids Metadata field that identifies unique samples.
#' @param res_range Range of resolutions to test.
#' @param verbose Output messages.
#' @param within_batch Batch variable, for a given sample only those with the
#' same value for the batch variable will be used for training.
#' @param num_trees Number of trees to use in the random forest.
#' @param train_with Either odd or even PCs for clustering and training.
#' Default is odd.
#' @param min_cells Minimum cells per subject, default is 50
#' @return A data.frame containing a distribution of silhouette scores for each
#' resolution.
#'
#' @export
#' @importFrom Seurat DefaultAssay DietSeurat RunPCA FindNeighbors FindClusters
#' @importFrom Seurat ScaleData FindVariableFeatures
#' @importFrom dplyr select contains
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
                      verbose = FALSE,
                      num_trees = 1000,
                      train_with = "odd",
                      min_cells = 50) {
  # Make sure a seed is set, only setting if the user has not set one
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    t <- as.integer(Sys.time())
    message("Setting seed: ", t)
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
      message("Sketching input data")
      input <- leverage_sketch(input, sketch_size, dtype)
    } else {
      message("Input is small enough to run with all cells")
    }
  }


  sample_names <- get_valid_samples(input, subject_ids, min_cells)
  if (is.null(sample_names)) {
    stop(paste0(
      "Unable to perform cluster resolution optimization for this data. ",
      "There are less than 3 subjects with at least", min_cells, "cells."
    ))
  }

  # Get every combination of test sample and resolution
  runs <- expand.grid(sample_names, res_range)
  message(paste0(
    "Found ", nrow(runs),
    " combinations of test subject and resolution"
  ))

  # Set up progress logging
  progressr::handlers("progress")
  p <- progressr::progressor(along = unique(runs[, 1]))

  res <- NULL
  for (sam in unique(runs[, 1])) {
    message(paste0("Holdout subject: ", sam))
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

    if (dtype == "scRNA") {
      train <- Seurat::RunPCA(train,
        npcs = ndim,
        verbose = FALSE,
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

      if (verbose) {
        message(sprintf("Clustering with %s", clust_pcs))
      }

      train <- Seurat::FindNeighbors(
        object = train,
        dims = seq_len(ncol(train@reductions[[clust_pcs]]@cell.embeddings)),
        verbose = FALSE,
        reduction = clust_pcs
      )

      train <- Seurat::FindClusters(
        object = train,
        resolution = res_range,
        verbose = FALSE
      )
    } else {
      train <- Seurat::ScaleData(train, features = NULL, verbose = verbose)
      train <- Seurat::FindVariableFeatures(train,
        selection.method = "vst", nfeatures = ndim
      )
      train <- Seurat::RunPCA(train,
        npcs = ndim,
        approx = FALSE,
        verbose = verbose
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
      train <- Seurat::FindNeighbors(
        object = train,
        dims = seq_len(ncol(train@reductions[[clust_pcs]]@cell.embeddings)),
        verbose = FALSE,
        reduction = clust_pcs
      )

      train <- Seurat::FindClusters(
        object = train,
        resolution = res_range,
        verbose = FALSE
      )
    }
    if (verbose) {
      message("Clustering complete..")
    }

    if (dtype == "scRNA") {
      train_clusters <- train@meta.data |>
        dplyr::select(dplyr::contains("SCT_snn_res"))


      df_list <- project_pca(train_seurat = train,
                             test_seurat = test,
                             train_with_pcs = train_with_pcs,
                             clust_pcs = clust_pcs,
                             dtype = dtype,
                             verbose = verbose
      )

      rm(train, test)
    } else {
      train_clusters <- train@meta.data |>
        dplyr::select(dplyr::contains("RNA_snn_res"))

      df_list <- project_pca(train_seurat = train,
                             test_seurat = test,
        train_with_pcs = train_with_pcs,
        clust_pcs = clust_pcs,
        dtype = dtype,
        verbose = verbose
      )
      rm(train, test)
    }

    this_result <- future.apply::future_lapply(
      rownames(runs[runs$Var1 == sam, ]),
      function(i) {
        train_random_forest(
          res = runs[i, 2],
          df_list = df_list,
          train_clusters = train_clusters,
          sam = runs[i, 1],
          num_trees = num_trees
        )
      },
      future.seed = TRUE,
      future.packages = c("SeuratObject", "Seurat", "ranger", "cluster", "dplyr")
    )
    res <- c(res, this_result)
    p()
  }

  purrr::map_df(res, .f = as.data.frame)
}


#' @title prep_train
#' @description
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
#' @importFrom dplyr filter pull
#' @import SeuratObject
#' @importFrom Seurat SCTransform DefaultAssay DietSeurat Idents WhichCells
prep_train <- function(input,
                       subject_ids,
                       dtype = "scRNA",
                       within_batch = NA,
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
        stop(paste0("More than one batch found for this sample: ", test_id))
      }

      train_cells <- input@meta.data |>
        dplyr::filter(get(subject_ids) != test_id &
          get(within_batch) == this_batch) |>
        rownames()
      train_seurat <- subset(input, cells = train_cells)

      # Normalize the training samples
      train_seurat <- Seurat::SCTransform(train_seurat,
        assay = "RNA",
        verbose = FALSE
      )
      train_seurat <- Seurat::DietSeurat(train_seurat, assays = "SCT")
      return(train_seurat)
    } else {
      # Return all other samples
      Seurat::Idents(input) <- subject_ids
      train_cells <- Seurat::WhichCells(
        object = input,
        idents = test_id,
        invert = TRUE
      )
      train_seurat <- subset(input, cells = train_cells)
      # Normalize the training samples
      train_seurat <- Seurat::SCTransform(train_seurat,
        assay = Seurat::DefaultAssay(train_seurat),
        verbose = FALSE
      )
      train_seurat <- Seurat::DietSeurat(train_seurat, assays = "SCT")
      return(train_seurat)
    }
  } else {
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
      return(train_seurat)
    } else {
      # Return all other samples
      Seurat::Idents(input) <- subject_ids
      train_cells <- Seurat::WhichCells(
        object = input,
        idents = test_id,
        invert = TRUE
      )
      train_seurat <- subset(input, cells = train_cells)
      return(train_seurat)
    }
  }
}

#' @title prep_test
#' @description
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
#' @import SeuratObject
#' @importFrom Seurat SCTransform DefaultAssay DietSeurat Idents WhichCells
prep_test <- function(input,
                      subject_ids,
                      dtype = "scRNA",
                      test_id) {
  if (dtype == "scRNA") {
    Seurat::Idents(input) <- subject_ids
    test_cells <- Seurat::WhichCells(object = input, idents = test_id)
    test_seurat <- subset(input, cells = test_cells)

    test_seurat <- Seurat::SCTransform(test_seurat,
      assay = "RNA",
      verbose = FALSE,
      variable.features.n = length(rownames(test_seurat)),
      return.only.var.genes = FALSE,
      min_cells = 1
    )
    test_seurat <- Seurat::DietSeurat(test_seurat, assays = "SCT")

    return(test_seurat)
  } else {
    Seurat::Idents(input) <- subject_ids
    test_cells <- Seurat::WhichCells(object = input, idents = test_id)
    test_seurat <- subset(input, cells = test_cells)
    return(test_seurat)
  }
}

#' @title project_pca
#' @description
#' Project Training and Test Seurat Objects onto Principal Components
#'
#' This function projects both training and test Seurat objects onto a set of
#' principal components derived from the training data.
#'
#' @param train_seurat A Seurat object representing the training data set.
#' @param test_seurat A Seurat object representing the test data set.
#' @param train_with_pcs Which reduction should be used for training 
#' "odd_pca" or "even_pca"
#' @param clust_pcs Which reduction was used for clustering 
#' @param dtype Type of data in the Seurat object "scRNA" or "CyTOF"
#' @param verbose print messages
#'
#' @details
#' Identifies features that are common between the training
#' and test data sets. Extracts the PCA loadings from the training data
#' for common features. Both training and test data are projected onto these
#' loadings. Data is projected onto loadings for 2 PC sets (odd and even)
#'
#'
#' @return
#' A list containing two elements: `train_df` and `test_df`. Each is a matrix
#' of the projected data for the training and test sets, respectively, using
#' the specified number of principal components.
#'
#'
#' @export
#'
#' @importFrom Seurat VariableFeatures Loadings
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter
project_pca <- function(train_seurat,
                        test_seurat,
                        train_with_pcs,
                        clust_pcs,
                        dtype,
                        verbose) {
  
  if(clust_pcs == train_with_pcs) {
    stop("clust_pcs and train_with_pcs must be independent")
  }
  
  # Validate input
  if (!inherits(train_seurat, "Seurat") || !inherits(test_seurat, "Seurat")) {
    stop("Both train_seurat and test_seurat must be Seurat objects")
  }

  if (verbose) {
    message(sprintf("Training with data projected onto %s", train_with_pcs))
  }

  assay_id <- switch(dtype,
    scRNA = "SCT",
    CyTOF = "RNA"
  )
  # Features present in both the training variable features and test sample
  common_features <- base::intersect(
    rownames(test_seurat@assays[[assay_id]]@scale.data),
    Seurat::VariableFeatures(train_seurat)
  )
  n_shared_genes <- length(common_features)
  total_genes <- length(rownames(train_seurat@assays[[assay_id]]@scale.data))
  message(sprintf(
    "Found %d (%.2f%%) shared genes used for projecting test data",
    n_shared_genes,
    (n_shared_genes / total_genes) * 100
  ))



  project_data_train <- function(seurat_obj, assay_id) {
    loadings_common_features <- Seurat::Loadings(
      train_seurat[[train_with_pcs]]
    ) |>
      tibble::as_tibble(rownames = "features") |>
      dplyr::filter(features %in% common_features) |>
      as.matrix()

    rownames(loadings_common_features) <- loadings_common_features[, 1]
    loadings_common_features <- loadings_common_features[, -1]
    class(loadings_common_features) <- "numeric"
    scale_data <- as.matrix(seurat_obj[[assay_id]]@scale.data
                            )[common_features, ]

    t(scale_data) %*% loadings_common_features
  }
  if (verbose) {
    message(sprintf("Evaluating test data projected onto %s", clust_pcs))
  }
  project_data_for_eval <- function(seurat_obj, assay_id) {
    loadings_common_features <- Seurat::Loadings(
      train_seurat[[clust_pcs]]) |>
      tibble::as_tibble(rownames = "features") |>
      dplyr::filter(features %in% common_features) |>
      as.matrix()

    rownames(loadings_common_features) <- loadings_common_features[, 1]
    loadings_common_features <- loadings_common_features[, -1]
    class(loadings_common_features) <- "numeric"
    scale_data <- as.matrix(seurat_obj[[assay_id]]@scale.data)[common_features, ]

    t(scale_data) %*% loadings_common_features
  }

  # Project the cells in the training and test data onto the opposite PCs
  pca_train_data <- project_data_train(train_seurat, assay_id)
  pca_test_data <- project_data_train(test_seurat, assay_id)
  pca_test_eval_data <- project_data_for_eval(test_seurat, assay_id)
  pca_train_eval_data <- project_data_for_eval(train_seurat, assay_id)

  list(
    train_proj_train_with_pcs = pca_train_data |>
      as.data.frame(),
    test_proj_train_with_pcs = pca_test_data |>
      as.data.frame(),
    test_proj_clust_pcs = pca_test_eval_data |>
      as.data.frame(),
    train_proj_clust_pcs  = pca_train_eval_data |>
      as.data.frame()
  )
}

#' @title train_random_forest
#' @description
#' Train the random forest and predict on the test sample
#' @param res Resolution to train on
#' @param df_list List containing training and test data
#' @param train_clusters Cluster assignments for the training data
#' @param sam Test sample
#' @param num_trees Number of trees for the random forest
#' @return A list containing the resolution, silhouette score, and number of
#' predicted clusters.
#'
#' @export
#' @importFrom dplyr mutate select contains pull
#' @importFrom ranger ranger predictions
#' @importFrom stats predict
train_random_forest <- function(res, df_list, train_clusters,
                                sam, num_trees) {
  # Get cluster assignments for this res
  train_df <- df_list[["train_proj_train_with_pcs"]] |>
    dplyr::mutate(clusters = train_clusters |>
      dplyr::select(dplyr::contains(as.character(res))) |>
      dplyr::pull())

  # Train model
  rf <- ranger::ranger(as.factor(clusters) ~ .,
    data = train_df,
    num.trees = num_trees,
    write.forest = TRUE,
    num.threads = 1
  )
  rm(train_df)

  # Predict on the hold out sample
  predicted <- stats::predict(rf, df_list[["test_proj_train_with_pcs"]])
  predicted <- ranger::predictions(predicted)
  predicted_clusters_table <- base::table(predicted)
  rm(rf)
  # Evaluate clustering on data project on to the opposite PCs
  sil <- calculate_silhouette_score(
    predicted,
    df_list[["test_proj_clust_pcs"]]
  )

  list(
    resolution = res,
    test_sample = sam,
    avg_width = sil$avg_width,
    cluster_median_widths = sil$group_median_width,
    n_predicted_clusters = length(unique(as.character(predicted))),
    min_predicted_cell_per_cluster = min(predicted_clusters_table),
    max_predicted_cell_per_cluster = max(predicted_clusters_table)
  )
}
#' @title calculate_silhouette_score
#' @description
#' Calculate silhouette score
#' @param predicted Cluster assignments
#' @param data_frame Data frame containing the data
#' @return A list containing the average silhouette score and the average
#' silhouette score for each cluster.
#'
#' @export
#' @importFrom stats dist
#' @importFrom cluster silhouette
calculate_silhouette_score <- function(predicted, data_frame) {
  # Check if there's only one unique cluster
  unique_clusters <- length(unique(predicted))
  if (unique_clusters <= 1) {
    return(list(avg_width = NA, group_median_width = NA))
  }
  
  # Calculate silhouette scores
  sil <- cluster::silhouette(
    as.numeric(as.character(predicted)),
    dist(data_frame)
  )
  
  sil_summary <- summary(sil)
  return(list(
    avg_width = sil_summary$avg.width,
    group_median_width = median(sil_summary$clus.avg.widths)
  ))
}
