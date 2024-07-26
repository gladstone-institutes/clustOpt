#' @include utils.R
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
#'  slot) and sketching is not implemented for CyTOF.
#' @param sketch_size Number of cells to use for sketching.
#' @param subject_ids Metadata field that identifies unique samples.
#' @param res_range Range of resolutions to test.
#' @param verbose print messages.
#' @param within_batch Batch variable, for a given sample only those with the
#' same value for the batch variable will be used for training.
#' @param num.trees Number of trees to use in the random forest.
#' @return A data.frame containing a distribution of silhouette scores for each
#' resolution.
#'
#' @export
#'
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
                      num_trees = 1000) {
  sample_names <- as.vector((unique(input@meta.data[[subject_ids]])))

  if (!(dtype %in% c("CyTOF", "scRNA"))) {
    stop("dtype is not one of 'CyTOF' or 'scRNA'")
  }

  if (dtype == "scRNA") {
    # Clear any previous normalizations
    Seurat::DefaultAssay(input) <- "RNA"
    input <- Seurat::DietSeurat(input, assays = "RNA")
    if (!skip_sketch) {
      if (check_size(input) || !is.null(sketch_size)) {
        message("Sketching input data")
        input <- leverage_sketch(input, sketch_size)
      } else {
        message("Input is small enough to run with all cells")
      }
    }
  }

  # Make sure a seed is set
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    t <- as.integer(Sys.time())
    message("Setting seed: ", t)
    set.seed(t)
  }

  # Get every combination of test sample and resolution
  runs <- expand.grid(sample_names, res_range)
  message(paste0(
    "Found ", nrow(runs),
    " combinations of test sample and resolution"
  ))

  # Set up progress logging
  progressr::handlers("progress")
  p <- progressr::progressor(along = unique(runs[, 1]))

  res <- NULL
  for (sam in unique(runs[, 1])) {
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
    if (dtype == "scRNA") {
      message(sprintf(
        "Found %d (%.2f%%) shared genes between testing and training data",
        length(intersect(
          rownames(test@assays[["SCT"]]@scale.data),
          Seurat::VariableFeatures(train)
        )),
        length(intersect(
          rownames(test@assays[["SCT"]]@scale.data),
          Seurat::VariableFeatures(train)
        )) / length(rownames(train@assays[["SCT"]]@scale.data)) * 100
      ))
    }
    # Add switch here to select odd or even
    if (dtype == "scRNA") {
      train <- Seurat::RunPCA(train,
        npcs = ndim,
        verbose = FALSE,
        assay = "SCT"
      )
      # Separate even or odd PCA dimensions
      # 2 reductions
      # Even PCA and Odd PCA
      
      train <- Seurat::FindNeighbors(
        object = train,
        dims = 1:ndim,
        verbose = verbose
        # reduction = "[switch]_pca"
      )

      train <- Seurat::FindClusters(
        object = train,
        resolution = res_range,
        verbose = verbose
        # reduction = "[switch]_pca"
      )
    } else {
      train <- Seurat::ScaleData(train, features = NULL, verbose = verbose)
      train <- Seurat::FindVariableFeatures(train,
        selection.method = "vst", nfeatures = ndim
      )
      train <- Seurat::RunPCA(train,
        npcs = ndim, approx = FALSE,
        verbose = verbose
      )

      train <- Seurat::FindNeighbors(
        object = train, dims = 1:ndim,
        verbose = verbose
        # reduction = "[switch]_pca"
      )
      train <- Seurat::FindClusters(
        object = train, resolution = res_range,
        verbose = verbose
        # reduction = "[switch]_pca"
      )
    }
    if (verbose) {
      message("Clustering complete..")
    }

    if (dtype == "scRNA") {
      train_clusters <- train@meta.data |>
        dplyr::select(dplyr::contains("SCT_snn_res"))
      
      if (verbose) {
        message("Project the cells in the test data onto the train PCs..")
      }
      df_list <- project_pca(train, test, ndim) # Add switch here to select opposite of PCs that were used for clustering
      rm(train, test)
    } else {
      train_clusters <- train@meta.data |>
        dplyr::select(dplyr::contains("RNA_snn_res"))

      df_list <- prepare_cytof(train, test) # Add switch here to select opposite of PCs that were used for clustering
      rm(train, test)
    }

    this_result <- future.apply::future_lapply(
<<<<<<< HEAD
      seq_len(nrow(runs)),
||||||| d571bce
      1:nrow(runs),
=======
      rownames(runs[runs$Var1 == sam,]),
>>>>>>> main
      function(i) {
        train_random_forest(
          res = runs[i, 2],
          df_list = df_list,
          train_clusters = train_clusters,
          sam = runs[i, 1],
          num.trees = num_trees
          # Add switch here to select odd or even
        )
      },
      future.seed = TRUE
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
        stop("More than one batch found for this sample")
      }

      train_cells <- input@meta.data |>
        dplyr::filter(get(subject_ids) != test_id &
          get(within_batch) == this_batch) |>
        rownames()
      train_seurat <- subset(input, cells = train_cells)

      # Normalize the training samples
      train_seurat <- Seurat::SCTransform(train_seurat,
        assay = Seurat::DefaultAssay(train_seurat),
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
#'
prep_test <- function(input,
                      subject_ids,
                      dtype = "scRNA",
                      test_id) {
  if (dtype == "scRNA") {
    Seurat::Idents(input) <- subject_ids
    test_cells <- Seurat::WhichCells(object = input, idents = test_id)
    test_seurat <- subset(input, cells = test_cells)

    test_seurat <- Seurat::SCTransform(test_seurat,
      assay = Seurat::DefaultAssay(test_seurat),
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
#' @param ndim The number of principal components to use for projection.
#'
#' @details
#' Identifies features that are common between the training
#' and test data sets. It then extracts the PCA loadings from the training data
#' for these common features. Both training and test data are projected onto
#' these loadings. Finally, the function selects the specified number of
#' dimensions (principal components) for the output.
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
project_pca <- function(train_seurat, test_seurat, ndim) {
  # Validate input
  if (!("Seurat" %in% class(train_seurat)) ||
    !("Seurat" %in% class(test_seurat))) {
    stop("Both train_seurat and test_seurat must be Seurat objects")
  }
  if (!is.numeric(ndim) || ndim <= 0) {
    stop("ndim must be a positive integer")
  }

  # Features present in both the training variable features and test sample
  common_features <- base::intersect(
    rownames(test_seurat@assays[["SCT"]]@scale.data),
    Seurat::VariableFeatures(train_seurat)
  )

  # Extract the loadings for common features from the training data
  # ** Project the cells in the training and test data onto the opposite PCs (if even is used for clustering, use odd PCs)
  # Subset to ndim, then select odd or even
  loadings_common_features <- Seurat::Loadings(train_seurat[["pca"]]) |>
    tibble::as_tibble(rownames = "features") |>
    dplyr::filter(features %in% common_features) |> # nolint
    as.matrix()

  rownames(loadings_common_features) <- loadings_common_features[, 1]
  loadings_common_features <- loadings_common_features[, -1]
  class(loadings_common_features) <- "numeric"

  # Function to project data
  project_data <- function(seurat_obj) {
    scale_data <- as.matrix(seurat_obj[["SCT"]]@scale.data)[common_features, ]
    t(scale_data) %*% loadings_common_features
  }

  # Project the cells in the training and test data onto the opposite PCs (if even is used for clustering, use odd PCs)
  pca_train_data <- project_data(train_seurat)
  pca_test_data <- project_data(test_seurat)

  rm(loadings_common_features)


  # Select the number of dimensions to train with
  list(
    train_df = pca_train_data[, 1:ndim, drop = FALSE] |>
      as.data.frame(),
    test_df = pca_test_data[, 1:ndim, drop = FALSE] |>
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
#' @return A list containing the resolution, silhouette score, and number of
#' predicted clusters.
#' @export
#'
train_random_forest <- function(res, df_list, train_clusters, sam, num.trees) {
  # Get cluster assignments for this res
  train_df <- df_list[[1]] |>
    dplyr::mutate(clusters = train_clusters |>
      dplyr::select(dplyr::contains(as.character(res))) |>
      dplyr::pull())
  
  # Train model
  rf <- ranger::ranger(as.factor(clusters) ~ .,
    data = train_df,
    num.trees = num.trees,
    write.forest = TRUE,
    num.threads = 1
  )
  rm(train_df)

  # Predict on the hold out sample
  predicted <- stats::predict(rf, df_list[[2]])
  predicted <- ranger::predictions(predicted)
  predicted_clusters_table <- base::table(predicted)
  rm(rf)
  sil <- calculate_silhouette_score(predicted, df_list[[2]]) # Replace df_list[[2]] with an independent PCA using total dimensions totally independent of the training

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
#' @export
#'
calculate_silhouette_score <- function(predicted, data_frame) {
  sil <- cluster::silhouette(
    as.numeric(as.character(predicted)),
    dist(data_frame)
  )

  # Check if there is only one cluster
  if (class(sil) == "logical") {
    return(list(avg_width = NA, group_median_width = NA))
  } else {
    sil_summary <- summary(sil)
    return(list(
      avg_width = sil_summary$avg.width,
      group_median_width = median(sil_summary$clus.avg.widths)
    ))
  }
}

#' @title prepare_cytof
#' @description
#' Prepare CyTOF data for random forest
#'
#' @param train_seurat A Seurat object representing the training data set.
#' @param test_seurat A Seurat object representing the test data set.
#'
#' @return A list of the training and test data formatted for RF
#' @export
#'
prepare_cytof <- function(train_seurat, test_seurat) {
  list(
    Seurat::GetAssayData(train_seurat, slot = "counts", assay = "RNA") |>
      as.data.frame() |>
      t() |>
      tibble::as_tibble(rownames = "CellID") |>
      dplyr::select(CellID, everything()),
    Seurat::GetAssayData(test_seurat, slot = "counts", assay = "RNA") |>
      as.data.frame() |>
      t() |>
      tibble::as_tibble(rownames = "CellID") |>
      dplyr::select(CellID, everything())
  )
}
