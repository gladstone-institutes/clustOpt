#' @include metrics.R
#' @importFrom rlang .data
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data Preparation Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title split_pca_dimensions
#' @description
#' Takes a Seurat object with an existing PCA reduction and splits the
#' dimensions into 2 sets odd and even PCs.
#'
#' @param input Seurat object
#' @param verbose Integer verbosity level (0 = silent, 1 = milestones,
#' 2 = detailed, 3 = includes Seurat output)
#'
#' @return Seurat object with new PCA reductions
#' @keywords internal
#' @import Seurat
split_pca_dimensions <- function(input,
                                 verbose = 0) {
  if (!inherits(input, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (!"pca" %in% names(input@reductions)) {
    stop("Seurat object must have a PCA reduction")
  }


  if (verbose >= 2) {
    message("Splitting PCs into 2 sets odd and even")
  }
  pca <- input@reductions$pca
  dims <- ncol(pca@cell.embeddings)

  even_dims <- seq(2, dims, by = 2)
  odd_dims <- seq(1, dims, by = 2)

  even_pca <- pca
  even_pca@cell.embeddings <- even_pca@cell.embeddings[,
    even_dims,
    drop = FALSE
  ]
  even_pca@feature.loadings <- even_pca@feature.loadings[,
    even_dims,
    drop = FALSE
  ]
  even_pca@stdev <- even_pca@stdev[even_dims]
  even_pca@key <- "even_pca_"

  odd_pca <- pca
  odd_pca@cell.embeddings <- odd_pca@cell.embeddings[, odd_dims, drop = FALSE]
  odd_pca@feature.loadings <- odd_pca@feature.loadings[, odd_dims, drop = FALSE]
  odd_pca@stdev <- odd_pca@stdev[odd_dims]
  odd_pca@key <- "odd_pca_"

  input@reductions$even_pca <- even_pca
  input@reductions$odd_pca <- odd_pca
  return(input)
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
                       test_id,
                       verbose = 0) {
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
        stop(paste0("More than one batch found for this subject: ", test_id))
      }

      train_cells <- input@meta.data |>
        dplyr::filter(get(subject_ids) != test_id &
          get(within_batch) == this_batch) |>
        rownames()
      train_seurat <- subset(input, cells = train_cells)

      # Normalize the training subjects
      train_seurat <- Seurat::SCTransform(train_seurat,
        assay = "RNA",
        verbose = (verbose >= 3)
      )
      train_seurat <- Seurat::DietSeurat(train_seurat, assays = "SCT")
      return(train_seurat)
    } else {
      # Return all other subjects
      Seurat::Idents(input) <- subject_ids
      train_cells <- Seurat::WhichCells(
        object = input,
        idents = test_id,
        invert = TRUE
      )
      train_seurat <- subset(input, cells = train_cells)
      # Normalize the training subjects
      train_seurat <- Seurat::SCTransform(train_seurat,
        assay = Seurat::DefaultAssay(train_seurat),
        verbose = (verbose >= 3)
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
      # Return all other subjects
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
#' @param subject_ids Metadata field that identifies unique subjects.
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
                      test_id,
                      verbose = 0) {
  if (dtype == "scRNA") {
    Seurat::Idents(input) <- subject_ids
    test_cells <- Seurat::WhichCells(object = input, idents = test_id)
    test_seurat <- subset(input, cells = test_cells)

    test_seurat <- Seurat::SCTransform(test_seurat,
      assay = "RNA",
      verbose = (verbose >= 3),
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
#' @param verbose Integer verbosity level (0 = silent, 1 = milestones,
#' 2 = detailed, 3 = includes Seurat output)
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
#' @importFrom Seurat VariableFeatures Loadings GetAssayData
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter
project_pca <- function(train_seurat,
                        test_seurat,
                        train_with_pcs,
                        clust_pcs,
                        dtype,
                        verbose) {
  if (clust_pcs == train_with_pcs) {
    stop("clust_pcs and train_with_pcs must be independent")
  }

  # Validate input
  if (!inherits(train_seurat, "Seurat") || !inherits(test_seurat, "Seurat")) {
    stop("Both train_seurat and test_seurat must be Seurat objects")
  }

  if (verbose >= 2) {
    message(sprintf("Training with data projected onto %s", train_with_pcs))
  }

  assay_id <- switch(dtype,
    scRNA = "SCT",
    CyTOF = "RNA"
  )
  # Features present in both the training variable features and test sample
  common_features <- base::intersect(
    rownames(Seurat::GetAssayData(test_seurat,
      assay = assay_id,
      layer = "scale.data"
    )),
    Seurat::VariableFeatures(train_seurat)
  )
  n_shared_genes <- length(common_features)
  total_genes <- length(rownames(Seurat::GetAssayData(train_seurat,
    assay = assay_id,
    layer = "scale.data"
  )))
  if (verbose >= 2) {
    message(sprintf(
      "Found %d (%.2f%%) shared genes used for projecting test data",
      n_shared_genes,
      (n_shared_genes / total_genes) * 100
    ))
  }

  if ((n_shared_genes / total_genes) * 100 < 80) {
    warning("Less than 80% of genes available for projection.")
  }



  project_data_train <- function(seurat_obj, assay_id) {
    loadings_common_features <- Seurat::Loadings(
      train_seurat[[train_with_pcs]]
    ) |>
      tibble::as_tibble(rownames = "features") |>
      dplyr::filter(.data$features %in% common_features) |>
      as.matrix()

    rownames(loadings_common_features) <- loadings_common_features[, 1]
    loadings_common_features <- loadings_common_features[, -1]
    class(loadings_common_features) <- "numeric"
    scale_data <- as.matrix(Seurat::GetAssayData(seurat_obj,
      assay = assay_id,
      layer = "scale.data"
    ))[common_features, ]

    t(scale_data) %*% loadings_common_features
  }
  if (verbose >= 2) {
    message(sprintf("Evaluating test data projected onto %s", clust_pcs))
  }
  project_data_for_eval <- function(seurat_obj, assay_id) {
    loadings_common_features <- Seurat::Loadings(
      train_seurat[[clust_pcs]]
    ) |>
      tibble::as_tibble(rownames = "features") |>
      dplyr::filter(.data$features %in% common_features) |>
      as.matrix()

    rownames(loadings_common_features) <- loadings_common_features[, 1]
    loadings_common_features <- loadings_common_features[, -1]
    class(loadings_common_features) <- "numeric"
    scale_data <- as.matrix(Seurat::GetAssayData(seurat_obj,
      assay = assay_id,
      layer = "scale.data"
    ))[common_features, ]

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
    train_proj_clust_pcs = pca_train_eval_data |>
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

  # proportion of clusters in train_df
  qx <- data.frame(base::table(train_df$clusters))
  colnames(qx) <- c("clusters", "Freq_q")
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

  px <- data.frame(base::table(predicted))
  colnames(px) <- c("clusters", "Freq_p")

  probs <- base::merge(qx, px, all = TRUE)
  # assign count of 1 to clusters missing in the predicted test set, so that KL divergence does not blow up
  probs[is.na(probs)] <- 0
  probs[,2:3] <- probs[,2:3] + 1
  probs[,2] <- probs[,2]/base::sum(probs[,2])
  probs[,3] <- probs[,3]/base::sum(probs[,3])

  rm(rf)

  # Calculate distribution divergence metrics
  KLdivergence <- calculate_kl_divergence(probs$Freq_q, probs$Freq_p)
  Hellinger <- calculate_hellinger_distance(probs$Freq_q, probs$Freq_p)
  # Evaluate clustering on data project on to the opposite PCs
  sil <- calculate_silhouette_score(
    predicted,
    df_list[["test_proj_clust_pcs"]]
  )

  mse_value <- calculate_mse_score(
    predicted,
    df_list[["test_proj_clust_pcs"]]
  )

  modularity_value <- calculate_modularity_from_coords(
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
    max_predicted_cell_per_cluster = max(predicted_clusters_table),
    mse = mse_value$mse,
    mad = mse_value$mad,
    KLdivergence = KLdivergence,
    Hellinger = Hellinger,
    modularity = modularity_value$modularity
  )
}
