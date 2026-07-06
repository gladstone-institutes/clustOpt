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
    cli::cli_alert_info("Splitting PCs into 2 sets odd and even")
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
#' @param verbose Integer verbosity level (0 = silent, 1 = milestones,
#' 2 = detailed, 3 = includes Seurat output)
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

      # Drop any pre-existing SCT assay so SCTransform doesn't warn about
      # mismatched cells/features when overwriting one built on a different
      # cell subset
      Seurat::DefaultAssay(train_seurat) <- "RNA"
      train_seurat <- Seurat::DietSeurat(train_seurat, assays = "RNA")

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

      # Drop any pre-existing SCT assay so SCTransform doesn't warn about
      # mismatched cells/features when overwriting one built on a different
      # cell subset
      Seurat::DefaultAssay(train_seurat) <- "RNA"
      train_seurat <- Seurat::DietSeurat(train_seurat, assays = "RNA")

      # Normalize the training subjects
      train_seurat <- Seurat::SCTransform(train_seurat,
        assay = "RNA",
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
#' @param verbose Integer verbosity level (0 = silent, 1 = milestones,
#' 2 = detailed, 3 = includes Seurat output)
#' @param residual_features Character vector of features to restrict the SCT
#' residual computation to (passed to \code{SCTransform(residual.features=)}).
#' Only used for scRNA data. Default \code{NULL} preserves prior behavior
#' (residuals computed for all genes).
#' @return Training data formatted for sil_score, format depends on dtype
#'
#' @export
#' @import SeuratObject
#' @importFrom Seurat SCTransform DefaultAssay DietSeurat Idents WhichCells
prep_test <- function(input,
                      subject_ids,
                      dtype = "scRNA",
                      test_id,
                      verbose = 0,
                      residual_features = NULL) {
  if (dtype == "scRNA") {
    Seurat::Idents(input) <- subject_ids
    test_cells <- Seurat::WhichCells(object = input, idents = test_id)
    test_seurat <- subset(input, cells = test_cells)

    # Drop any pre-existing SCT assay so SCTransform doesn't warn about
    # mismatched cells/features when overwriting one built on a different
    # cell subset
    Seurat::DefaultAssay(test_seurat) <- "RNA"
    test_seurat <- Seurat::DietSeurat(test_seurat, assays = "RNA")

    if (!is.null(residual_features)) {
      # Restrict SCT residuals to the requested features (typically the
      # training variable features that project_pca() actually uses),
      # intersecting first to avoid asking SCTransform for absent genes.
      residual_features <- intersect(residual_features, rownames(test_seurat))
      test_seurat <- Seurat::SCTransform(test_seurat,
        assay = "RNA",
        verbose = (verbose >= 3),
        residual.features = residual_features,
        return.only.var.genes = FALSE,
        min_cells = 1
      )
    } else {
      test_seurat <- Seurat::SCTransform(test_seurat,
        assay = "RNA",
        verbose = (verbose >= 3),
        variable.features.n = length(rownames(test_seurat)),
        return.only.var.genes = FALSE,
        min_cells = 1
      )
    }
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
#' @param compute_train_eval Logical; if \code{TRUE}, also compute the
#'   training data projected onto clustering PCs (\code{train_proj_clust_pcs}).
#'   Default is \code{FALSE} because this projection is not used in the
#'   standard \code{\link{clust_opt}} pipeline.
#'
#' @details
#' Identifies features that are common between the training
#' and test data sets. Extracts the PCA loadings from the training data
#' for common features. Both training and test data are projected onto these
#' loadings. Data is projected onto loadings for 2 PC sets (odd and even)
#'
#'
#' @return
#' A list containing projected data frames/matrices:
#' \describe{
#'   \item{train_proj_train_with_pcs}{Training data projected onto training PCs (data.frame)}
#'   \item{test_proj_train_with_pcs}{Test data projected onto training PCs (data.frame)}
#'   \item{test_proj_clust_pcs}{Test data projected onto clustering PCs (matrix)}
#'   \item{train_proj_clust_pcs}{Training data projected onto clustering PCs (data.frame); only present when \code{compute_train_eval = TRUE}}
#' }
#'
#'
#' @export
#'
#' @importFrom Seurat VariableFeatures Loadings GetAssayData
project_pca <- function(train_seurat,
                        test_seurat,
                        train_with_pcs,
                        clust_pcs,
                        dtype,
                        verbose = 0,
                        compute_train_eval = FALSE) {
  if (clust_pcs == train_with_pcs) {
    stop("clust_pcs and train_with_pcs must be independent")
  }

  # Validate input
  if (!inherits(train_seurat, "Seurat") || !inherits(test_seurat, "Seurat")) {
    stop("Both train_seurat and test_seurat must be Seurat objects")
  }

  if (verbose >= 2) {
    cli::cli_alert_info("Training with data projected onto {.field {train_with_pcs}}")
  }

  assay_id <- switch(dtype,
    scRNA = "SCT",
    CyTOF = "RNA"
  )
  # Extract scale.data ONCE per Seurat object
  train_scale_full <- as.matrix(Seurat::GetAssayData(
    train_seurat, assay = assay_id, layer = "scale.data"
  ))
  test_scale_full <- as.matrix(Seurat::GetAssayData(
    test_seurat, assay = assay_id, layer = "scale.data"
  ))

  # Features present in both the training variable features and test sample
  common_features <- base::intersect(
    rownames(test_scale_full),
    Seurat::VariableFeatures(train_seurat)
  )
  n_shared_genes <- length(common_features)
  total_genes <- nrow(train_scale_full)
  if (verbose >= 2) {
    pct <- round((n_shared_genes / total_genes) * 100, 2)
    cli::cli_alert_info(
      "Found {n_shared_genes} ({pct}%) shared genes for projecting test data"
    )
  }

  if ((n_shared_genes / total_genes) * 100 < 80) {
    cli::cli_warn("Less than 80% of genes available for projection.")
  }



  # Subset to common features
  train_scale <- train_scale_full[common_features, ]
  test_scale <- test_scale_full[common_features, ]
  rm(train_scale_full, test_scale_full)

  # Extract loadings ONCE per reduction, using direct matrix subsetting
  loadings_train <- Seurat::Loadings(
    train_seurat[[train_with_pcs]]
  )[common_features, ]
  loadings_eval <- Seurat::Loadings(
    train_seurat[[clust_pcs]]
  )[common_features, ]

  if (verbose >= 2) {
    cli::cli_alert_info("Evaluating test data projected onto {.field {clust_pcs}}")
  }

  # Projections using crossprod (avoids materializing transposed matrices)
  result <- list(
    train_proj_train_with_pcs = as.data.frame(crossprod(train_scale, loadings_train)),
    test_proj_train_with_pcs  = as.data.frame(crossprod(test_scale, loadings_train)),
    test_proj_clust_pcs       = crossprod(test_scale, loadings_eval)
  )

  if (compute_train_eval) {
    result$train_proj_clust_pcs <- as.data.frame(crossprod(train_scale, loadings_eval))
  }

  result
}

#' @title train_random_forest
#' @description
#' Train the random forest and predict on the test sample
#' @param res Resolution to train on
#' @param df_list List containing training and test data
#' @param cluster_by_res Named list of cluster assignment vectors, keyed by
#'   resolution (as character). Pre-extracted from training metadata using
#'   exact column names.
#' @param sam Test sample
#' @param num_trees Number of trees for the random forest
#' @param verbose Integer verbosity level (0 = silent, 1 = milestones,
#' 2 = detailed, 3 = includes Seurat output)
#' @param snn_graph Precomputed SNN graph (sparse matrix). Required.
#' @param precomputed_dist Optional precomputed distance matrix (output of
#'   \code{\link[stats]{dist}}). Passed through to
#'   \code{\link{calculate_silhouette_score}}.
#' @param rf_num_threads Number of threads for ranger (default 1).
#' @return A list containing the resolution, silhouette score, and number of
#' predicted clusters.
#'
#' @export
#' @importFrom ranger ranger predictions
#' @importFrom stats predict
train_random_forest <- function(res, df_list, cluster_by_res,
                                sam, num_trees, verbose = 0,
                                snn_graph = NULL,
                                precomputed_dist = NULL,
                                rf_num_threads = 1) {
  # Cluster assignments for this res (already a factor from Seurat metadata)
  train_clusters <- as.factor(cluster_by_res[[as.character(res)]])

  # proportion of clusters in the training assignments
  qx <- data.frame(base::table(train_clusters))
  colnames(qx) <- c("clusters", "Freq_q")
  # Train model on the projected training PCs. The x/y interface avoids the
  # formula/model.frame overhead and the per-call data.frame copy.
  rf <- ranger::ranger(
    x = df_list[["train_proj_train_with_pcs"]],
    y = train_clusters,
    num.trees = num_trees,
    write.forest = TRUE,
    num.threads = rf_num_threads
  )

  # Predict on the hold out sample
  predicted <- stats::predict(rf, df_list[["test_proj_train_with_pcs"]])
  predicted <- ranger::predictions(predicted)
  predicted_char <- as.character(predicted)
  predicted_clusters_table <- base::table(predicted)

  px <- data.frame(predicted_clusters_table)
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
    df_list[["test_proj_clust_pcs"]],
    precomputed_dist = precomputed_dist
  )

  mse_value <- calculate_mse_score(
    predicted,
    df_list[["test_proj_clust_pcs"]]
  )

  if (is.null(snn_graph)) {
    stop("snn_graph must be provided (precomputed via FindNeighbors)")
  }
  modularity_value <- calculate_modularity(snn_graph, predicted)

  list(
    resolution = res,
    test_sample = sam,
    avg_width = sil$avg_width,
    cluster_median_widths = sil$group_median_width,
    n_predicted_clusters = length(unique(predicted_char)),
    min_predicted_cell_per_cluster = min(predicted_clusters_table),
    max_predicted_cell_per_cluster = max(predicted_clusters_table),
    mse = mse_value$mse,
    mad = mse_value$mad,
    KLdivergence = KLdivergence,
    Hellinger = Hellinger,
    modularity = modularity_value
  )
}
