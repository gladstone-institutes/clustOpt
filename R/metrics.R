#' @importFrom stats median var sd dist
#' @importFrom dplyr group_by summarize
#' @importFrom rlang .data
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clustering Quality Metrics Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Distribution Distance Metrics
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title calculate_kl_divergence
#' @description
#' Calculate Kullback-Leibler divergence between two probability distributions.
#' D_KL(Q || P) measures how much information is lost when P is used to
#' approximate Q. The result is always non-negative, with 0 indicating
#' identical distributions.
#'
#' @param q Numeric vector representing the "true" probability distribution.
#'   Must be non-negative and sum to 1.
#' @param p Numeric vector representing the "approximating" probability
#'   distribution. Must be non-negative, sum to 1, and have the same length as q.
#'
#' @return Numeric value representing D_KL(Q || P), always >= 0.
#'   Returns 0 when distributions are identical.
#'
#' @details
#' The KL divergence is calculated as: D_KL(Q || P) = sum(q * log(q / p))
#'
#' Note that KL divergence is asymmetric: D_KL(Q || P) != D_KL(P || Q).
#' When q[i] > 0 but p[i] = 0, the divergence is infinite. This implementation
#' requires all elements of p to be positive when corresponding elements of q
#' are positive.
#'
#' @examples
#' # Identical distributions
#' p <- c(0.25, 0.25, 0.25, 0.25)
#' calculate_kl_divergence(p, p)  # Returns 0
#'
#' # Different distributions
#' q <- c(0.5, 0.3, 0.1, 0.1)
#' p <- c(0.25, 0.25, 0.25, 0.25)
#' calculate_kl_divergence(q, p)  # Returns positive value
#'
#' @export
calculate_kl_divergence <- function(q, p) {
  # Input validation
  if (!is.numeric(q) || !is.numeric(p)) {
    stop("q and p must be numeric vectors")
  }

  if (length(q) != length(p)) {
    stop("q and p must have the same length")
  }

  if (any(q < 0) || any(p < 0)) {
    stop("q and p must contain non-negative values")
  }

  # Check that distributions sum to 1 (with tolerance for floating point)
  if (abs(sum(q) - 1) > 1e-6) {
    stop("q must sum to 1 (current sum: ", sum(q), ")")
  }

  if (abs(sum(p) - 1) > 1e-6) {
    stop("p must sum to 1 (current sum: ", sum(p), ")")
  }

  # Handle case where q[i] > 0 but p[i] = 0 (would give infinite divergence)
  if (any(q > 0 & p == 0)) {
    stop("p must be positive wherever q is positive (KL divergence undefined)")
  }

  # Calculate KL divergence, only for indices where q > 0
  # (0 * log(0/p) = 0 by convention)
  positive_q <- q > 0
  kl <- sum(q[positive_q] * log(q[positive_q] / p[positive_q]))

  return(kl)
}

#' @title calculate_hellinger_distance
#' @description
#' Calculate Hellinger distance between two probability distributions.
#' The Hellinger distance is a symmetric measure bounded between 0 and 1,
#' where 0 indicates identical distributions and 1 indicates distributions
#' with no overlap.
#'
#' @param q Numeric vector representing the first probability distribution.
#'   Must be non-negative and sum to 1.
#' @param p Numeric vector representing the second probability distribution.
#'   Must be non-negative, sum to 1, and have the same length as q.
#'
#' @return Numeric value between 0 and 1 representing the Hellinger distance.
#'
#' @details
#' The Hellinger distance is calculated as:
#' H(P, Q) = sqrt(0.5 * sum((sqrt(p) - sqrt(q))^2))
#'
#' Unlike KL divergence, Hellinger distance is symmetric: H(P, Q) = H(Q, P).
#'
#' @examples
#' # Identical distributions
#' p <- c(0.25, 0.25, 0.25, 0.25)
#' calculate_hellinger_distance(p, p)  # Returns 0
#'
#' # Completely different distributions (no overlap)
#' q <- c(1, 0, 0, 0)
#' p <- c(0, 1, 0, 0)
#' calculate_hellinger_distance(q, p)  # Returns 1
#'
#' @export
calculate_hellinger_distance <- function(q, p) {
  # Input validation
  if (!is.numeric(q) || !is.numeric(p)) {
    stop("q and p must be numeric vectors")
  }

  if (length(q) != length(p)) {
    stop("q and p must have the same length")
  }

  if (any(q < 0) || any(p < 0)) {
    stop("q and p must contain non-negative values")
  }

  # Check that distributions sum to 1 (with tolerance for floating point)
  if (abs(sum(q) - 1) > 1e-6) {
    stop("q must sum to 1 (current sum: ", sum(q), ")")
  }

  if (abs(sum(p) - 1) > 1e-6) {
    stop("p must sum to 1 (current sum: ", sum(p), ")")
  }

  # Calculate Hellinger distance
  hellinger <- sqrt(0.5 * sum((sqrt(q) - sqrt(p))^2))

  return(hellinger)
}

#' @title calculate_silhouette_score
#' @description
#' Calculate silhouette score
#' @param predicted Cluster assignments
#' @param data_frame Data frame containing the data
#' @param precomputed_dist Optional precomputed distance matrix (output of
#'   \code{\link[stats]{dist}}). If \code{NULL} (default), the distance matrix
#'   is computed from \code{data_frame}.
#' @return A list containing the average silhouette score and the average
#' silhouette score for each cluster.
#'
#' @export
#' @importFrom stats dist
#' @importFrom cluster silhouette
calculate_silhouette_score <- function(predicted, data_frame,
                                       precomputed_dist = NULL) {
  # Check if there's only one unique cluster
  unique_clusters <- length(unique(predicted))
  if (unique_clusters <= 1) {
    return(list(avg_width = NA, group_median_width = NA))
  }

  d <- if (!is.null(precomputed_dist)) precomputed_dist else dist(data_frame)

  # Calculate silhouette scores
  sil <- tryCatch({
    cluster::silhouette(
      as.numeric(as.character(predicted)),
      d
    )
  }, error = function(e) {
    warning(sprintf("Silhouette calculation failed: %s", e$message))
    return(NULL)
  })

  if (is.null(sil)) {
    return(list(avg_width = NA, group_median_width = NA))
  }

  sil_summary <- tryCatch(summary(sil), error = function(e) NULL)
  if (is.null(sil_summary) || is.atomic(sil_summary)) {
    return(list(avg_width = NA, group_median_width = NA))
  }
  return(list(
    avg_width = sil_summary$avg.width,
    group_median_width = median(sil_summary$clus.avg.widths)
  ))
}

#' @title calculate_mse_score
#' @description
#' Calculate mean squared error score
#' @param predicted Cluster assignments
#' @param data_frame Data frame containing the data
#' @return A numeric value representing the mse
#'
#' @export
calculate_mse_score <- function(predicted, data_frame) {
  n <- nrow(data_frame)
  data_mat <- as.matrix(data_frame)
  # Use character labels for safe named-row indexing (avoids factor level issues)
  pred_char <- as.character(predicted)
  unique_clusters <- unique(pred_char)

  # MSE: vectorized via rowsum + character-based row indexing
  cluster_sums <- rowsum(data_mat, pred_char)
  cluster_sizes <- as.numeric(table(pred_char))
  centroids <- cluster_sums / cluster_sizes
  centroid_expanded <- centroids[pred_char, , drop = FALSE]
  diff_mat <- data_mat - centroid_expanded
  mse <- sum(diff_mat^2) / n

  # MAD: per-cluster medians (vapply avoids intermediate list allocation)
  cluster_medians <- t(vapply(unique_clusters, function(cl) {
    apply(data_mat[pred_char == cl, , drop = FALSE], 2, median)
  }, numeric(ncol(data_mat))))
  rownames(cluster_medians) <- unique_clusters
  median_expanded <- cluster_medians[pred_char, , drop = FALSE]
  mad <- sum(abs(data_mat - median_expanded)) / n

  list(mse = mse, mad = mad)
}

#' @title calculate_modularity
#' @description
#' Calculate modularity score for a graph given an adjacency matrix and cluster
#' assignments
#'
#' @param adj_matrix A sparse adjacency matrix (Matrix package format)
#' @param clusters Vector of cluster assignments for each node
#' @param directed Logical, whether the graph is directed (default: FALSE)
#' @return Modularity value (numeric)
#'
#' @export
calculate_modularity <- function(adj_matrix, clusters, directed = FALSE) {

  # Ensure adjacency matrix is sparse
  if (!inherits(adj_matrix, "sparseMatrix")) {
    adj_matrix <- as(adj_matrix, "dgCMatrix")
  }

  n_nodes <- nrow(adj_matrix)
  if (length(clusters) != n_nodes) {
    stop("Number of cluster assignments must equal number of nodes")
  }

  # Calculate total edges
  if (directed) {
    m <- sum(adj_matrix)
  } else {
    m <- sum(adj_matrix) / 2
  }

  if (m == 0) return(0)

  # Create indicator matrix for clusters
  unique_clusters <- unique(clusters)
  n_clusters <- length(unique_clusters)

  # Create cluster indicator matrix (sparse)
  B <- Matrix::sparseMatrix(
    i = 1:n_nodes,
    j = match(clusters, unique_clusters),
    x = 1,
    dims = c(n_nodes, n_clusters)
  )

  # Calculate modularity using matrix operations (crossprod avoids dense t(B))
  AB <- adj_matrix %*% B
  if (directed) {
    k_in <- Matrix::rowSums(adj_matrix)
    k_out <- Matrix::colSums(adj_matrix)

    # Edges within communities
    trace_term <- sum(Matrix::diag(Matrix::crossprod(B, AB)))

    # Expected edges under null model
    expected_term <- sum(Matrix::crossprod(B, k_out) *
                         Matrix::crossprod(B, k_in))

  } else {
    k <- Matrix::rowSums(adj_matrix)

    # Edges within communities
    trace_term <- sum(Matrix::diag(Matrix::crossprod(B, AB))) / 2

    # Expected edges under null model
    k_communities <- as.vector(Matrix::crossprod(B, k))
    expected_term <- sum(k_communities^2) / (4 * m)
  }

  Q <- trace_term / m - expected_term / m
  return(Q)
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
    dplyr::group_by(.data$resolution) |>
    dplyr::summarize(
      median_score = stats::median(.data$avg_width, na.rm = TRUE),
      variance_score = stats::var(.data$avg_width, na.rm = TRUE),
      standard_error_score = stats::sd(.data$avg_width,
        na.rm = TRUE
      ) / sqrt(length(.data$avg_width)),
      cluster_median_score = stats::median(.data$cluster_median_widths, na.rm = TRUE),
      cluster_variance_score = stats::var(.data$cluster_median_widths, na.rm = TRUE),
      cluster_standard_error_score = stats::sd(.data$cluster_median_widths,
        na.rm = TRUE
      ) /
        sqrt(length(.data$cluster_median_widths))
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
    return(1)
  }

  # This is the ARI equation and will output our value
  ari_result <- (sum_comb_ij - expected_index) / (max_index - expected_index)
  return(ari_result)
}
