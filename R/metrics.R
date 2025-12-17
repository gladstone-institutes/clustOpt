#' @importFrom stats median var sd dist
#' @importFrom dplyr group_by summarize
#' @importFrom rlang .data
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clustering Quality Metrics Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#' @title calculate_mse_score
#' @description
#' Calculate mean squared error score
#' @param predicted Cluster assignments
#' @param data_frame Data frame containing the data
#' @return A numeric value representing the mse
#'
#' @export
calculate_mse_score <- function(predicted, data_frame) {

  k <- length(unique(predicted))  # number of clusters
  mse <- 0
  mad <- 0
  n <- nrow(data_frame)

  for (c in unique(predicted)) {
    # Subset points in cluster c
    cluster_points <- data_frame[predicted == c, , drop = FALSE]
    # Compute centroid
    centroid <- colMeans(cluster_points)
    # Squared distances to centroid
    dists <- rowSums((cluster_points - matrix(centroid,
                                              nrow = nrow(cluster_points),
                                              ncol = ncol(data_frame),
                                              byrow = TRUE))^2)
    mse <- mse + sum(dists)

    median_centroid <- apply(cluster_points, 2, median)
    mad_dists <- rowSums(abs(cluster_points - matrix(median_centroid,
                                              nrow = nrow(cluster_points),
                                              ncol = ncol(data_frame),
                                              byrow = TRUE)))

    mad <- mad + sum(mad_dists)


  }

  # Mean squared error
  mse <- mse / n
  # Mean absolute deviation from median
  mad <- mad /n
  return(list(mse=mse, mad=mad))

}

#' @title compute_knn_graph
#' @description
#' Compute k-nearest neighbor graph from coordinate matrix
#'
#' @param coords Numeric matrix where rows are nodes and columns are coordinates
#' @param k Number of nearest neighbors
#' @param mutual Logical, whether to create mutual kNN graph (default: FALSE)
#' @param distance_metric Distance metric to use ("euclidean", "manhattan", "cosine")
#' @return Sparse adjacency matrix of the kNN graph
#'
#' @keywords internal
compute_knn_graph <- function(coords, k, mutual = FALSE, distance_metric = "euclidean") {

  if (!is.matrix(coords) && !is.data.frame(coords)) {
    stop("coords must be a matrix or data.frame")
  }

  coords <- as.matrix(coords)
  n_nodes <- nrow(coords)

  if (k >= n_nodes) {
    stop("k must be less than the number of nodes")
  }

  # Calculate distance matrix
  if (distance_metric == "euclidean") {
    dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
  } else if (distance_metric == "manhattan") {
    dist_matrix <- as.matrix(dist(coords, method = "manhattan"))
  } else if (distance_metric == "cosine") {
    # Cosine distance = 1 - cosine similarity
    coords_norm <- coords / sqrt(rowSums(coords^2))
    cosine_sim <- coords_norm %*% t(coords_norm)
    dist_matrix <- 1 - cosine_sim
  } else {
    stop("Unsupported distance metric. Use 'euclidean', 'manhattan', or 'cosine'")
  }

  # Set diagonal to infinity to avoid self-loops
  diag(dist_matrix) <- Inf

  # Find k nearest neighbors for each node
  knn_indices <- t(apply(dist_matrix, 1, function(x) order(x)[1:k]))

  # Create adjacency matrix
  adj_matrix <- Matrix::sparseMatrix(
    i = rep(1:n_nodes, each = k),
    j = as.vector(t(knn_indices)),
    x = 1,
    dims = c(n_nodes, n_nodes)
  )

  if (mutual) {
    # Keep only mutual connections
    adj_matrix <- adj_matrix * t(adj_matrix)
  }

  return(adj_matrix)
}

#' Calculate modularity for a graph given coordinate matrix and cluster assignments
#'
#' @param coords Numeric matrix where rows are nodes and columns are coordinates
#' @param clusters Vector of cluster assignments for each node
#' @param k Number of nearest neighbors for kNN graph construction
#' @param mutual Logical, whether to create mutual kNN graph (default: FALSE)
#' @param distance_metric Distance metric to use ("euclidean", "manhattan", "cosine")
#' @param directed Logical, whether the graph is directed (default: FALSE)
#' @return List containing modularity value and the kNN adjacency matrix
#'
#' @keywords internal
calculate_modularity_from_coords <- function(clusters, coords, k=20,
                                             mutual = FALSE,
                                             distance_metric = "euclidean",
                                             directed = FALSE) {

  # Compute kNN graph
  adj_matrix <- compute_knn_graph(coords, k, mutual, distance_metric)

  # Calculate modularity
  Q <- calculate_modularity(adj_matrix, clusters, directed)

  return(list(modularity = Q, adjacency_matrix = adj_matrix))
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

  # Calculate modularity using matrix operations
  if (directed) {
    k_in <- Matrix::rowSums(adj_matrix)
    k_out <- Matrix::colSums(adj_matrix)

    # Edges within communities
    trace_term <- sum(Matrix::diag(t(B) %*% adj_matrix %*% B))

    # Expected edges under null model
    expected_term <- sum((t(B) %*% k_out) * (t(B) %*% k_in))

  } else {
    k <- Matrix::rowSums(adj_matrix)

    # Edges within communities
    trace_term <- sum(Matrix::diag(t(B) %*% adj_matrix %*% B)) / 2

    # Expected edges under null model
    k_communities <- as.vector(t(B) %*% k)
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
