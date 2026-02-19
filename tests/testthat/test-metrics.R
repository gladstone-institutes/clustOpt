library(testthat)
library(clustOpt)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# KL Divergence Tests
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("calculate_kl_divergence returns 0 for identical distributions", {
  p <- c(0.25, 0.25, 0.25, 0.25)
  expect_equal(calculate_kl_divergence(p, p), 0)

  # Different distribution
  q <- c(0.5, 0.3, 0.15, 0.05)
  expect_equal(calculate_kl_divergence(q, q), 0)
})

test_that("calculate_kl_divergence returns positive value for different distributions", {
  q <- c(0.5, 0.3, 0.1, 0.1)
  p <- c(0.25, 0.25, 0.25, 0.25)

  result <- calculate_kl_divergence(q, p)
  expect_gt(result, 0)
})

test_that("calculate_kl_divergence is asymmetric", {
  q <- c(0.7, 0.2, 0.1)
  p <- c(0.3, 0.4, 0.3)

  kl_qp <- calculate_kl_divergence(q, p)
  kl_pq <- calculate_kl_divergence(p, q)

  # KL divergence is asymmetric

  expect_false(isTRUE(all.equal(kl_qp, kl_pq)))
})

test_that("calculate_kl_divergence handles zeros in q correctly", {
  # When q[i] = 0, that term contributes 0 to the sum (0 * log(0/p) = 0)
  q <- c(0.5, 0.5, 0, 0)
  p <- c(0.25, 0.25, 0.25, 0.25)

  result <- calculate_kl_divergence(q, p)
  expect_gt(result, 0)
})

test_that("calculate_kl_divergence errors when p is 0 where q is positive", {
  q <- c(0.5, 0.5, 0, 0)
  p <- c(0.5, 0, 0.25, 0.25)  # p[2] = 0 but q[2] = 0.5

  expect_error(
    calculate_kl_divergence(q, p),
    "p must be positive wherever q is positive"
  )
})

test_that("calculate_kl_divergence validates input lengths", {
  q <- c(0.5, 0.5)
  p <- c(0.25, 0.25, 0.25, 0.25)

  expect_error(
    calculate_kl_divergence(q, p),
    "same length"
  )
})

test_that("calculate_kl_divergence validates distributions sum to 1", {
  q <- c(0.5, 0.3)  # sums to 0.8
  p <- c(0.5, 0.5)

  expect_error(
    calculate_kl_divergence(q, p),
    "must sum to 1"
  )
})

test_that("calculate_kl_divergence validates non-negative values", {
  q <- c(0.5, 0.5)
  p <- c(0.7, -0.2)  # negative value

  expect_error(
    calculate_kl_divergence(q, p),
    "non-negative"
  )
})

test_that("calculate_kl_divergence validates numeric input", {
  expect_error(
    calculate_kl_divergence("a", c(0.5, 0.5)),
    "numeric"
  )
})

test_that("calculate_kl_divergence computes known analytic value", {
  # KL divergence from Bernoulli(0.5) to Bernoulli(0.25)
  # D_KL = 0.5 * log(0.5/0.25) + 0.5 * log(0.5/0.75)
  q <- c(0.5, 0.5)
  p <- c(0.25, 0.75)

  expected <- 0.5 * log(0.5 / 0.25) + 0.5 * log(0.5 / 0.75)
  expect_equal(calculate_kl_divergence(q, p), expected, tolerance = 1e-10)
})

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hellinger Distance Tests
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("calculate_hellinger_distance returns 0 for identical distributions", {
  p <- c(0.25, 0.25, 0.25, 0.25)
  expect_equal(calculate_hellinger_distance(p, p), 0)

  q <- c(0.5, 0.3, 0.15, 0.05)
  expect_equal(calculate_hellinger_distance(q, q), 0)
})

test_that("calculate_hellinger_distance returns 1 for disjoint distributions", {
  q <- c(1, 0, 0, 0)
  p <- c(0, 1, 0, 0)

  expect_equal(calculate_hellinger_distance(q, p), 1)
})

test_that("calculate_hellinger_distance is symmetric", {
  q <- c(0.7, 0.2, 0.1)
  p <- c(0.3, 0.4, 0.3)

  h_qp <- calculate_hellinger_distance(q, p)
  h_pq <- calculate_hellinger_distance(p, q)

  expect_equal(h_qp, h_pq)
})

test_that("calculate_hellinger_distance is bounded between 0 and 1", {
  q <- c(0.5, 0.3, 0.1, 0.1)
  p <- c(0.1, 0.2, 0.3, 0.4)

  result <- calculate_hellinger_distance(q, p)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("calculate_hellinger_distance validates input lengths", {
  q <- c(0.5, 0.5)
  p <- c(0.25, 0.25, 0.25, 0.25)

  expect_error(
    calculate_hellinger_distance(q, p),
    "same length"
  )
})

test_that("calculate_hellinger_distance validates distributions sum to 1", {
  q <- c(0.5, 0.3)  # sums to 0.8
  p <- c(0.5, 0.5)

  expect_error(
    calculate_hellinger_distance(q, p),
    "must sum to 1"
  )
})

test_that("calculate_hellinger_distance validates non-negative values", {
  q <- c(0.5, 0.5)
  p <- c(0.7, -0.2)

  expect_error(
    calculate_hellinger_distance(q, p),
    "non-negative"
  )
})

test_that("calculate_hellinger_distance validates numeric input", {
  expect_error(
    calculate_hellinger_distance("a", c(0.5, 0.5)),
    "numeric"
  )
})

test_that("calculate_hellinger_distance computes known analytic value", {
  # Hellinger distance between two 2-element distributions
  # H = sqrt(0.5 * ((sqrt(0.3) - sqrt(0.7))^2 + (sqrt(0.7) - sqrt(0.3))^2))
  q <- c(0.3, 0.7)
  p <- c(0.7, 0.3)

  expected <- sqrt(0.5 * ((sqrt(0.3) - sqrt(0.7))^2 + (sqrt(0.7) - sqrt(0.3))^2))
  expect_equal(calculate_hellinger_distance(q, p), expected, tolerance = 1e-10)
})

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Silhouette Score Tests
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("calculate_silhouette_score returns NA for single cluster", {
  predicted <- rep(1, 10)
  data_frame <- data.frame(x = rnorm(10), y = rnorm(10))

  result <- calculate_silhouette_score(predicted, data_frame)
  expect_true(is.na(result$avg_width))
  expect_true(is.na(result$group_median_width))
})

test_that("calculate_silhouette_score returns valid values for well-separated clusters", {
  set.seed(42)
  # Create two well-separated clusters
  cluster1 <- data.frame(x = rnorm(20, mean = 0), y = rnorm(20, mean = 0))
  cluster2 <- data.frame(x = rnorm(20, mean = 10), y = rnorm(20, mean = 10))
  data_frame <- rbind(cluster1, cluster2)
  predicted <- c(rep(1, 20), rep(2, 20))

  result <- calculate_silhouette_score(predicted, data_frame)

  expect_false(is.na(result$avg_width))
  expect_false(is.na(result$group_median_width))
  # Well-separated clusters should have high silhouette score
  expect_gt(result$avg_width, 0.5)
})

test_that("calculate_silhouette_score returns list with expected elements", {
  set.seed(42)
  data_frame <- data.frame(x = rnorm(30), y = rnorm(30))
  predicted <- c(rep(1, 10), rep(2, 10), rep(3, 10))

  result <- calculate_silhouette_score(predicted, data_frame)

  expect_type(result, "list")
  expect_true("avg_width" %in% names(result))
  expect_true("group_median_width" %in% names(result))
})

test_that("calculate_silhouette_score with precomputed dist gives identical results", {
  set.seed(42)
  cluster1 <- data.frame(x = rnorm(20, mean = 0), y = rnorm(20, mean = 0))
  cluster2 <- data.frame(x = rnorm(20, mean = 10), y = rnorm(20, mean = 10))
  data_frame <- rbind(cluster1, cluster2)
  predicted <- c(rep(1, 20), rep(2, 20))

  result_default <- calculate_silhouette_score(predicted, data_frame)
  precomputed <- dist(data_frame)
  result_precomputed <- calculate_silhouette_score(
    predicted, data_frame, precomputed_dist = precomputed
  )

  expect_equal(result_default$avg_width, result_precomputed$avg_width)
  expect_equal(result_default$group_median_width,
               result_precomputed$group_median_width)
})

test_that("calculate_silhouette_score NULL default still works", {
  set.seed(42)
  data_frame <- data.frame(x = rnorm(30), y = rnorm(30))
  predicted <- c(rep(1, 10), rep(2, 10), rep(3, 10))

  result <- calculate_silhouette_score(predicted, data_frame,
                                       precomputed_dist = NULL)
  expect_type(result, "list")
  expect_false(is.na(result$avg_width))
})

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MSE Score Tests
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("calculate_mse_score returns 0 for single-point clusters", {
  # Each point is its own cluster - distance to centroid is 0
  data_frame <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
  predicted <- c(1, 2, 3)

  result <- calculate_mse_score(predicted, data_frame)

  expect_equal(result$mse, 0)
  expect_equal(result$mad, 0)
})

test_that("calculate_mse_score returns list with mse and mad", {
  set.seed(42)
  data_frame <- data.frame(x = rnorm(30), y = rnorm(30))
  predicted <- c(rep(1, 10), rep(2, 10), rep(3, 10))

  result <- calculate_mse_score(predicted, data_frame)

  expect_type(result, "list")
  expect_true("mse" %in% names(result))
  expect_true("mad" %in% names(result))
  expect_gte(result$mse, 0)
  expect_gte(result$mad, 0)
})

test_that("calculate_mse_score increases with more spread clusters", {
  # Tight cluster
  data_tight <- data.frame(x = rnorm(20, sd = 0.1), y = rnorm(20, sd = 0.1))
  # Spread cluster
  data_spread <- data.frame(x = rnorm(20, sd = 5), y = rnorm(20, sd = 5))

  predicted <- rep(1, 20)

  mse_tight <- calculate_mse_score(predicted, data_tight)$mse
  mse_spread <- calculate_mse_score(predicted, data_spread)$mse

  expect_lt(mse_tight, mse_spread)
})

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Modularity Tests
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("calculate_modularity returns 0 for empty graph", {
  adj_matrix <- matrix(0, nrow = 5, ncol = 5)
  clusters <- c(1, 1, 2, 2, 2)

  result <- calculate_modularity(adj_matrix, clusters)
  expect_equal(result, 0)
})

test_that("calculate_modularity validates cluster length matches nodes", {
  adj_matrix <- matrix(0, nrow = 5, ncol = 5)
  clusters <- c(1, 1, 2)  # Wrong length

  expect_error(
    calculate_modularity(adj_matrix, clusters),
    "Number of cluster assignments must equal number of nodes"
  )
})

test_that("calculate_modularity returns higher value for well-clustered graph", {
  # Create a graph with clear community structure
  # Two communities of 3 nodes each, densely connected within, sparse between
  adj_matrix <- matrix(0, nrow = 6, ncol = 6)
  # Community 1: nodes 1-3
  adj_matrix[1, 2] <- adj_matrix[2, 1] <- 1
  adj_matrix[1, 3] <- adj_matrix[3, 1] <- 1
  adj_matrix[2, 3] <- adj_matrix[3, 2] <- 1
  # Community 2: nodes 4-6
  adj_matrix[4, 5] <- adj_matrix[5, 4] <- 1
  adj_matrix[4, 6] <- adj_matrix[6, 4] <- 1
  adj_matrix[5, 6] <- adj_matrix[6, 5] <- 1

  # Correct clustering
  correct_clusters <- c(1, 1, 1, 2, 2, 2)
  # Wrong clustering
  wrong_clusters <- c(1, 2, 1, 2, 1, 2)

  mod_correct <- calculate_modularity(adj_matrix, correct_clusters)
  mod_wrong <- calculate_modularity(adj_matrix, wrong_clusters)

  expect_gt(mod_correct, mod_wrong)
})
