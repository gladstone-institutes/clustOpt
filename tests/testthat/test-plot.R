library(testthat)
library(clustOpt)

# Example test file for create_sil_plots function
test_that("create_sil_plots returns a list of ggplot objects", {
  # Create a mock sil_dist object (ensure this structure matches what
  # create_sil_plots expects)
  sil_dist <- data.frame(
    resolution = c(1, 1, 2, 2),
    avg_width = runif(4),
    cluster_median_widths = runif(4)
  )

  # Call the function
  result <- create_sil_plots(sil_dist)

  # Check if the result is a list
  expect_type(result, "list")

  # Check if the list contains four elements
  expect_length(result, 4)

  # Check if each element in the list is a ggplot object
  expect_s3_class(result[[1]], "ggplot")
  expect_s3_class(result[[2]], "ggplot")
  expect_s3_class(result[[3]], "ggplot")
  expect_s3_class(result[[4]], "ggplot")
})

# Test for error handling (e.g., passing an invalid sil_dist structure)
test_that("create_sil_plots handles incorrect input gracefully", {
  # Provide an invalid sil_dist object
  sil_dist_invalid <- data.frame(
    invalid_column = c(1, 2, 3),
    some_other_column = c(4, 5, 6)
  )

  # Expect an error
  expect_error(create_sil_plots(sil_dist_invalid))
})

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tests for second_order_diff (internal function)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_that("second_order_diff computes correct values for known input", {
  # Linear input: second differences should be 0
  x_linear <- c(1, 2, 3, 4, 5)
  result <- clustOpt:::second_order_diff(x_linear)
  expect_equal(result[2:4], c(0, 0, 0))

  # Quadratic input: x = c(1, 4, 9, 16, 25) -> second diff should be constant
  x_quad <- c(1, 4, 9, 16, 25)
  result_quad <- clustOpt:::second_order_diff(x_quad)
  # x[i] - (x[i+1] + x[i-1]) / 2 for quadratic = constant
  # i=2: 4 - (9 + 1)/2 = 4 - 5 = -1
  # i=3: 9 - (16 + 4)/2 = 9 - 10 = -1
  # i=4: 16 - (25 + 9)/2 = 16 - 17 = -1
  expect_equal(result_quad[2:4], c(-1, -1, -1))
})

test_that("second_order_diff returns NA at endpoints", {
  x <- c(1, 2, 3, 4, 5)
  result <- clustOpt:::second_order_diff(x)
  expect_true(is.na(result[1]))
  expect_true(is.na(result[5]))
})

test_that("second_order_diff handles length-2 and length-1 vectors", {
  result_2 <- clustOpt:::second_order_diff(c(1, 2))
  expect_equal(length(result_2), 2)
  expect_true(all(is.na(result_2)))

  result_1 <- clustOpt:::second_order_diff(c(1))
  expect_equal(length(result_1), 1)
  expect_true(all(is.na(result_1)))
})

test_that("second_order_diff returns correct output length", {
  for (n in c(3, 5, 10)) {
    x <- seq_len(n)
    result <- clustOpt:::second_order_diff(x)
    expect_equal(length(result), n)
  }
})

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tests for summarize_cv_metrics
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_that("summarize_cv_metrics returns correct structure", {
  cv_results <- data.frame(
    resolution = rep(c(0.1, 0.2, 0.3), each = 3),
    avg_width = c(0.6, 0.65, 0.62, 0.5, 0.52, 0.48, 0.4, 0.42, 0.38),
    KLdivergence = c(0.1, 0.12, 0.11, 0.2, 0.22, 0.18, 0.3, 0.32, 0.28),
    Hellinger = c(0.1, 0.12, 0.11, 0.15, 0.17, 0.14, 0.2, 0.22, 0.19),
    modularity = c(0.7, 0.72, 0.68, 0.6, 0.62, 0.58, 0.5, 0.52, 0.48)
  )

  result <- summarize_cv_metrics(cv_results)

  # One row per resolution
  expect_equal(nrow(result), 3)

  # Expected columns
  expected_cols <- c(
    "resolution", "median_score", "median_KLD_score", "median_Hell_score",
    "median_modularity_score", "standard_error_Hell_score"
  )
  expect_true(all(expected_cols %in% names(result)))
})

test_that("summarize_cv_metrics computes correct median values", {
  cv_results <- data.frame(
    resolution = rep(0.1, 3),
    avg_width = c(0.6, 0.65, 0.62),
    KLdivergence = c(0.1, 0.12, 0.11),
    Hellinger = c(0.1, 0.12, 0.11),
    modularity = c(0.7, 0.72, 0.68)
  )

  result <- summarize_cv_metrics(cv_results)
  expect_equal(result$median_score, median(c(0.6, 0.65, 0.62)))
  expect_equal(result$median_KLD_score, median(c(0.1, 0.12, 0.11)))
  expect_equal(result$median_Hell_score, median(c(0.1, 0.12, 0.11)))
  expect_equal(result$median_modularity_score, median(c(0.7, 0.72, 0.68)))
})

test_that("summarize_cv_metrics computes correct standard error", {
  cv_results <- data.frame(
    resolution = rep(0.1, 5),
    avg_width = runif(5),
    KLdivergence = runif(5),
    Hellinger = c(0.1, 0.2, 0.15, 0.12, 0.18),
    modularity = runif(5)
  )

  result <- summarize_cv_metrics(cv_results)
  expected_se <- sd(c(0.1, 0.2, 0.15, 0.12, 0.18)) / sqrt(5)
  expect_equal(result$standard_error_Hell_score, expected_se)
})

test_that("summarize_cv_metrics handles NA values via na.rm", {
  cv_results <- data.frame(
    resolution = rep(0.1, 3),
    avg_width = c(0.6, NA, 0.62),
    KLdivergence = c(0.1, 0.12, NA),
    Hellinger = c(0.1, 0.12, 0.11),
    modularity = c(NA, 0.72, 0.68)
  )

  result <- summarize_cv_metrics(cv_results)
  expect_equal(nrow(result), 1)
  expect_equal(result$median_score, median(c(0.6, NA, 0.62), na.rm = TRUE))
  expect_equal(result$median_KLD_score, median(c(0.1, 0.12, NA), na.rm = TRUE))
})

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tests for suggest_resolution
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Helper: create cv_results with n resolutions
make_cv_data <- function(n_res = 5, n_folds = 5) {
  resolutions <- seq(0.1, by = 0.1, length.out = n_res)
  data.frame(
    resolution = rep(resolutions, each = n_folds),
    avg_width = rep(seq(0.7, 0.3, length.out = n_res), each = n_folds) +
      rnorm(n_res * n_folds, 0, 0.02),
    KLdivergence = rep(seq(0.05, 0.25, length.out = n_res), each = n_folds) +
      rnorm(n_res * n_folds, 0, 0.01),
    Hellinger = rep(seq(0.05, 0.15, length.out = n_res), each = n_folds) +
      rnorm(n_res * n_folds, 0, 0.01),
    modularity = rep(seq(0.8, 0.4, length.out = n_res), each = n_folds) +
      rnorm(n_res * n_folds, 0, 0.02)
  )
}

test_that("suggest_resolution returns a data frame with expected columns", {
  set.seed(42)
  cv_results <- make_cv_data(5, 5)
  result <- suggest_resolution(cv_results)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)

  # Rank-based columns
  rank_cols <- c("rank_sil", "rank_kl", "rank_hellinger",
                 "rank_modularity", "mean_rank")
  expect_true(all(rank_cols %in% names(result)))

  # Curvature-based columns
  curv_cols <- c("curvature_rank_sil", "curvature_rank_kl",
                 "curvature_rank_hellinger", "curvature_rank_modularity",
                 "curvature_mean_rank")
  expect_true(all(curv_cols %in% names(result)))

  # Upper Hellinger CI

  expect_true("upper_Hell_95ci" %in% names(result))
})

test_that("suggest_resolution is sorted by mean_rank", {
  set.seed(42)
  cv_results <- make_cv_data(5, 5)
  result <- suggest_resolution(cv_results)

  expect_true(all(diff(result$mean_rank) >= 0))
})

test_that("suggest_resolution rank-based columns have no NAs", {
  set.seed(42)
  cv_results <- make_cv_data(5, 5)
  result <- suggest_resolution(cv_results)

  rank_cols <- c("rank_sil", "rank_kl", "rank_hellinger",
                 "rank_modularity", "mean_rank")
  for (col in rank_cols) {
    expect_false(any(is.na(result[[col]])), info = paste("NA found in", col))
  }
})

test_that("suggest_resolution curvature endpoints are NA", {
  set.seed(42)
  cv_results <- make_cv_data(5, 5)
  result <- suggest_resolution(cv_results)

  # Re-sort by resolution to find endpoints
  by_res <- result[order(result$resolution), ]
  curv_cols <- c("curvature_rank_sil", "curvature_rank_kl",
                 "curvature_rank_hellinger", "curvature_rank_modularity")
  for (col in curv_cols) {
    expect_true(is.na(by_res[[col]][1]),
                info = paste("First endpoint not NA in", col))
    expect_true(is.na(by_res[[col]][nrow(by_res)]),
                info = paste("Last endpoint not NA in", col))
  }
})

test_that("suggest_resolution has one row per resolution", {
  set.seed(42)
  cv_results <- make_cv_data(7, 3)
  result <- suggest_resolution(cv_results)
  expect_equal(nrow(result), 7)
  expect_equal(length(unique(result$resolution)), 7)
})

test_that("suggest_resolution handles NA values in input", {
  set.seed(42)
  cv_results <- make_cv_data(5, 5)
  cv_results$avg_width[c(3, 7, 12)] <- NA
  cv_results$KLdivergence[c(5, 15)] <- NA

  result <- suggest_resolution(cv_results)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
})

test_that("suggest_resolution works with 2 resolutions", {
  set.seed(42)
  cv_results <- make_cv_data(2, 5)
  result <- suggest_resolution(cv_results)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  # All curvature ranks should be NA with only 2 resolutions
  expect_true(all(is.na(result$curvature_rank_sil)))
})

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tests for plot_rank_metrics
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_that("plot_rank_metrics returns ggplot with rank method (default)", {
  set.seed(42)
  cv_results <- make_cv_data(5, 5)
  rankings <- suggest_resolution(cv_results)

  p <- plot_rank_metrics(rankings)
  expect_s3_class(p, "ggplot")

  p2 <- plot_rank_metrics(rankings, highlight_best = FALSE)
  expect_s3_class(p2, "ggplot")
})

test_that("plot_rank_metrics works with curvature method", {
  set.seed(42)
  cv_results <- make_cv_data(5, 5)
  rankings <- suggest_resolution(cv_results)

  p <- plot_rank_metrics(rankings, method = "curvature")
  expect_s3_class(p, "ggplot")
})

test_that("plot_rank_metrics curvature errors with fewer than 3 resolutions", {
  set.seed(42)
  cv_results <- make_cv_data(2, 5)
  rankings <- suggest_resolution(cv_results)

  expect_error(
    plot_rank_metrics(rankings, method = "curvature"),
    "At least 3 resolutions"
  )
})

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tests for plot_mean_rank
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_that("plot_mean_rank returns ggplot with rank method (default)", {
  set.seed(42)
  cv_results <- make_cv_data(5, 5)
  rankings <- suggest_resolution(cv_results)

  p <- plot_mean_rank(rankings)
  expect_s3_class(p, "ggplot")

  p2 <- plot_mean_rank(rankings, highlight_best = FALSE)
  expect_s3_class(p2, "ggplot")
})

test_that("plot_mean_rank works with curvature method", {
  set.seed(42)
  cv_results <- make_cv_data(5, 5)
  rankings <- suggest_resolution(cv_results)

  p <- plot_mean_rank(rankings, method = "curvature")
  expect_s3_class(p, "ggplot")
})

test_that("plot_mean_rank curvature errors with fewer than 3 resolutions", {
  set.seed(42)
  cv_results <- make_cv_data(2, 5)
  rankings <- suggest_resolution(cv_results)

  expect_error(
    plot_mean_rank(rankings, method = "curvature"),
    "At least 3 resolutions"
  )
})
