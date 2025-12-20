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

# Tests for suggest_resolution function
test_that("suggest_resolution returns correct structure", {
  cv_results <- data.frame(
    resolution = rep(c(0.1, 0.2, 0.3), each = 3),
    avg_width = c(0.6, 0.65, 0.62, 0.5, 0.52, 0.48, 0.4, 0.42, 0.38),
    KLdivergence = c(0.1, 0.12, 0.11, 0.2, 0.22, 0.18, 0.3, 0.32, 0.28),
    Hellinger = c(0.1, 0.12, 0.11, 0.15, 0.17, 0.14, 0.2, 0.22, 0.19),
    modularity = c(0.7, 0.72, 0.68, 0.6, 0.62, 0.58, 0.5, 0.52, 0.48)
  )

  result <- suggest_resolution(cv_results)

  # Check output is a data frame

  expect_s3_class(result, "data.frame")

  # Check required columns exist
  expected_cols <- c(
    "resolution", "median_sil", "median_kl", "median_hellinger",
    "median_modularity", "rank_sil", "rank_kl", "rank_hellinger",
    "rank_modularity", "mean_rank"
  )
  expect_true(all(expected_cols %in% names(result)))

  # Check result is sorted by mean_rank
  expect_equal(result$mean_rank, sort(result$mean_rank))

  # Check we have one row per resolution

  expect_equal(nrow(result), 3)
})

test_that("suggest_resolution handles NA values", {
  cv_results <- data.frame(
    resolution = rep(c(0.1, 0.2), each = 3),
    avg_width = c(0.6, NA, 0.62, 0.5, 0.52, 0.48),
    KLdivergence = c(0.1, 0.12, 0.11, 0.2, 0.22, 0.18),
    Hellinger = c(0.1, 0.12, 0.11, 0.15, 0.17, 0.14),
    modularity = c(0.7, 0.72, 0.68, 0.6, 0.62, 0.58)
  )

  result <- suggest_resolution(cv_results)

  # Should complete without error and drop NA rows
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})

# Tests for plot_rank_metrics function
test_that("plot_rank_metrics returns ggplot object", {
  rankings <- data.frame(
    resolution = c(0.1, 0.2, 0.3),
    median_sil = c(0.6, 0.5, 0.4),
    median_kl = c(0.1, 0.2, 0.3),
    median_hellinger = c(0.1, 0.15, 0.2),
    median_modularity = c(0.7, 0.6, 0.5),
    rank_sil = c(1, 2, 3),
    rank_kl = c(1, 2, 3),
    rank_hellinger = c(1, 2, 3),
    rank_modularity = c(1, 2, 3),
    mean_rank = c(1, 2, 3)
  )

  p <- plot_rank_metrics(rankings)
  expect_s3_class(p, "ggplot")

  # Test with highlight_best = FALSE
  p2 <- plot_rank_metrics(rankings, highlight_best = FALSE)
  expect_s3_class(p2, "ggplot")
})

# Tests for plot_mean_rank function
test_that("plot_mean_rank returns ggplot object", {
  rankings <- data.frame(
    resolution = c(0.1, 0.2, 0.3),
    mean_rank = c(1, 2, 3)
  )

  p <- plot_mean_rank(rankings)
  expect_s3_class(p, "ggplot")

  # Test with highlight_best = FALSE
  p2 <- plot_mean_rank(rankings, highlight_best = FALSE)
  expect_s3_class(p2, "ggplot")
})
