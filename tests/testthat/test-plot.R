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
