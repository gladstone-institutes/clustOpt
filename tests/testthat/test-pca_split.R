library(testthat)
library(Seurat)
library(clustOpt)


library(testthat)
library(Seurat)

test_that("split_pca_dimensions works with 'odd_even' method", {
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_2_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  # Ensure PCA is present
  seurat_obj <- NormalizeData(seurat_obj) |>
    ScaleData() |>
    RunPCA(features = VariableFeatures(seurat_obj),
           verbose = FALSE
    )

  # Split PCA dimensions by odd/even
  result_obj <- split_pca_dimensions(seurat_obj,
    split_method = "odd_even",
    verbose = FALSE
  )

  # Check the result is still a Seurat object
  expect_true(inherits(result_obj, "Seurat"))

  # Check new reductions exist
  expect_true("A_pca" %in% names(result_obj@reductions))
  expect_true("B_pca" %in% names(result_obj@reductions))

  # Original dimension count
  orig_dims <- ncol(seurat_obj@reductions$pca)

  # Check dimension counts in new reductions:
  # A_pca should contain the odd dimensions, B_pca the even.
  odd_count <- length(seq(1, orig_dims, by = 2))
  even_count <- length(seq(2, orig_dims, by = 2))

  expect_equal(ncol(result_obj@reductions$A_pca), odd_count)
  expect_equal(ncol(result_obj@reductions$B_pca), even_count)
})

test_that("split_pca_dimensions works with 'var_balanced' method", {
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_2_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  # Ensure PCA is present
  seurat_obj <- NormalizeData(seurat_obj) |>
    ScaleData() |>
    RunPCA(features = VariableFeatures(seurat_obj),
    verbose = FALSE
  )

  # Split PCA dimensions by variance balancing
  result_obj <- split_pca_dimensions(seurat_obj,
    split_method = "var_balanced",
    verbose = FALSE
  )

  # Check the result is still a Seurat object
  expect_true(inherits(result_obj, "Seurat"))

  # Check new reductions exist
  expect_true("A_pca" %in% names(result_obj@reductions))
  expect_true("B_pca" %in% names(result_obj@reductions))

  # Compare the total variance of each split to ensure it is reasonably balanced
  var_A <- sum(result_obj@reductions$A_pca@stdev^2)
  var_B <- sum(result_obj@reductions$B_pca@stdev^2)
  var_total <- sum(seurat_obj@reductions$pca@stdev^2)

  # Variance should be split so that A and B are not drastically different.
  expect_lt(abs(var_A - var_B), 0.25 * var_total)

  # Also check that we have at least 1 PC in each set
  expect_gt(ncol(result_obj@reductions$A_pca), 0)
  expect_gt(ncol(result_obj@reductions$B_pca), 0)
})
