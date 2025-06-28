library(testthat)
library(Seurat)
library(clustOpt)

test_that("split_pca_dimensions works", {
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_3_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  # Ensure PCA is present
  seurat_obj <- NormalizeData(seurat_obj) |>
    ScaleData() |>
    RunPCA(
      features = VariableFeatures(seurat_obj),
      verbose = FALSE
    )

  # Split PCA dimensions by odd/even
  result_obj <- split_pca_dimensions(seurat_obj,
    verbose = FALSE
  )

  # Check the result is still a Seurat object
  expect_true(inherits(result_obj, "Seurat"))

  # Check new reductions exist
  expect_true("odd_pca" %in% names(result_obj@reductions))
  expect_true("even_pca" %in% names(result_obj@reductions))

  # Original dimension count
  orig_dims <- ncol(seurat_obj@reductions$pca)

  # Check dimension counts in new reductions:
  # odd_pca should contain the odd dimensions, even_pca the even.
  odd_count <- length(seq(1, orig_dims, by = 2))
  even_count <- length(seq(2, orig_dims, by = 2))

  expect_equal(ncol(result_obj@reductions$odd_pca), odd_count)
  expect_equal(ncol(result_obj@reductions$even_pca), even_count)
})

