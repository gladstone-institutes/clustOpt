library(testthat)
library(Seurat)
library(clustOpt)

test_that("split_pca_dimensions works correctly", {
  # Create a Seurat object with random data (with non-uniform variance)
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_2_celltypes_AIDA.rds",
    package = "clustOpt"
  ))
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)


  # Perform PCA on the Seurat object
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

  # Apply the split_pca_dimensions function
  seurat_obj <- split_pca_dimensions(seurat_obj)

  # Check that the new reductions are added
  expect_true("even_pca" %in% names(seurat_obj@reductions))
  expect_true("odd_pca" %in% names(seurat_obj@reductions))

  # Get the dimensions of the original PCA
  pca_dims <- ncol(seurat_obj@reductions$pca@cell.embeddings)

  # Get the dimensions of the even and odd PCAs
  even_pca_dims <- ncol(seurat_obj@reductions$even_pca@cell.embeddings)
  odd_pca_dims <- ncol(seurat_obj@reductions$odd_pca@cell.embeddings)

  # Check that the dimensions are correctly split
  if (pca_dims %% 2 == 0) {
    expect_equal(even_pca_dims, odd_pca_dims)
  } else {
    expect_equal(even_pca_dims, odd_pca_dims + 1)
  }

  # Check that the keys are set correctly
  expect_equal(seurat_obj@reductions$even_pca@key, "even_pca_")
  expect_equal(seurat_obj@reductions$odd_pca@key, "odd_pca_")
})

test_that("Projection respects chosen PCs and avoids data leakage", {
  # Mock Seurat objects to test the function
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_2_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  train_seurat <- prep_train(seurat_obj, subject_ids = "donor_id", test_id = "SG_HEL_H099") |>
    RunPCA(verbose = FALSE) |>
    split_pca_dimensions()
  test_seurat <- prep_test(seurat_obj, subject_ids = "donor_id", test_id = "SG_HEL_H099")

  # Perform the projection
  result <- project_pca(train_seurat, test_seurat, train_with = "even")

  # Verify that the PCs used for clustering (odd PCs in this case) are NOT used for training and testing
  # This assumes that even PCs are used for clustering, so odd PCs should be used for training/prediction.

  train_pcs_used <- colnames(result$train_projected_data)
  test_pcs_used <- colnames(result$test_projected_data)

  # You might customize these checks depending on how the projection is implemented
  expect_true(all(!grepl("PC[2468]", train_pcs_used))) # Check no even PCs are used
  expect_true(all(!grepl("PC[2468]", test_pcs_used))) # Same for test data
})
