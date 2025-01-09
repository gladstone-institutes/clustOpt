library(testthat)
library(Seurat)
library(clustOpt)


test_that("Projection respects chosen PCs and avoids data leakage", {
  # Tutorial data set
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_2_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  train_seurat <- prep_train(seurat_obj,
    subject_ids = "donor_id",
    test_id = "JP_RIK_H079"
  ) |>
    RunPCA(verbose = FALSE) |>
    split_pca_dimensions()
  test_seurat <- prep_test(seurat_obj,
    subject_ids = "donor_id",
    test_id = "JP_RIK_H079"
  )

  # Perform the projection
  result <- project_pca(train_seurat, test_seurat,
    train_with = "B",
    verbose = TRUE
  )

  # Verify that the PCs used for clustering (odd PCs in this case) are NOT used
  # for training and testing
  train_pcs_used <- colnames(result$train_projected_data)
  test_pcs_used <- colnames(result$test_projected_data)

  # Check no even PCs are used
  expect_true(all(!grepl("PC[2468]", train_pcs_used)))
  expect_true(all(!grepl("PC[2468]", test_pcs_used))) # Same for test data
})
