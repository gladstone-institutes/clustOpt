library(testthat)
library(clustOpt)

# --- normalize_verbose() ---

test_that("normalize_verbose converts TRUE to 1L", {
  expect_identical(clustOpt:::normalize_verbose(TRUE), 1L)
})

test_that("normalize_verbose converts FALSE to 0L", {
  expect_identical(clustOpt:::normalize_verbose(FALSE), 0L)
})

test_that("normalize_verbose passes through integers", {
  expect_identical(clustOpt:::normalize_verbose(0), 0L)
  expect_identical(clustOpt:::normalize_verbose(1), 1L)
})

test_that("normalize_verbose coerces numeric to integer", {
  expect_identical(clustOpt:::normalize_verbose(2.7), 2L)
})

test_that("normalize_verbose errors on invalid input", {
  expect_error(clustOpt:::normalize_verbose("high"))
  expect_error(clustOpt:::normalize_verbose(-1))
  expect_error(clustOpt:::normalize_verbose(NA))
  expect_error(clustOpt:::normalize_verbose(c(1, 2)))
})

# --- get_valid_samples() verbose behavior ---

test_that("get_valid_samples shows messages at verbose=1", {
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_3_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  expect_message(
    get_valid_samples(seurat_obj, "donor_id", 50, verbose = 1),
    "Using"
  )
})

test_that("get_valid_samples is silent at verbose=0", {
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_3_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  expect_silent(
    get_valid_samples(seurat_obj, "donor_id", 50, verbose = 0)
  )
})

# --- split_pca_dimensions() verbose behavior ---

test_that("split_pca_dimensions messages at verbose=2, silent at 0 and 1", {
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_3_celltypes_AIDA.rds",
    package = "clustOpt"
  ))
  seurat_obj <- Seurat::SCTransform(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = 10, verbose = FALSE)

  expect_message(
    clustOpt:::split_pca_dimensions(seurat_obj, verbose = 2),
    "Splitting PCs"
  )
  expect_silent(
    clustOpt:::split_pca_dimensions(seurat_obj, verbose = 0)
  )
  expect_silent(
    clustOpt:::split_pca_dimensions(seurat_obj, verbose = 1)
  )
})

# --- Backward compatibility: verbose = TRUE works in clust_opt ---

test_that("clust_opt accepts verbose = TRUE without error", {
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_3_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  expect_no_error(
    clust_opt(
      input = seurat_obj,
      subject_ids = "donor_id",
      ndim = 10,
      res_range = c(0.01, 0.02),
      num_trees = 10,
      verbose = TRUE
    )
  )
})
