library(testthat)
library(Seurat)
library(clustOpt)
library(dplyr)

test_that("clust_opt runs", {
  # Tutorial data set
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_3_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  df <- clust_opt(input = seurat_obj,
                           subject_ids = "donor_id",
                           ndim = 10,
                           res_range = c(0.01,0.02),
                           num_trees = 10)
  
  expect_true(exists("df"))
})
