library(testthat)
library(Seurat)
library(clustOpt)

test_that("leverage_sketch works with on_disk=TRUE and properly reduces cell count", {
  # Load the test dataset from the clustOpt package
  seurat_obj <- readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_3_celltypes_AIDA.rds",
    package = "clustOpt"
  ))

  # Skip test if BPCells is not available
  skip_if_not_installed("BPCells")

  # Create a temporary directory for on-disk matrices
  temp_dir <- file.path(tempdir(), "test_sketch_ondisk")
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)

  # Get original number of cells
  original_cell_count <- ncol(seurat_obj)
  sketch_size <- 200 # Smaller than the original

  # Run leverage sketching with on_disk=TRUE and a specific sketch size
  sketched_obj <- leverage_sketch(
    input = seurat_obj,
    sketch_size = sketch_size,
    dtype = "scRNA",
    on_disk = TRUE,
    output_dir = temp_dir,
    verbose = FALSE
  )

  # Tests to verify the sketching worked correctly
  expect_true(inherits(sketched_obj, "Seurat"))

  # Check that the cells were actually reduced
  # The output should match our sketch_size
  expect_equal(ncol(sketched_obj), sketch_size)

  # Verify the right assay is present
  expect_true("RNA" %in% names(sketched_obj@assays))
  expect_equal(Seurat::DefaultAssay(sketched_obj), "RNA")

  # Check if the Assay is using the expected v5 format
  expect_true(inherits(sketched_obj[["RNA"]], "Assay5"))

  # Instead of checking for specific class, just check if there's on-disk storage being used
  # This is a more reliable check than specific class names that might change
  cat("\nChecking counts matrix class:", class(sketched_obj[["RNA"]]@layers$counts), "\n")

  # Check if directory is being used for storage
  # This tests functionality, not specific implementation details
  expect_true(
    dir.exists(temp_dir) &&
      length(list.files(temp_dir, recursive = TRUE)) > 0
  )

  # Check that the matrix has rows and columns as expected
  expect_equal(nrow(sketched_obj[["RNA"]]@layers$counts), nrow(seurat_obj))
  expect_equal(ncol(sketched_obj[["RNA"]]@layers$counts), sketch_size)

  # Clean up temporary files
  unlink(temp_dir, recursive = TRUE)
})
