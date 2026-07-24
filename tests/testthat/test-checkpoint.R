library(testthat)
library(Seurat)
library(clustOpt)

fixture <- function() {
  readRDS(system.file(
    "extdata", "1000_cell_sketch_10_donors_3_celltypes_AIDA.rds",
    package = "clustOpt"
  ))
}

run_co <- function(dir, ndim = 10) {
  clust_opt(
    input = fixture(),
    subject_ids = "donor_id",
    ndim = ndim,
    res_range = c(0.01, 0.02),
    num_trees = 10,
    checkpoint_dir = dir
  )
}

test_that("checkpoint_dir writes a manifest and one file per subject", {
  dir <- tempfile("ckpt")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  set.seed(1)
  full <- run_co(dir)

  expect_true(file.exists(file.path(dir, "manifest.rds")))
  subj_files <- list.files(dir, pattern = "^subject-")
  expect_equal(length(subj_files), length(unique(full$test_sample)))
  # no stray temp files left behind by the atomic writer
  expect_length(list.files(dir, pattern = "\\.tmp$"), 0)
})

test_that("a resumed run is bit-identical to an uninterrupted one", {
  dir <- tempfile("ckpt")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  set.seed(1)
  full <- run_co(dir)

  # Drop half the subject checkpoints, then resume from a *different* upstream
  # RNG state: the persisted seed must drive the recompute, not the caller's.
  subj_files <- list.files(dir, pattern = "^subject-", full.names = TRUE)
  file.remove(subj_files[seq_len(floor(length(subj_files) / 2))])

  set.seed(999)
  resumed <- run_co(dir)

  expect_identical(resumed, full)
})

test_that("reusing a checkpoint_dir with a different config errors", {
  dir <- tempfile("ckpt")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  set.seed(1)
  run_co(dir, ndim = 10)

  expect_error(
    run_co(dir, ndim = 8),
    "different clust_opt"
  )
})

test_that("the sketched input is persisted and reused on resume", {
  dir <- tempfile("ckpt")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  # Force the sketch branch without depending on real sketch output: pretend the
  # input is large and sketch to identity so the downstream pipeline stays valid.
  local_mocked_bindings(
    check_size = function(input) TRUE,
    leverage_sketch = function(input, ...) input
  )

  set.seed(1)
  full <- run_co(dir)
  expect_true(file.exists(file.path(dir, "sketch.rds")))

  # Clear the per-subject checkpoints but keep the manifest and sketch, then
  # resume with leverage_sketch booby-trapped: it must not be called.
  file.remove(list.files(dir, pattern = "^subject-", full.names = TRUE))
  local_mocked_bindings(
    leverage_sketch = function(input, ...) stop("leverage_sketch called on resume")
  )

  set.seed(999)
  resumed <- run_co(dir)
  expect_identical(resumed, full)
})
