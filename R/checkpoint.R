#' @importFrom Seurat DefaultAssay
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Checkpoint / resume helpers
#
# All file I/O for checkpointing runs in the main clust_opt() process, never in
# a future worker, so there is no connection serialization or write contention.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Save an object to RDS atomically
#'
#' Writes to a sibling \code{.tmp} file then renames, so a killed job never
#' leaves a half-written checkpoint. Only called for paths that do not yet
#' exist, so the rename never has to clobber a destination.
#'
#' @param obj Object to serialize.
#' @param path Destination path.
#' @return \code{path}, invisibly.
#' @keywords internal
.save_rds_atomic <- function(obj, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Failed to write checkpoint to ", path)
  }
  invisible(path)
}

#' Build the run configuration fingerprint stored in the manifest
#'
#' Captures the run-defining arguments plus a cheap input signature (dims,
#' assay names, and a hash of the subject id column). The expression matrix is
#' deliberately not hashed.
#'
#' @keywords internal
.checkpoint_config <- function(input, ndim, dtype, sketch_size, skip_sketch,
                               subject_ids, res_range, within_batch, num_trees,
                               train_with, min_cells) {
  list(
    ndim = ndim,
    dtype = dtype,
    sketch_size = sketch_size,
    skip_sketch = skip_sketch,
    subject_ids = subject_ids,
    res_range = res_range,
    within_batch = within_batch,
    num_trees = num_trees,
    train_with = train_with,
    min_cells = min_cells,
    input_dim = c(nrow(input), ncol(input)),
    default_assay = Seurat::DefaultAssay(input),
    assay_names = sort(names(input@assays)),
    subject_hash = digest::digest(input@meta.data[[subject_ids]])
  )
}

#' Installed clustOpt version as a string
#'
#' Wrapped in a helper so it can be stubbed in tests.
#'
#' @keywords internal
.clustopt_version <- function() {
  as.character(utils::packageVersion("clustOpt"))
}

#' Initialize or validate the checkpoint manifest
#'
#' On first run, writes a manifest holding the clustOpt version, the config
#' fingerprint, and the resolved seed. On resume, validates that the version and
#' config both match (erroring on drift, version first for a clear diagnostic)
#' and returns the persisted seed so the sketch and per-subject seeding are
#' identical to the original run.
#'
#' @param dir Checkpoint directory.
#' @param config Config list from \code{.checkpoint_config()}.
#' @param seed Resolved integer seed for this run.
#' @param version Installed clustOpt version string.
#' @return The seed to actually use (persisted seed on resume, else \code{seed}).
#' @keywords internal
.checkpoint_init_manifest <- function(dir, config, seed, version) {
  path <- file.path(dir, "manifest.rds")
  if (file.exists(path)) {
    prev <- readRDS(path)
    if (!identical(prev$version, version)) {
      stop(
        "checkpoint_dir '", dir, "' was created with clustOpt ", prev$version,
        " but this session is running clustOpt ", version, ". Resuming across ",
        "versions could mix results from different algorithm versions. Use a ",
        "fresh checkpoint_dir, or reinstall clustOpt ", prev$version,
        " to resume this run.",
        call. = FALSE
      )
    }
    if (!identical(prev$config, config)) {
      stop(
        "checkpoint_dir '", dir, "' was created with a different clust_opt() ",
        "configuration or input. Reusing it could mix incompatible results. ",
        "Point checkpoint_dir at a fresh directory, or delete this one to ",
        "start over.",
        call. = FALSE
      )
    }
    return(prev$seed)
  }
  .save_rds_atomic(list(version = version, config = config, seed = seed), path)
  seed
}

#' Deterministic per-subject seed
#'
#' Derives a stable integer seed from the run seed and the subject name, so a
#' subject reproduces exactly regardless of which subjects ran before it in the
#' process (the basis for bit-identical resume).
#'
#' @keywords internal
.subject_seed <- function(seed, sam) {
  digest::digest2int(paste0(seed, "-", as.character(sam)))
}

#' Canonical checkpoint path for one holdout subject
#'
#' Index is zero-padded to the width of the subject count so files sort in run
#' order; the (sanitized) subject name is embedded only for readability. Both
#' index and name are deterministic given an identical config and persisted
#' sketch, so the exact path is reproducible on resume.
#'
#' @keywords internal
.subject_checkpoint_path <- function(dir, idx, sam, n_total) {
  width <- nchar(as.character(n_total))
  padded <- formatC(idx, width = width, flag = "0")
  safe <- gsub("[^A-Za-z0-9._-]", "_", as.character(sam))
  file.path(dir, sprintf("subject-%s_%s.rds", padded, safe))
}
