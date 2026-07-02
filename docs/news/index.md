# Changelog

## clustOpt 1.2.3

### Internal

- [`train_random_forest()`](https://gladstone-institutes.github.io/clustOpt/reference/train_random_forest.md)
  now fits `ranger` through the x/y interface instead of the formula
  interface, avoiding the per-call `model.frame` construction and the
  copy of the projected training matrix, and computes `table(predicted)`
  and `as.character(predicted)` once each instead of twice. Output is
  verified bit-identical; the avoided per-call allocation grows with the
  number of training cells and PCs, so it is neutral on small inputs and
  helps on large ones.
- Verbose per-step timings are now printed in a human-readable format:
  sub-second durations show as milliseconds and durations of a minute or
  more are split into minutes and hours (e.g. `142ms`, `6.3s`, `2m 0s`,
  `1h 3m 7s`) instead of always reporting raw seconds.

## clustOpt 1.2.2

### Performance improvements

- [`prep_test()`](https://gladstone-institutes.github.io/clustOpt/reference/prep_test.md)
  now accepts a `residual_features` argument that restricts the
  test-side
  [`SCTransform()`](https://satijalab.org/seurat/reference/SCTransform.html)
  residual computation to a supplied set of features.
  [`clust_opt()`](https://gladstone-institutes.github.io/clustOpt/reference/clust_opt.md)
  passes the training SCT variable features, which are the only features
  [`project_pca()`](https://gladstone-institutes.github.io/clustOpt/reference/project_PCA.md)
  ever uses, so per-fold residuals are no longer materialized for the
  full transcriptome. The projected output is unchanged (verified
  bit-identical); the redundant work avoided grows with the number of
  cells per held-out subject, so the saving is negligible on small
  subjects and matters on large ones.
- [`clust_opt()`](https://gladstone-institutes.github.io/clustOpt/reference/clust_opt.md)
  no longer serializes the dense O(n^2) silhouette distance matrix to
  every parallel worker. Under a socket-cluster plan (`multisession` or
  remote) the per-subject [`dist()`](https://rdrr.io/r/stats/dist.html)
  is skipped and each worker recomputes it in parallel from the
  already-shipped low-dimensional coordinates, cutting peak memory and
  inter-process transfer. Sequential and forked (`multicore`) plans keep
  the single shared precompute, so there is no regression on the default
  path. The per-resolution `future_lapply()` closure now also captures
  only the scalar hold-out subject id and resolution vector instead of
  the full run grid.

## clustOpt 1.2.1

### Bug fixes

- [`prep_train()`](https://gladstone-institutes.github.io/clustOpt/reference/prep_train.md)
  and
  [`prep_test()`](https://gladstone-institutes.github.io/clustOpt/reference/prep_test.md)
  now drop any pre-existing SCT assay (reset to `RNA` and
  [`DietSeurat()`](https://satijalab.org/seurat/reference/DietSeurat.html))
  before calling
  [`SCTransform()`](https://satijalab.org/seurat/reference/SCTransform.html),
  avoiding warnings about mismatched cells/features when an SCT assay
  built on a different cell subset would otherwise be overwritten.

### Dependencies

- Raised major dependency floors to current releases:
  `ggplot2 (>= 4.0.0)`, `purrr (>= 1.0.0)`, `Seurat (>= 5.4.0)`,
  `dplyr (>= 1.2.0)`, and added explicit floors for `cli (>= 3.4.0)` and
  `SeuratObject (>= 5.3.0)`. The previous `ggplot2 (>= 3.3.5)` floor was
  too low: the plotting code uses
  [`scale_linewidth_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
  and the `linewidth` aesthetic, which require ggplot2 3.4.0+.

### Internal

- Migrated remaining base
  [`warning()`](https://rdrr.io/r/base/warning.html) calls to
  [`cli::cli_warn()`](https://cli.r-lib.org/reference/cli_abort.html)
  for consistent cli-based messaging.
- Replaced superseded
  [`purrr::map_df()`](https://purrr.tidyverse.org/reference/map_dfr.html)
  with `purrr::list_rbind(purrr::map(...))`.

## clustOpt 1.2

### New features

- Added KL divergence
  ([`calculate_kl_divergence()`](https://gladstone-institutes.github.io/clustOpt/reference/calculate_kl_divergence.md))
  and Hellinger distance
  ([`calculate_hellinger_distance()`](https://gladstone-institutes.github.io/clustOpt/reference/calculate_hellinger_distance.md))
  metrics for evaluating cluster distribution consistency between
  training and held-out subjects.
- Added modularity score
  ([`calculate_modularity()`](https://gladstone-institutes.github.io/clustOpt/reference/calculate_modularity.md))
  computed on the precomputed SNN graph.
- Added MSE and MAD scores
  ([`calculate_mse_score()`](https://gladstone-institutes.github.io/clustOpt/reference/calculate_mse_score.md))
  for centroid-based cluster quality evaluation.
- New
  [`suggest_resolution()`](https://gladstone-institutes.github.io/clustOpt/reference/suggest_resolution.md)
  function that ranks resolutions using two complementary methods:
  direct rank aggregation across four metrics, and curvature-based local
  optima detection via second-order finite differences.
- New
  [`summarize_cv_metrics()`](https://gladstone-institutes.github.io/clustOpt/reference/summarize_cv_metrics.md)
  for per-resolution metric summaries.
- New visualization functions
  [`plot_rank_metrics()`](https://gladstone-institutes.github.io/clustOpt/reference/plot_rank_metrics.md)
  and
  [`plot_mean_rank()`](https://gladstone-institutes.github.io/clustOpt/reference/plot_mean_rank.md).
- Logging improvements to make clustOpt messages distinct from its
  dependencies.

### Performance improvements

- Pre-allocate the results list in the main cross-validation loop
  instead of growing it with [`c()`](https://rdrr.io/r/base/c.html),
  reducing memory allocation overhead for many subjects.
- Use [`crossprod()`](https://rdrr.io/r/base/crossprod.html) instead of
  explicit `t() %*%` in PCA projections, avoiding materialization of
  transposed gene-by-cell matrices.
- Use
  [`Matrix::crossprod()`](https://rdrr.io/pkg/Matrix/man/matmult-methods.html)
  in modularity calculation to stay in sparse matrix space and avoid
  dense transposition of the cluster indicator matrix.
- Replace `do.call(rbind, lapply(...))` with `t(vapply(...))` in MAD
  calculation to avoid intermediate list allocation.
- Precompute SNN graph and distance matrix once per held-out subject
  (shared across resolutions) instead of recomputing per resolution.
- Pre-extract cluster assignments from metadata by resolution to avoid
  repeated column lookups.

### Logging and verbosity

- `verbose` parameter now accepts integer levels (0-3) for fine-grained
  control: 0 = silent, 1 = key milestones, 2 = detailed progress, 3 =
  Seurat output.
- Backward-compatible: `verbose = TRUE` maps to level 1, `FALSE` to 0.
- Added per-step timing via `[step_name] Xs` log messages at
  `verbose >= 1`.

### Bug fixes

- Fixed flipped sign in KL divergence calculation.
- Fixed ARI estimate for comparing singleton clusters.
- Handled edge cases in sample validation and metric computation.

### Package reorganization

- Split monolithic `clustOpt.R` and `utils.R` into focused modules:
  `clust_opt.R`, `data_preparation.R`, `metrics.R`, `sketching.R`,
  `validation.R`, `visualization.R`.
- Removed vignette build cache from version control.

## clustOpt 1.1

- Added additional clustering metrics.

## clustOpt 1.0

- Initial public release with subject-wise cross-validation, leverage
  score sketching, BPCells on-disk support, and silhouette-based
  evaluation.
