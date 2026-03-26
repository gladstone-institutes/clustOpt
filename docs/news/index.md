# Changelog

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
