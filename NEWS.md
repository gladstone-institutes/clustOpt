# clustOpt 1.2

## New features

- Added KL divergence (`calculate_kl_divergence()`) and Hellinger distance
  (`calculate_hellinger_distance()`) metrics for evaluating cluster distribution
  consistency between training and held-out subjects.
- Added modularity score (`calculate_modularity()`) computed on the precomputed
  SNN graph.
- Added MSE and MAD scores (`calculate_mse_score()`) for centroid-based cluster
  quality evaluation.
- New `suggest_resolution()` function that ranks resolutions using two
  complementary methods: direct rank aggregation across four metrics, and
  curvature-based local optima detection via second-order finite differences.
- New `summarize_cv_metrics()` for per-resolution metric summaries.
- New visualization functions `plot_rank_metrics()` and `plot_mean_rank()`.
- Logging improvements to make clustOpt messages distinct from its dependencies.

## Performance improvements

- Pre-allocate the results list in the main cross-validation loop instead of
  growing it with `c()`, reducing memory allocation overhead for many subjects.
- Use `crossprod()` instead of explicit `t() %*%` in PCA projections, avoiding
  materialization of transposed gene-by-cell matrices.
- Use `Matrix::crossprod()` in modularity calculation to stay in sparse matrix
  space and avoid dense transposition of the cluster indicator matrix.
- Replace `do.call(rbind, lapply(...))` with `t(vapply(...))` in MAD
  calculation to avoid intermediate list allocation.
- Precompute SNN graph and distance matrix once per held-out subject (shared
  across resolutions) instead of recomputing per resolution.
- Pre-extract cluster assignments from metadata by resolution to avoid repeated
  column lookups.

## Logging and verbosity

- `verbose` parameter now accepts integer levels (0-3) for fine-grained control:
  0 = silent, 1 = key milestones, 2 = detailed progress, 3 = Seurat output.
- Backward-compatible: `verbose = TRUE` maps to level 1, `FALSE` to 0.
- Added per-step timing via `[step_name] Xs` log messages at `verbose >= 1`.

## Bug fixes

- Fixed flipped sign in KL divergence calculation.
- Fixed ARI estimate for comparing singleton clusters.
- Handled edge cases in sample validation and metric computation.

## Package reorganization

- Split monolithic `clustOpt.R` and `utils.R` into focused modules:
  `clust_opt.R`, `data_preparation.R`, `metrics.R`, `sketching.R`,
  `validation.R`, `visualization.R`.
- Removed vignette build cache from version control.

# clustOpt 1.1

- Added additional clustering metrics.


# clustOpt 1.0

- Initial public release with subject-wise cross-validation, leverage score
  sketching, BPCells on-disk support, and silhouette-based evaluation.
