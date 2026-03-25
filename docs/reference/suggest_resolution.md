# suggest_resolution

Rank clustering resolutions using one of two methods: `"rank"` (default)
ranks all resolutions by their median metric values, while
`"local_optima"` first filters to reproducible resolutions (low upper
Hellinger CI), then ranks the filtered set using both direct rank
aggregation and curvature-based local optima detection.

## Usage

``` r
suggest_resolution(
  cv_results,
  method = c("rank", "local_optima"),
  upper_Hell_score_thresh = NULL,
  upper_Hell_score_thresh_relaxed = NULL,
  min_resolutions = 4
)
```

## Arguments

- cv_results:

  Data frame from
  [`clust_opt`](https://gladstone-institutes.github.io/clustOpt/reference/clust_opt.md)
  cross-validation.

- method:

  Character; `"rank"` (default) ranks all resolutions by median metric
  values without filtering, `"local_optima"` applies the Hellinger
  reproducibility filter and ranks the filtered set by both rank
  aggregation and curvature.

- upper_Hell_score_thresh:

  Numeric. Initial threshold on the upper 95% Hellinger CI for filtering
  reproducible resolutions. Default is scaled dynamically based on the
  number of subjects (calibrated at 0.1 for 11 subjects). Only used when
  `method = "local_optima"`.

- upper_Hell_score_thresh_relaxed:

  Numeric. Relaxed threshold used when fewer than `min_resolutions` pass
  the initial filter. Default is scaled dynamically (calibrated at 0.2
  for 11 subjects). Only used when `method = "local_optima"`.

- min_resolutions:

  Integer. Minimum number of resolutions required to pass the initial
  threshold before relaxing (default 4). Capped at the total number of
  resolutions tested. Only used when `method = "local_optima"`.

## Value

A data frame with one row per resolution (filtered for `"local_optima"`,
all for `"rank"`), containing summarized median metrics and
`upper_Hell_95ci`. For `"local_optima"`: both rank-based (`rank_sil`,
`rank_kl`, `rank_hellinger`, `rank_modularity`, `mean_rank`) and
curvature-based (`curvature_rank_sil`, etc., `curvature_mean_rank`)
columns, sorted by `mean_rank`. For `"rank"`: rank-based columns only,
sorted by `mean_rank`.

## Details

For `method = "local_optima"`, resolutions are filtered to those with
`upper_Hell_95ci < upper_Hell_score_thresh`. If fewer than
`min_resolutions` pass, the threshold is relaxed to
`upper_Hell_score_thresh_relaxed`. An error is raised if no resolutions
pass the relaxed threshold. Curvature values (`second_order_diff`) are
computed on the full, unfiltered resolution sequence to preserve
consecutive spacing, then ranked on only the filtered subset. Both rank
aggregation (`mean_rank`) and curvature (`curvature_mean_rank`) columns
are included in the output for comparison.

The base thresholds (0.1 / 0.2) were validated with 11-subject
leave-one-out cross-validation simulations. When not explicitly set,
thresholds are scaled by `sqrt(11 / n_subjects)` to account for wider
confidence intervals with fewer subjects.

For `method = "rank"`, all resolutions are ranked directly by their
median metric values with no filtering applied.

## See also

[`clust_opt`](https://gladstone-institutes.github.io/clustOpt/reference/clust_opt.md),
[`plot_rank_metrics`](https://gladstone-institutes.github.io/clustOpt/reference/plot_rank_metrics.md),
[`plot_mean_rank`](https://gladstone-institutes.github.io/clustOpt/reference/plot_mean_rank.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cv_results <- clust_opt(seurat_obj, subject_ids = "donor_id")
rankings <- suggest_resolution(cv_results)
rankings_rank <- suggest_resolution(cv_results, method = "rank")
} # }
```
