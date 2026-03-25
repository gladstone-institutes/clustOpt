# plot_mean_rank

Create a publication-quality plot showing the mean rank across
resolutions. Lower mean rank indicates better overall performance across
all metrics.

## Usage

``` r
plot_mean_rank(
  rank_results,
  method = "rank",
  highlight_best = TRUE,
  base_size = 12
)
```

## Arguments

- rank_results:

  A data.frame output from
  [`suggest_resolution`](https://gladstone-institutes.github.io/clustOpt/reference/suggest_resolution.md).

- method:

  Character; `"rank"` (default) for direct rank aggregation, or
  `"curvature"` for curvature-based ranking. The curvature method
  requires at least 3 resolutions.

- highlight_best:

  Logical; if TRUE (default), highlight the best resolution with a
  vertical dashed line.

- base_size:

  Numeric; base font size for theme. Default is 12.

## Value

A ggplot object

## See also

[`suggest_resolution`](https://gladstone-institutes.github.io/clustOpt/reference/suggest_resolution.md),
[`plot_rank_metrics`](https://gladstone-institutes.github.io/clustOpt/reference/plot_rank_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cv_results <- clust_opt(seurat_obj, subject_ids = "donor_id")
rankings <- suggest_resolution(cv_results)
plot_mean_rank(rankings)
plot_mean_rank(rankings, method = "curvature")
} # }
```
