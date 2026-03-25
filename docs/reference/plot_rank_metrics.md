# plot_rank_metrics

Create a publication-quality plot showing individual metric ranks and
mean rank across clustering resolutions.

## Usage

``` r
plot_rank_metrics(
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

  Logical; if TRUE (default), highlight the resolution with the lowest
  mean rank with a vertical dashed line.

- base_size:

  Numeric; base font size for theme. Default is 12.

## Value

A ggplot object

## Details

The plot shows:

- Individual metric ranks as colored lines/points

- Mean rank as a bold black line

- Lower ranks indicate better performance (y-axis is reversed)

- The best resolution (lowest mean rank) is optionally highlighted

Uses a colorblind-friendly palette (Okabe-Ito) for metric colors.

## See also

[`suggest_resolution`](https://gladstone-institutes.github.io/clustOpt/reference/suggest_resolution.md),
[`plot_mean_rank`](https://gladstone-institutes.github.io/clustOpt/reference/plot_mean_rank.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cv_results <- clust_opt(seurat_obj, subject_ids = "donor_id")
rankings <- suggest_resolution(cv_results)
plot_rank_metrics(rankings)
plot_rank_metrics(rankings, method = "curvature")
} # }
```
