#' @include metrics.R
#' @importFrom rlang .data
#' @importFrom tidyr drop_na pivot_longer
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw labs geom_errorbar
#' @importFrom ggplot2 geom_point geom_line scale_color_manual scale_y_reverse
#' @importFrom ggplot2 geom_vline theme_classic element_text theme
NULL


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualization Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title create_sil_plots
#' @description
#' Silhouette score distribution plots
#'
#' @param sil_dist the output of clustOpt
#' @return A list of ggplot objects
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw
#' @importFrom ggplot2 labs geom_errorbar geom_point geom_line
create_sil_plots <- function(sil_dist) {
  sil_summary_data <- sil_summary(sil_dist)

  plot1 <- sil_dist |>
    ggplot2::ggplot(ggplot2::aes(x = as.factor(.data$resolution), y = .data$avg_width)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Resolution",
      y = "Avg. Silhouette Score Across All Cells"
    )

  plot2 <- sil_dist |>
    ggplot2::ggplot(ggplot2::aes(
      x = as.factor(.data$resolution),
      y = .data$cluster_median_widths
    )) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Resolution",
      y = "Median Silhouette Score Across Clusters"
    )

  plot3 <- ggplot2::ggplot(
    sil_summary_data,
    ggplot2::aes(x = as.factor(.data$resolution), y = .data$median_score, group = 1)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$median_score - (1.96 * .data$standard_error_score),
        ymax = .data$median_score + (1.96 * .data$standard_error_score),
        width = .3
      ),
      color = "red"
    ) +
    ggplot2::geom_point(colour = "#619CFF") +
    ggplot2::geom_line(colour = "#619CFF") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Resolution",
      y = "Median Silhouette Score Across All Cells"
    )

  plot4 <- ggplot2::ggplot(
    sil_summary_data,
    ggplot2::aes(
      x = as.factor(.data$resolution),
      y = .data$cluster_median_score,
      group = 1
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$cluster_median_score - (1.96 * .data$standard_error_score),
        ymax = .data$cluster_median_score + (1.96 * .data$standard_error_score),
        width = .3
      ),
      color = "red"
    ) +
    ggplot2::geom_point(colour = "#619CFF") +
    ggplot2::geom_line(colour = "#619CFF") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Resolution",
      y = "Median Silhouette Score Across Clusters"
    )

  list(plot1, plot2, plot3, plot4)
}

#' @title suggest_resolution
#' @description
#' Aggregate multiple clustering quality metrics using rank-based analysis to
#' suggest optimal clustering resolution(s). This function computes median values
#' for each metric across cross-validation folds, ranks resolutions by each metric,
#' and calculates a mean rank for comprehensive comparison.
#'
#' @param cv_results A data.frame output from \code{\link{clust_opt}} containing
#'   cross-validation results with columns: resolution, avg_width, KLdivergence,
#'   Hellinger, and modularity.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{resolution}: The resolution parameter values
#'     \item \code{median_sil}: Median silhouette score (avg_width)
#'     \item \code{median_kl}: Median KL divergence
#'     \item \code{median_hellinger}: Median Hellinger distance
#'     \item \code{median_modularity}: Median modularity score
#'     \item \code{rank_sil}: Rank by silhouette (higher is better, ranked descending)
#'     \item \code{rank_kl}: Rank by KL divergence (lower is better)
#'     \item \code{rank_hellinger}: Rank by Hellinger distance (lower is better)
#'     \item \code{rank_modularity}: Rank by modularity (higher is better, ranked descending)
#'     \item \code{mean_rank}: Average of the individual ranks
#'   }
#'   The tibble is sorted by mean_rank (ascending), with the best resolution first.
#'
#' @details
#' The rank aggregation approach helps identify resolutions that perform well across
#' multiple metrics. When individual metrics conflict, the mean rank provides a
#' balanced recommendation.
#'
#' Metric interpretations:
#' \itemize{
#'   \item Silhouette: Higher is better (well-separated clusters)
#'   \item KL Divergence: Lower is better (reproducible cluster proportions)
#'   \item Hellinger Distance: Lower is better (reproducible cluster proportions)
#'   \item Modularity: Higher is better (well-defined graph clusters)
#' }
#'
#' @examples
#' \dontrun{
#' cv_results <- clust_opt(seurat_obj, subject_ids = "donor_id")
#'
#' # Get resolution recommendations
#' rankings <- suggest_resolution(cv_results)
#'
#' # View top recommendations
#' head(rankings)
#' }
#'
#' @seealso \code{\link{clust_opt}}, \code{\link{plot_rank_metrics}}
#' @export
suggest_resolution <- function(cv_results) {
  cv_results |>
    tidyr::drop_na() |>
    dplyr::group_by(.data$resolution) |>
    dplyr::summarize(
      median_sil = stats::median(.data$avg_width),
      median_kl = stats::median(.data$KLdivergence),
      median_hellinger = stats::median(.data$Hellinger),
      median_modularity = stats::median(.data$modularity)
    ) |>
    dplyr::mutate(
      rank_sil = rank(-.data$median_sil),
      rank_kl = rank(.data$median_kl),
      rank_hellinger = rank(.data$median_hellinger),
      rank_modularity = rank(-.data$median_modularity)
    ) |>
    dplyr::mutate(
      mean_rank = (.data$rank_sil + .data$rank_kl +
                     .data$rank_hellinger + .data$rank_modularity) / 4
    ) |>
    dplyr::arrange(.data$mean_rank)
}

#' @title plot_rank_metrics
#' @description
#' Create a publication-quality plot showing individual metric ranks and mean rank
#' across clustering resolutions. This visualization helps identify optimal resolutions
#' by comparing how each resolution ranks across different quality metrics.
#'
#' @param rank_results A data.frame output from \code{\link{suggest_resolution}}
#'   containing rank columns (rank_sil, rank_kl, rank_hellinger, rank_modularity,
#'   mean_rank).
#' @param highlight_best Logical; if TRUE (default), highlight the resolution with
#'   the lowest mean rank with a vertical dashed line.
#' @param base_size Numeric; base font size for theme. Default is 12.
#'
#' @return A ggplot object
#'
#' @details
#' The plot shows:
#' \itemize{
#'   \item Individual metric ranks as colored lines/points
#'   \item Mean rank as a bold black line
#'   \item Lower ranks indicate better performance (y-axis is reversed)
#'   \item The best resolution (lowest mean rank) is optionally highlighted
#' }
#'
#' Uses a colorblind-friendly palette (Okabe-Ito) for metric colors.
#'
#' @examples
#' \dontrun{
#' rankings <- suggest_resolution(cv_results)
#' plot_rank_metrics(rankings)
#' }
#'
#' @seealso \code{\link{suggest_resolution}}, \code{\link{plot_mean_rank}}
#' @export
plot_rank_metrics <- function(rank_results, highlight_best = TRUE, base_size = 12) {
  # Reshape data for plotting - include mean_rank as a metric
  plot_data <- rank_results |>
    tidyr::pivot_longer(
      cols = c("rank_sil", "rank_kl", "rank_hellinger", "rank_modularity", "mean_rank"),
      names_to = "metric",
      values_to = "rank"
    ) |>
    dplyr::mutate(
      resolution = as.factor(.data$resolution),
      metric = factor(.data$metric,
        levels = c("rank_sil", "rank_kl", "rank_hellinger", "rank_modularity", "mean_rank"),
        labels = c("Silhouette", "KL Divergence", "Hellinger", "Modularity", "Mean Rank")
      )
    )

  # Colorblind-friendly palette (Okabe-Ito) plus black for mean rank
  colors <- c(
    "Silhouette" = "#E69F00",
    "KL Divergence" = "#56B4E9",
    "Hellinger" = "#009E73",
    "Modularity" = "#CC79A7",
    "Mean Rank" = "black"
  )

  # Set line widths - thicker for mean rank
  linewidths <- c(
    "Silhouette" = 0.8,
    "KL Divergence" = 0.8,
    "Hellinger" = 0.8,
    "Modularity" = 0.8,
    "Mean Rank" = 1.2
  )

  # Set point sizes - larger for mean rank
  point_sizes <- c(
    "Silhouette" = 2.5,
    "KL Divergence" = 2.5,
    "Hellinger" = 2.5,
    "Modularity" = 2.5,
    "Mean Rank" = 3.5
  )

  # Set alpha - more opaque for mean rank
  alphas <- c(
    "Silhouette" = 0.7,
    "KL Divergence" = 0.7,
    "Hellinger" = 0.7,
    "Modularity" = 0.7,
    "Mean Rank" = 1.0
  )

  n_resolutions <- nrow(rank_results)

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$resolution,
      y = .data$rank,
      color = .data$metric,
      linewidth = .data$metric,
      size = .data$metric,
      alpha = .data$metric,
      group = .data$metric
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = colors, name = "Metric") +
    ggplot2::scale_linewidth_manual(values = linewidths, guide = "none") +
    ggplot2::scale_size_manual(values = point_sizes, guide = "none") +
    ggplot2::scale_alpha_manual(values = alphas, guide = "none") +
    ggplot2::scale_y_reverse(breaks = seq_len(n_resolutions)) +
    ggplot2::labs(
      x = "Resolution",
      y = "Rank (lower is better)",
      title = "Resolution Quality Ranking by Metric"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (highlight_best) {
    best_res <- as.factor(rank_results$resolution[which.min(rank_results$mean_rank)])
    p <- p + ggplot2::geom_vline(
      xintercept = best_res,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    )
  }

  p
}

#' @title plot_mean_rank
#' @description
#' Create a publication-quality plot showing the mean rank across resolutions.
#' Lower mean rank indicates better overall performance across all metrics.
#'
#' @param rank_results A data.frame output from \code{\link{suggest_resolution}}.
#' @param highlight_best Logical; if TRUE (default), highlight the best resolution
#'   with a larger outlined point.
#' @param base_size Numeric; base font size for theme. Default is 12.
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' rankings <- suggest_resolution(cv_results)
#' plot_mean_rank(rankings)
#' }
#'
#' @seealso \code{\link{suggest_resolution}}, \code{\link{plot_rank_metrics}}
#' @export
plot_mean_rank <- function(rank_results, highlight_best = TRUE, base_size = 12) {
  plot_data <- rank_results |>
    dplyr::mutate(resolution = as.factor(.data$resolution))

  n_resolutions <- nrow(rank_results)

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$resolution, y = .data$mean_rank, group = 1)
  ) +
    ggplot2::geom_line(linewidth = 1, color = "#0072B2") +
    ggplot2::geom_point(size = 4, color = "#0072B2") +
    ggplot2::scale_y_reverse(breaks = seq_len(n_resolutions)) +
    ggplot2::labs(
      x = "Resolution",
      y = "Mean Rank (lower is better)",
      title = "Overall Resolution Quality"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (highlight_best) {
    best_res <- as.factor(rank_results$resolution[which.min(rank_results$mean_rank)])
    p <- p + ggplot2::geom_vline(
      xintercept = best_res,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    )
  }

  p
}
