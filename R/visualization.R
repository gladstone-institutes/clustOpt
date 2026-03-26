#' @include metrics.R
#' @importFrom rlang .data
#' @importFrom dplyr all_of
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
 
#' Compute second-order finite differences for a numeric vector
#'
#' @param x Numeric vector, ordered by resolution.
#' @return Numeric vector of same length, with NA at endpoints.
#' @keywords internal
second_order_diff <- function(x) {
  n <- length(x)
  if (n < 3) return(rep(NA_real_, n))
  c(NA_real_,
    x[2:(n - 1)] - (x[3:n] + x[1:(n - 2)]) / 2,
    NA_real_)
}

#' @title summarize_cv_metrics
#' @description Summarize cross-validation metrics across resolutions.
#' @param cv_results Data frame from clust_opt cross-validation.
#' @return Summarized data frame with one row per resolution.
#' @export
summarize_cv_metrics <- function(cv_results) {
  cv_results |>
    dplyr::group_by(.data$resolution) |>
    dplyr::summarize(
      median_score            = median(.data$avg_width, na.rm = TRUE),
      median_KLD_score        = median(.data$KLdivergence, na.rm = TRUE),
      median_Hell_score       = median(.data$Hellinger, na.rm = TRUE),
      median_modularity_score = median(.data$modularity, na.rm = TRUE),
      standard_error_Hell_score = sd(.data$Hellinger, na.rm = TRUE) /
        sqrt(length(.data$Hellinger)),
      .groups = "drop"
    )
}

#' @title suggest_resolution
#' @description
#' Rank clustering resolutions using one of two methods: \code{"rank"}
#' (default) ranks all resolutions by their median metric values, while
#' \code{"local_optima"} first filters to reproducible resolutions (low upper
#' Hellinger CI), then ranks the filtered set using both direct rank
#' aggregation and curvature-based local optima detection.
#'
#' @param cv_results Data frame from \code{\link{clust_opt}} cross-validation.
#' @param method Character; \code{"rank"} (default) ranks all resolutions by
#'   median metric values without filtering, \code{"local_optima"} applies
#'   the Hellinger reproducibility filter and ranks the filtered set by both
#'   rank aggregation and curvature.
#' @param upper_Hell_score_thresh Numeric. Initial threshold on the upper 95\%
#'   Hellinger CI for filtering reproducible resolutions (default 0.1).
#'   Only used when \code{method = "local_optima"}.
#' @param upper_Hell_score_thresh_relaxed Numeric. Relaxed threshold used when
#'   fewer than \code{min_resolutions} pass the initial filter (default 0.2).
#'   Only used when \code{method = "local_optima"}.
#' @param min_resolutions Integer. Minimum number of resolutions required to
#'   pass the initial threshold before relaxing (default 4). Capped at the
#'   total number of resolutions tested. Only used when
#'   \code{method = "local_optima"}.
#'
#' @details
#' For \code{method = "local_optima"}, resolutions are filtered to those with
#' \code{upper_Hell_95ci < upper_Hell_score_thresh}. If fewer than
#' \code{min_resolutions} pass, the threshold is relaxed to
#' \code{upper_Hell_score_thresh_relaxed}. An error is raised if no resolutions
#' pass the relaxed threshold. Curvature values (\code{second_order_diff}) are
#' computed on the full, unfiltered resolution sequence to preserve consecutive
#' spacing, then ranked on only the filtered subset. Both rank aggregation
#' (\code{mean_rank}) and curvature (\code{curvature_mean_rank}) columns are
#' included in the output for comparison.
#'
#' The default thresholds (0.1 / 0.2) have straightforward interpretation given
#' that the Hellinger distance is bounded between 0 and 1.
#'
#' For \code{method = "rank"}, all resolutions are ranked directly by their
#' median metric values with no filtering applied.
#'
#' @return A data frame with one row per resolution (filtered for
#'   \code{"local_optima"}, all for \code{"rank"}), containing summarized
#'   median metrics and \code{upper_Hell_95ci}. For \code{"local_optima"}:
#'   both rank-based (\code{rank_sil}, \code{rank_kl}, \code{rank_hellinger},
#'   \code{rank_modularity}, \code{mean_rank}) and curvature-based
#'   (\code{curvature_rank_sil}, etc., \code{curvature_mean_rank}) columns,
#'   sorted by \code{mean_rank}. For \code{"rank"}: rank-based columns only,
#'   sorted by \code{mean_rank}.
#'
#' @examples
#' \dontrun{
#' cv_results <- clust_opt(seurat_obj, subject_ids = "donor_id")
#' rankings <- suggest_resolution(cv_results)
#' rankings_rank <- suggest_resolution(cv_results, method = "rank")
#' }
#'
#' @seealso \code{\link{clust_opt}}, \code{\link{plot_rank_metrics}},
#'   \code{\link{plot_mean_rank}}
#' @export
suggest_resolution <- function(cv_results,
                               method = c("rank", "local_optima"),
                               upper_Hell_score_thresh = NULL,
                               upper_Hell_score_thresh_relaxed = NULL,
                               min_resolutions = 4) {

  method <- match.arg(method)

  # --- Summarize CV folds per resolution ---
  summ <- summarize_cv_metrics(cv_results)

  # --- Upper confidence bound on Hellinger ---
  summ$upper_Hell_95ci <- summ$median_Hell_score +
    1.96 * summ$standard_error_Hell_score

  if (method == "rank") {
    # --- Rank all resolutions by median metric values ---
    summ$rank_sil        <- rank(-summ$median_score)
    summ$rank_kl         <- rank(summ$median_KLD_score)
    summ$rank_hellinger  <- rank(summ$median_Hell_score)
    summ$rank_modularity <- rank(-summ$median_modularity_score)
    summ$mean_rank <- rowMeans(
      cbind(summ$rank_sil, summ$rank_kl,
            summ$rank_hellinger, summ$rank_modularity)
    )
    summ <- summ[order(summ$mean_rank), ]
    rownames(summ) <- NULL
    return(summ)
  }

  # --- method = "local_optima" ---

  if (is.null(upper_Hell_score_thresh))
    upper_Hell_score_thresh <- 0.1
  if (is.null(upper_Hell_score_thresh_relaxed))
    upper_Hell_score_thresh_relaxed <- 0.2

  # --- Curvature on full (unfiltered) resolution sequence ---
  summ$curv_sil  <- second_order_diff(summ$median_score)
  summ$curv_kl   <- second_order_diff(summ$median_KLD_score)
  summ$curv_hell <- second_order_diff(summ$median_Hell_score)
  summ$curv_mod  <- second_order_diff(summ$median_modularity_score)

  # --- Filter for reproducible resolutions ---
  min_res <- min(min_resolutions, nrow(summ))
  keep <- summ$upper_Hell_95ci < upper_Hell_score_thresh
  if (sum(keep) < min_res)
    keep <- summ$upper_Hell_95ci < upper_Hell_score_thresh_relaxed

  if (sum(keep) == 0) {
    stop("No resolutions passed the Hellinger reproducibility filter ",
         "(upper_Hell_95ci < ", round(upper_Hell_score_thresh_relaxed, 3),
         "). Consider increasing upper_Hell_score_thresh_relaxed or adding ",
         "more subjects to reduce the confidence interval width.")
  }
  summ <- summ[keep, ]

  # --- Rank-based ranking (direct ranking of median values) ---
  summ$rank_sil        <- rank(-summ$median_score)
  summ$rank_kl         <- rank(summ$median_KLD_score)
  summ$rank_hellinger  <- rank(summ$median_Hell_score)
  summ$rank_modularity <- rank(-summ$median_modularity_score)
  summ$mean_rank <- rowMeans(
    cbind(summ$rank_sil, summ$rank_kl,
          summ$rank_hellinger, summ$rank_modularity)
  )

  # --- Curvature-based ranking (rank the pre-filter curvature values) ---
  summ$curvature_rank_sil        <- rank(-summ$curv_sil, na.last = "keep")
  summ$curvature_rank_kl         <- rank(-summ$curv_kl, na.last = "keep")
  summ$curvature_rank_hellinger  <- rank(summ$curv_hell, na.last = "keep")
  summ$curvature_rank_modularity <- rank(-summ$curv_mod, na.last = "keep")
  summ$curvature_mean_rank <- rowMeans(
    cbind(summ$curvature_rank_sil, summ$curvature_rank_kl,
          summ$curvature_rank_hellinger, summ$curvature_rank_modularity),
    na.rm = TRUE
  )

  # --- Clean up intermediate columns and sort ---
  summ$curv_sil  <- NULL
  summ$curv_kl   <- NULL
  summ$curv_hell <- NULL
  summ$curv_mod  <- NULL
  summ <- summ[order(summ$mean_rank), ]
  rownames(summ) <- NULL
  summ
}

#' @title plot_rank_metrics
#' @description
#' Create a publication-quality plot showing individual metric ranks
#' and mean rank across clustering resolutions.
#'
#' @param rank_results A data.frame output from \code{\link{suggest_resolution}}.
#' @param method Character; \code{"rank"} (default) for direct rank aggregation,
#'   or \code{"curvature"} for curvature-based ranking. The curvature method
#'   requires at least 3 resolutions.
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
#' cv_results <- clust_opt(seurat_obj, subject_ids = "donor_id")
#' rankings <- suggest_resolution(cv_results)
#' plot_rank_metrics(rankings)
#' plot_rank_metrics(rankings, method = "curvature")
#' }
#'
#' @seealso \code{\link{suggest_resolution}}, \code{\link{plot_mean_rank}}
#' @export
plot_rank_metrics <- function(rank_results, method = "rank",
                              highlight_best = TRUE, base_size = 12) {
  method <- match.arg(method, c("rank", "curvature"))

  if (method == "curvature" && nrow(rank_results) < 3) {
    stop("At least 3 resolutions are required for curvature-based ranking.")
  }

  # Select column prefix based on method
  if (method == "rank") {
    rank_cols <- c("rank_sil", "rank_kl", "rank_hellinger",
                   "rank_modularity", "mean_rank")
    mean_col <- "mean_rank"
    n_ranks <- nrow(rank_results)
  } else {
    rank_cols <- c("curvature_rank_sil", "curvature_rank_kl",
                   "curvature_rank_hellinger", "curvature_rank_modularity",
                   "curvature_mean_rank")
    mean_col <- "curvature_mean_rank"
    n_ranks <- nrow(rank_results) - 2
  }

  # Reshape data for plotting
  plot_data <- rank_results |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(rank_cols),
      names_to = "metric",
      values_to = "rank"
    ) |>
    dplyr::mutate(
      resolution = as.factor(.data$resolution),
      metric = factor(.data$metric,
        levels = rank_cols,
        labels = c("Silhouette", "KL Divergence", "Hellinger",
                   "Modularity", "Mean Rank")
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

  linewidths <- c(
    "Silhouette" = 0.8, "KL Divergence" = 0.8,
    "Hellinger" = 0.8, "Modularity" = 0.8, "Mean Rank" = 1.2
  )

  point_sizes <- c(
    "Silhouette" = 2.5, "KL Divergence" = 2.5,
    "Hellinger" = 2.5, "Modularity" = 2.5, "Mean Rank" = 3.5
  )

  alphas <- c(
    "Silhouette" = 0.7, "KL Divergence" = 0.7,
    "Hellinger" = 0.7, "Modularity" = 0.7, "Mean Rank" = 1.0
  )

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
    ggplot2::scale_y_reverse(breaks = seq_len(n_ranks)) +
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
    best_idx <- which.min(rank_results[[mean_col]])
    best_res <- as.factor(rank_results$resolution[best_idx])
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
#' @param method Character; \code{"rank"} (default) for direct rank aggregation,
#'   or \code{"curvature"} for curvature-based ranking. The curvature method
#'   requires at least 3 resolutions.
#' @param highlight_best Logical; if TRUE (default), highlight the best resolution
#'   with a vertical dashed line.
#' @param base_size Numeric; base font size for theme. Default is 12.
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' cv_results <- clust_opt(seurat_obj, subject_ids = "donor_id")
#' rankings <- suggest_resolution(cv_results)
#' plot_mean_rank(rankings)
#' plot_mean_rank(rankings, method = "curvature")
#' }
#'
#' @seealso \code{\link{suggest_resolution}}, \code{\link{plot_rank_metrics}}
#' @export
plot_mean_rank <- function(rank_results, method = "rank",
                           highlight_best = TRUE, base_size = 12) {
  method <- match.arg(method, c("rank", "curvature"))

  if (method == "curvature" && nrow(rank_results) < 3) {
    stop("At least 3 resolutions are required for curvature-based ranking.")
  }

  mean_col <- if (method == "rank") "mean_rank" else "curvature_mean_rank"
  n_ranks <- if (method == "rank") nrow(rank_results) else nrow(rank_results) - 2

  plot_data <- rank_results |>
    dplyr::mutate(resolution = as.factor(.data$resolution))

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$resolution, y = .data[[mean_col]], group = 1)
  ) +
    ggplot2::geom_line(linewidth = 1, color = "#0072B2") +
    ggplot2::geom_point(size = 4, color = "#0072B2") +
    ggplot2::scale_y_reverse(breaks = seq_len(n_ranks)) +
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
    best_idx <- which.min(rank_results[[mean_col]])
    best_res <- as.factor(rank_results$resolution[best_idx])
    p <- p + ggplot2::geom_vline(
      xintercept = best_res,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    )
  }

  p
}
