#' @include utils.R
#' @importFrom rlang .data
NULL


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
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
  sil_summary <- sil_summary(sil_dist)

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
    sil_summary,
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
    sil_summary,
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
