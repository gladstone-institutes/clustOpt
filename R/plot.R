#' @include utils.R
#'
NULL


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Silhouette score distribution plots
#'
#' @param sil_dist the output of clustOpt
#' @return A list of ggplot objects
#' @export
#'
create_sil_plots <- function(sil_dist) {
  sil_summary <- sil_summary(sil_dist)

  plot1 <- sil_dist %>%
    ggplot2::ggplot(ggplot2::aes(x = as.factor(resolution), y = avg_width)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Resolution",
                  y = "Avg. Silhouette Score Across All Cells")

  plot2 <- sil_dist %>%
    ggplot2::ggplot(ggplot2::aes(x = as.factor(resolution),
                                 y = cluster_avg_widths)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Resolution", y = "Avg. Silhouette Score Across Clusters")

  plot3 <- ggplot2::ggplot(
    sil_summary,
    ggplot2::aes(x = as.factor(resolution), y = mean_score, group = 1)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_score - (1.96 * standard_error_score),
        ymax = mean_score + (1.96 * standard_error_score),
        width = .3
      ),
      color = "red"
    ) +
    ggplot2::geom_point(colour = "#619CFF") +
    ggplot2::geom_line(colour = "#619CFF") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Resolution", y = "Avg. Silhouette Score Across All Cells")

  plot4 <- ggplot2::ggplot(
    sil_summary,
    ggplot2::aes(x = as.factor(resolution), y = cluster_mean_score, group = 1) # scale should be the res range
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = cluster_mean_score - (1.96 * standard_error_score),
        ymax = cluster_mean_score + (1.96 * standard_error_score),
        width = .3
      ),
      color = "red"
    ) +
    ggplot2::geom_point(colour = "#619CFF") +
    ggplot2::geom_line(colour = "#619CFF") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Resolution", y = "Avg. Silhouette Score Across Clusters")

  list(plot1, plot2, plot3, plot4)
}
