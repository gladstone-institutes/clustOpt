#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title split_pca_dimensions
#' @description
#' Takes a Seurat object with an existing PCA reduction and splits the
#' dimensions into 2 sets A and B. The split_method is either "odd_even" or
#' "var_balanced". The last PC will be removed if there are an odd number of
#' PCs. The new reductions will be named "A_pca" and "B_pca" with A being odd
#' and B being even for split_method = "odd_even".
#'
#' @param input Seurat object
#' @param split_method How should the PCs be split into 2, options are odd even
#' split ("odd_even") or balancing variation explained between 2 sets
#' ("var_balanced").
#' @param verbose Ouput messages
#'
#' @return Seurat object with new PCA reductions
#' @export
#' @import Seurat
split_pca_dimensions <- function(input,
                                 split_method = "odd_even",
                                 verbose = FALSE) {
  if (!inherits(input, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (!"pca" %in% names(input@reductions)) {
    stop("Seurat object must have a PCA reduction")
  }

  if (split_method == "odd_even") {
    if (verbose) {
      message("Splitting PCs into 2 sets odd (A) and even (B)")
    }
    pca <- input@reductions$pca
    dims <- ncol(pca@cell.embeddings)

    even_dims <- seq(2, dims, by = 2)
    odd_dims <- seq(1, dims, by = 2)

    even_pca <- pca
    even_pca@cell.embeddings <- even_pca@cell.embeddings[, even_dims]
    even_pca@feature.loadings <- even_pca@feature.loadings[, even_dims]
    even_pca@key <- "B_pca_"

    odd_pca <- pca
    odd_pca@cell.embeddings <- odd_pca@cell.embeddings[, odd_dims]
    odd_pca@feature.loadings <- odd_pca@feature.loadings[, odd_dims]
    odd_pca@key <- "A_pca_"

    input@reductions$B_pca <- even_pca
    input@reductions$A_pca <- odd_pca

    return(input)
  } else if (split_method == "var_balanced") {
    if (verbose) {
      message(
        "Splitting PCs into 2 sets with approximately equal total",
        " variance explained"
      )
    }
    pca <- input@reductions$pca
    var_expl <- pca@stdev^2 # Variance explained by each PC
    dims <- seq_len(length(var_expl))

    # Sort PCs by variance explained in descending order
    order_desc <- order(var_expl, decreasing = TRUE)

    # Greedy partition: assign each PC to the set with smaller sum of variance
    set1 <- numeric(0)
    set2 <- numeric(0)
    sum1 <- 0
    sum2 <- 0

    for (i in order_desc) {
      if (sum1 <= sum2) {
        set1 <- c(set1, i)
        sum1 <- sum1 + var_expl[i]
      } else {
        set2 <- c(set2, i)
        sum2 <- sum2 + var_expl[i]
      }
    }

    # Reorder each set ascending
    set1 <- sort(set1)
    set2 <- sort(set2)

    # Create new reductions
    pca1 <- pca
    pca2 <- pca

    pca1@cell.embeddings <- pca1@cell.embeddings[, set1, drop = FALSE]
    pca1@feature.loadings <- pca1@feature.loadings[, set1, drop = FALSE]
    pca1@stdev <- pca1@stdev[set1]
    pca1@key <- "A_pca_"

    pca2@cell.embeddings <- pca2@cell.embeddings[, set2, drop = FALSE]
    pca2@feature.loadings <- pca2@feature.loadings[, set2, drop = FALSE]
    pca2@stdev <- pca2@stdev[set2]
    pca2@key <- "B_pca_"

    # Store new reductions in the Seurat object
    input@reductions$A_pca <- pca1
    input@reductions$B_pca <- pca2

    # Report how close the two sets are in total variance
    total_var <- sum(var_expl)
    if (verbose) {
      message(
        "Set A ", ncol(input@reductions$A_pca), " PCs", "\n",
        "Set A variance: ", round(sum1 / total_var * 100, 2), "%\n",
        "Set B ", ncol(input@reductions$B_pca), " PCs", "\n",
        "Set B variance: ", round(sum2 / total_var * 100, 2), "%"
      )
    }
    return(input)
  }
}
