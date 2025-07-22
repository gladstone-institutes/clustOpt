#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title split_pca_dimensions
#' @description
#' Takes a Seurat object with an existing PCA reduction and splits the
#' dimensions into 2 sets odd and even PCs.
#'
#' @param input Seurat object
#' @param verbose Output messages
#'
#' @return Seurat object with new PCA reductions
#' @export
#' @import Seurat
split_pca_dimensions <- function(input,
                                 verbose = FALSE) {
  if (!inherits(input, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (!"pca" %in% names(input@reductions)) {
    stop("Seurat object must have a PCA reduction")
  }


  if (verbose) {
    message("Splitting PCs into 2 sets odd and even")
  }
  pca <- input@reductions$pca
  dims <- ncol(pca@cell.embeddings)

  even_dims <- seq(2, dims, by = 2)
  odd_dims <- seq(1, dims, by = 2)

  even_pca <- pca
  even_pca@cell.embeddings <- even_pca@cell.embeddings[, even_dims, drop = FALSE]
  even_pca@feature.loadings <- even_pca@feature.loadings[, even_dims, drop = FALSE]
  even_pca@stdev <- even_pca@stdev[even_dims]
  even_pca@key <- "even_pca_"

  odd_pca <- pca
  odd_pca@cell.embeddings <- odd_pca@cell.embeddings[, odd_dims, drop = FALSE]
  odd_pca@feature.loadings <- odd_pca@feature.loadings[, odd_dims, drop = FALSE]
  odd_pca@stdev <- odd_pca@stdev[odd_dims]
  odd_pca@key <- "odd_pca_"

  input@reductions$even_pca <- even_pca
  input@reductions$odd_pca <- odd_pca
  return(input)
}
