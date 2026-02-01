test_that("compute_knn_graph handles large datasets (500k cells)", {
  skip_on_cran()
  skip_on_ci()

  n <- 500000
  set.seed(42)
  coords <- matrix(rnorm(n * 3), ncol = 3)

  adj_matrix <- compute_knn_graph(coords, k = 20, distance_metric = "euclidean")

  expect_s4_class(adj_matrix, "sparseMatrix")
  expect_equal(nrow(adj_matrix), n)
  expect_equal(ncol(adj_matrix), n)

  expect_true(all(Matrix::rowSums(adj_matrix > 0) >= 20))
})

test_that("calculate_modularity_from_coords works with large datasets (500k cells)", {
  skip_on_cran()
  skip_on_ci()

  n <- 500000
  set.seed(42)
  coords <- matrix(rnorm(n * 3), ncol = 3)
  clusters <- sample(1:10, n, replace = TRUE)

  result <- calculate_modularity_from_coords(clusters, coords, k = 20)

  expect_type(result, "list")
  expect_true("modularity" %in% names(result))
  expect_true(is.numeric(result$modularity))
  expect_true(!is.na(result$modularity))
})
