# Large data tests for modularity via Seurat's SNN graph are covered
# by the integration tests in test-metrics.R (calculate_modularity).
# The compute_knn_graph and calculate_modularity_from_coords functions
# were replaced by Seurat::FindNeighbors for SNN graph construction.
