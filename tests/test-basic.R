library(coziR)

cat("Running coziR base tests...\n")

coords <- matrix(c(
  0, 0,
  1, 0,
  0, 1,
  1, 1,
  0.5, 0.5,
  1.5, 0.5
), ncol = 2, byrow = TRUE)
labels <- c("A", "A", "B", "B", "A", "B")

# Test 1: KNN graph has no self-edges and expected out-degree
k <- 3L
adj_knn <- knn_graph(coords, n_neighbors = k)
row_degree <- tabulate(adj_knn@i + 1L, nbins = nrow(coords))
stopifnot(all(row_degree == min(k, nrow(coords) - 1L)))
stopifnot(all((adj_knn@i + 1L) != (adj_knn@j + 1L)))

# Test 2: End-to-end run_cozi output structure
res <- run_cozi(
  coords = coords,
  labels = labels,
  nbh_def = "knn",
  n_neighbors = 3,
  n_permutations = 20,
  random_state = 1
)
stopifnot(is.list(res))
stopifnot(all(c("cond_ratio", "zscore") %in% names(res)))
stopifnot(is.data.frame(res$cond_ratio))
stopifnot(is.data.frame(res$zscore))
stopifnot(identical(dim(res$cond_ratio), dim(res$zscore)))

# Test 3: n_permutations = 1 is numerically stable (no NA from SD calculation)
res_one <- run_cozi(
  coords = coords,
  labels = labels,
  nbh_def = "knn",
  n_neighbors = 3,
  n_permutations = 1,
  random_state = 1
)
stopifnot(!any(is.na(as.matrix(res_one$zscore))))

# Test 4: fixed_type by label works
res_fixed <- run_cozi(
  coords = coords,
  labels = labels,
  nbh_def = "radius",
  radius = 1.5,
  n_permutations = 10,
  random_state = 1,
  fixed_type = "A"
)
stopifnot(all(dim(res_fixed$zscore) == c(2, 2)))

# Test 5: Dotplot returns ggplot object when ggplot2 is available
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- cozi_dotplot(res)
  stopifnot(inherits(p, "ggplot"))
} else {
  cat("Skipping dotplot test: ggplot2 not installed.\n")
}

cat("All coziR base tests passed.\n")
