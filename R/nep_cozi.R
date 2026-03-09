#' Neighborhood Enrichment Preference Analysis
#'
#' Computes observed conditional neighborhood preference values, permutation-based
#' expected values, z-scores, and conditional cell ratios.
#'
#' @param adj Adjacency graph (`dgTMatrix`, other sparse matrix, or dense matrix).
#' @param labels Vector of cell type labels, one per cell.
#' @param n_permutations Number of permutations.
#' @param random_state Optional integer seed.
#' @param return_df Logical; if `TRUE`, returns data.frames with label names.
#' @param fixed_type Optional label name or integer label index to keep fixed during permutations.
#'
#' @return A list with `cond_ratio` and `zscore`.
#' @export
nep_analysis <- function(adj,
                         labels,
                         n_permutations = 1000L,
                         random_state = NULL,
                         return_df = TRUE,
                         fixed_type = NULL) {
  if (!is.null(random_state)) {
    set.seed(as.integer(random_state))
  }

  encoded <- .encode_labels(labels)
  labels_int <- encoded$labels_int
  label_names <- encoded$label_names
  n_types <- length(label_names)

  if (n_types == 0L) {
    if (isTRUE(return_df)) {
      empty <- data.frame()
      return(list(cond_ratio = empty, zscore = empty))
    }
    return(list(
      cond_ratio = matrix(numeric(0), nrow = 0L, ncol = 0L),
      zscore = matrix(numeric(0), nrow = 0L, ncol = 0L)
    ))
  }

  comp <- .compute_pair_counts_and_denominators(adj, labels_int, n_types)
  obs_counts <- comp$counts
  obs_denom <- comp$denom
  obs_norm <- obs_counts / pmax(obs_denom, 1)

  cond_ratio <- matrix(0, nrow = n_types, ncol = n_types)
  for (a in seq_len(n_types)) {
    total_a <- sum(labels_int == a)
    cond_ratio[a, ] <- obs_denom[a, ] / max(total_a, 1L)
  }

  n_permutations <- as.integer(n_permutations)
  if (n_permutations < 1L) {
    stop("`n_permutations` must be >= 1.", call. = FALSE)
  }

  fixed_type_idx <- NULL
  if (!is.null(fixed_type)) {
    if (is.character(fixed_type)) {
      if (!fixed_type %in% label_names) {
        stop("`fixed_type` (character) not found in labels.", call. = FALSE)
      }
      fixed_type_idx <- match(fixed_type, label_names)
    } else {
      fixed_type_idx <- as.integer(fixed_type)
      if (length(fixed_type_idx) != 1L || is.na(fixed_type_idx) ||
        fixed_type_idx < 1L || fixed_type_idx > n_types) {
        stop("`fixed_type` (integer) must be in 1..number of label types.", call. = FALSE)
      }
    }
  }

  perm_norm <- array(0, dim = c(n_permutations, n_types, n_types))

  for (k in seq_len(n_permutations)) {
    if (is.null(fixed_type_idx)) {
      perm <- sample(labels_int)
    } else {
      fixed_mask <- labels_int == fixed_type_idx
      perm <- labels_int
      perm[!fixed_mask] <- sample(labels_int[!fixed_mask])
    }

    perm_comp <- .compute_pair_counts_and_denominators(adj, perm, n_types)
    perm_norm[k, , ] <- perm_comp$counts / pmax(perm_comp$denom, 1)
  }

  expected <- apply(perm_norm, c(2, 3), mean)
  std <- apply(perm_norm, c(2, 3), .pop_sd) + 1e-6
  z <- (obs_norm - expected) / std

  dimnames(cond_ratio) <- list(label_names, label_names)
  dimnames(z) <- list(label_names, label_names)

  if (isTRUE(return_df)) {
    return(list(
      cond_ratio = as.data.frame(cond_ratio, check.names = FALSE),
      zscore = as.data.frame(z, check.names = FALSE)
    ))
  }

  list(cond_ratio = cond_ratio, zscore = z)
}

#' Run COZI workflow from coordinates and labels
#'
#' @param coords Numeric matrix/data.frame with two columns (`x`, `y`).
#' @param labels Vector of cell type labels.
#' @param nbh_def Neighborhood definition: `"knn"`, `"radius"`, or `"delaunay"`.
#' @param n_neighbors Number of neighbors for `nbh_def = "knn"`.
#' @param radius Radius threshold for `nbh_def = "radius"`.
#' @param n_permutations Number of permutations.
#' @param random_state Optional integer seed.
#' @param fixed_type Optional label kept fixed during permutations.
#'
#' @return List with `cond_ratio` and `zscore`.
#' @export
run_cozi <- function(coords,
                     labels,
                     nbh_def = "knn",
                     n_neighbors = 6L,
                     radius = 0.2,
                     n_permutations = 100L,
                     random_state = NULL,
                     fixed_type = NULL) {
  nbh_def <- match.arg(nbh_def, choices = c("knn", "radius", "delaunay"))

  adj <- switch(nbh_def,
    knn = knn_graph(coords, n_neighbors = n_neighbors),
    radius = radius_graph(coords, radius = radius),
    delaunay = delaunay_graph(coords)
  )

  nep_analysis(
    adj = adj,
    labels = labels,
    n_permutations = n_permutations,
    random_state = random_state,
    return_df = TRUE,
    fixed_type = fixed_type
  )
}
