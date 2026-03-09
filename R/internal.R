.validate_coords <- function(coords) {
  if (is.data.frame(coords)) {
    coords <- as.matrix(coords)
  }

  if (!is.matrix(coords) || !is.numeric(coords)) {
    stop("`coords` must be a numeric matrix or data.frame.", call. = FALSE)
  }

  if (ncol(coords) != 2L) {
    stop("`coords` must have exactly 2 columns (x, y).", call. = FALSE)
  }

  coords
}

.encode_labels <- function(labels) {
  labels <- as.character(labels)
  labels_factor <- factor(labels)

  list(
    labels_int = as.integer(labels_factor),
    label_names = levels(labels_factor)
  )
}

.pop_sd <- function(x) {
  m <- mean(x)
  sqrt(mean((x - m)^2))
}

.compute_pair_counts_and_denominators <- function(adj, labels_int, n_types) {
  n_cells <- length(labels_int)
  if (n_cells == 0L || n_types == 0L) {
    return(list(
      counts = matrix(0, nrow = n_types, ncol = n_types),
      denom = matrix(0, nrow = n_types, ncol = n_types)
    ))
  }

  adj <- .as_dgT_adj(adj, n_cells)

  if (length(adj@i) == 0L) {
    return(list(
      counts = matrix(0, nrow = n_types, ncol = n_types),
      denom = matrix(0, nrow = n_types, ncol = n_types)
    ))
  }

  src <- adj@i + 1L
  tgt <- adj@j + 1L

  src_labels <- labels_int[src]
  tgt_labels <- labels_int[tgt]

  pair_idx <- (src_labels - 1L) * n_types + tgt_labels
  counts <- matrix(
    tabulate(pair_idx, nbins = n_types * n_types),
    nrow = n_types,
    byrow = TRUE
  )

  cell_has <- Matrix::sparseMatrix(
    i = src,
    j = tgt_labels,
    x = 1L,
    dims = c(n_cells, n_types)
  )
  cell_has@x[] <- 1

  denom <- matrix(0, nrow = n_types, ncol = n_types)
  for (a in seq_len(n_types)) {
    mask_a <- labels_int == a
    if (any(mask_a)) {
      denom[a, ] <- Matrix::colSums(cell_has[mask_a, , drop = FALSE] > 0)
    }
  }

  list(counts = counts, denom = denom)
}

.as_dgT_adj <- function(adj, n_cells) {
  if (inherits(adj, "dgTMatrix")) {
    out <- adj
  } else if (inherits(adj, "sparseMatrix")) {
    out <- as(adj, "dgTMatrix")
  } else if (is.matrix(adj)) {
    nz <- which(adj != 0, arr.ind = TRUE)
    if (nrow(nz) == 0L) {
      out <- Matrix::sparseMatrix(
        i = integer(0), j = integer(0), x = numeric(0),
        dims = c(n_cells, n_cells), repr = "T"
      )
    } else {
      out <- Matrix::sparseMatrix(
        i = nz[, 1L], j = nz[, 2L], x = as.numeric(adj[nz]),
        dims = c(n_cells, n_cells), repr = "T"
      )
    }
  } else {
    stop("`adj` must be a matrix or sparse Matrix object.", call. = FALSE)
  }

  if (!all(dim(out) == c(n_cells, n_cells))) {
    stop("`adj` dimensions must match number of labels.", call. = FALSE)
  }

  out
}
