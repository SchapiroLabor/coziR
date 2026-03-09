#' Build a k-nearest-neighbor adjacency graph
#'
#' @param coords Numeric matrix/data.frame with two columns (`x`, `y`) and one row per cell.
#' @param n_neighbors Number of nearest neighbors per source cell.
#'
#' @return A sparse adjacency matrix (`dgTMatrix`) with directional edges.
#' @export
knn_graph <- function(coords, n_neighbors = 6L) {
  coords <- .validate_coords(coords)
  n_neighbors <- as.integer(n_neighbors)

  if (n_neighbors < 1L) {
    stop("`n_neighbors` must be >= 1.", call. = FALSE)
  }

  n_cells <- nrow(coords)
  if (n_cells == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
      dims = c(0L, 0L), repr = "T"
    ))
  }

  # FNN::get.knn() already excludes self-neighbors.
  k <- min(n_neighbors, max(n_cells - 1L, 0L))
  if (k <= 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_cells, n_cells), repr = "T"
    ))
  }

  nn <- FNN::get.knn(coords, k = k)$nn.index

  rows <- rep.int(seq_len(n_cells), times = k)
  cols <- as.integer(nn)

  Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = 1L,
    dims = c(n_cells, n_cells),
    repr = "T"
  )
}

#' Build a radius-based adjacency graph
#'
#' @param coords Numeric matrix/data.frame with two columns (`x`, `y`) and one row per cell.
#' @param radius Radius threshold for neighborhood inclusion.
#'
#' @return A sparse adjacency matrix (`dgTMatrix`) with directional edges.
#' @export
radius_graph <- function(coords, radius) {
  coords <- .validate_coords(coords)
  radius <- as.numeric(radius)

  if (length(radius) != 1L || is.na(radius) || radius <= 0) {
    stop("`radius` must be a single positive number.", call. = FALSE)
  }

  n_cells <- nrow(coords)
  if (n_cells == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
      dims = c(0L, 0L), repr = "T"
    ))
  }

  dist_m <- as.matrix(dist(coords))
  hits <- which(dist_m <= radius & dist_m > 0, arr.ind = TRUE)

  if (nrow(hits) == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_cells, n_cells), repr = "T"
    ))
  }

  Matrix::sparseMatrix(
    i = hits[, 1L],
    j = hits[, 2L],
    x = 1L,
    dims = c(n_cells, n_cells),
    repr = "T"
  )
}

#' Build a Delaunay triangulation adjacency graph
#'
#' @param coords Numeric matrix/data.frame with two columns (`x`, `y`) and one row per cell.
#'
#' @return A sparse adjacency matrix (`dgTMatrix`) with directional edges.
#' @export
delaunay_graph <- function(coords) {
  coords <- .validate_coords(coords)
  n_cells <- nrow(coords)

  if (n_cells == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
      dims = c(0L, 0L), repr = "T"
    ))
  }

  if (n_cells < 3L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_cells, n_cells), repr = "T"
    ))
  }

  tri <- .delaunay_triangles(coords)

  if (is.null(tri) || nrow(tri) == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_cells, n_cells), repr = "T"
    ))
  }

  pairs <- rbind(
    tri[, c(1L, 2L), drop = FALSE],
    tri[, c(2L, 1L), drop = FALSE],
    tri[, c(2L, 3L), drop = FALSE],
    tri[, c(3L, 2L), drop = FALSE],
    tri[, c(1L, 3L), drop = FALSE],
    tri[, c(3L, 1L), drop = FALSE]
  )

  pairs <- unique(pairs)

  Matrix::sparseMatrix(
    i = pairs[, 1L],
    j = pairs[, 2L],
    x = 1L,
    dims = c(n_cells, n_cells),
    repr = "T"
  )
}

.delaunay_triangles <- function(coords) {
  n <- nrow(coords)
  if (n < 3L) {
    return(matrix(integer(0), nrow = 0L, ncol = 3L))
  }

  coords <- .jitter_duplicate_points(coords)

  x <- coords[, 1L]
  y <- coords[, 2L]

  min_x <- min(x)
  min_y <- min(y)
  max_x <- max(x)
  max_y <- max(y)
  dx <- max_x - min_x
  dy <- max_y - min_y
  dmax <- max(dx, dy)
  mid_x <- (min_x + max_x) / 2
  mid_y <- (min_y + max_y) / 2

  super <- matrix(c(
    mid_x - 20 * dmax, mid_y - dmax,
    mid_x, mid_y + 20 * dmax,
    mid_x + 20 * dmax, mid_y - dmax
  ), ncol = 2, byrow = TRUE)

  all_pts <- rbind(coords, super)
  s1 <- n + 1L
  s2 <- n + 2L
  s3 <- n + 3L

  triangles <- matrix(c(s1, s2, s3), nrow = 1L)

  for (p in seq_len(n)) {
    bad <- logical(nrow(triangles))
    for (t in seq_len(nrow(triangles))) {
      bad[t] <- .circumcircle_contains(
        point = all_pts[p, ],
        tri_pts = all_pts[triangles[t, ], , drop = FALSE]
      )
    }

    bad_idx <- which(bad)
    if (length(bad_idx) == 0L) {
      next
    }

    bad_tris <- triangles[bad_idx, , drop = FALSE]
    edges <- rbind(
      bad_tris[, c(1L, 2L), drop = FALSE],
      bad_tris[, c(2L, 3L), drop = FALSE],
      bad_tris[, c(3L, 1L), drop = FALSE]
    )
    edges <- t(apply(edges, 1, sort))
    edge_key <- paste(edges[, 1L], edges[, 2L], sep = "-")
    edge_count <- table(edge_key)
    boundary_keys <- names(edge_count[edge_count == 1L])
    boundary_edges <- edges[edge_key %in% boundary_keys, , drop = FALSE]

    triangles <- triangles[!bad, , drop = FALSE]

    if (nrow(boundary_edges) > 0L) {
      new_tris <- cbind(boundary_edges, p)
      triangles <- rbind(triangles, new_tris)
    }
  }

  keep <- triangles[, 1L] <= n & triangles[, 2L] <= n & triangles[, 3L] <= n
  triangles <- triangles[keep, , drop = FALSE]
  if (nrow(triangles) == 0L) {
    return(matrix(integer(0), nrow = 0L, ncol = 3L))
  }
  unique(triangles)
}

.circumcircle_contains <- function(point, tri_pts) {
  ax <- tri_pts[1L, 1L]
  ay <- tri_pts[1L, 2L]
  bx <- tri_pts[2L, 1L]
  by <- tri_pts[2L, 2L]
  cx <- tri_pts[3L, 1L]
  cy <- tri_pts[3L, 2L]
  px <- point[1L]
  py <- point[2L]

  d <- 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
  if (abs(d) < .Machine$double.eps) {
    return(FALSE)
  }

  ux <- ((ax^2 + ay^2) * (by - cy) + (bx^2 + by^2) * (cy - ay) +
    (cx^2 + cy^2) * (ay - by)) / d
  uy <- ((ax^2 + ay^2) * (cx - bx) + (bx^2 + by^2) * (ax - cx) +
    (cx^2 + cy^2) * (bx - ax)) / d

  r2 <- (ux - ax)^2 + (uy - ay)^2
  p2 <- (ux - px)^2 + (uy - py)^2

  p2 <= r2 + 1e-12
}

.jitter_duplicate_points <- function(coords) {
  dup_idx <- which(duplicated(coords) | duplicated(coords, fromLast = TRUE))
  if (length(dup_idx) == 0L) {
    return(coords)
  }

  eps <- 1e-8
  offset <- seq_along(dup_idx) * eps
  coords[dup_idx, 1L] <- coords[dup_idx, 1L] + offset
  coords[dup_idx, 2L] <- coords[dup_idx, 2L] + offset
  coords
}
