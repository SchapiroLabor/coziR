\name{run_cozi}
\alias{run_cozi}
\title{Run COZI workflow from coordinates and labels}
\description{
Runs graph construction and NEP analysis in one function call.
}
\usage{
run_cozi(
  coords,
  labels,
  nbh_def = "knn",
  n_neighbors = 6L,
  radius = 0.2,
  n_permutations = 100L,
  random_state = NULL,
  fixed_type = NULL
)
}
\arguments{
\item{coords}{Numeric matrix/data.frame with two columns.}
\item{labels}{Vector of cell type labels.}
\item{nbh_def}{Neighborhood definition: \code{"knn"}, \code{"radius"}, or \code{"delaunay"}.}
\item{n_neighbors}{Neighbors used when \code{nbh_def = "knn"}.}
\item{radius}{Radius used when \code{nbh_def = "radius"}.}
\item{n_permutations}{Number of permutations.}
\item{random_state}{Optional integer seed.}
\item{fixed_type}{Optional label kept fixed during permutations.}
}
\value{List with \code{cond_ratio} and \code{zscore}.}
