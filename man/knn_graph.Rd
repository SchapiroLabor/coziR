\name{knn_graph}
\alias{knn_graph}
\title{Build a k-nearest-neighbor adjacency graph}
\description{
Constructs a directional sparse adjacency matrix from k-nearest neighbors.
}
\usage{
knn_graph(coords, n_neighbors = 6L)
}
\arguments{
\item{coords}{Numeric matrix/data.frame with two columns.}
\item{n_neighbors}{Number of nearest neighbors per source cell.}
}
\value{A sparse directional adjacency matrix (\code{dgTMatrix}).}
