\name{delaunay_graph}
\alias{delaunay_graph}
\title{Build a Delaunay triangulation adjacency graph}
\description{
Constructs a directional sparse adjacency matrix from Delaunay triangulation edges.
}
\usage{
delaunay_graph(coords)
}
\arguments{
\item{coords}{Numeric matrix/data.frame with two columns.}
}
\value{A sparse directional adjacency matrix (\code{dgTMatrix}).}
