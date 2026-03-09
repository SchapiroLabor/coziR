\name{radius_graph}
\alias{radius_graph}
\title{Build a radius-based adjacency graph}
\description{
Constructs a directional sparse adjacency matrix for all pairs within a radius.
}
\usage{
radius_graph(coords, radius)
}
\arguments{
\item{coords}{Numeric matrix/data.frame with two columns.}
\item{radius}{Radius threshold for neighborhood inclusion.}
}
\value{A sparse directional adjacency matrix (\code{dgTMatrix}).}
