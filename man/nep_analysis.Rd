\name{nep_analysis}
\alias{nep_analysis}
\title{Neighborhood Enrichment Preference Analysis}
\description{
Computes conditional normalized interaction scores and permutation z-scores.
}
\usage{
nep_analysis(
  adj,
  labels,
  n_permutations = 1000L,
  random_state = NULL,
  return_df = TRUE,
  fixed_type = NULL
)
}
\arguments{
\item{adj}{Adjacency graph matrix.}
\item{labels}{Vector of cell type labels.}
\item{n_permutations}{Number of permutations.}
\item{random_state}{Optional integer seed.}
\item{return_df}{Return data.frames if \code{TRUE}.}
\item{fixed_type}{Optional label name or index to keep fixed during permutations.}
}
\value{List with \code{cond_ratio} and \code{zscore}.}
