\name{cozi_dotplot}
\alias{cozi_dotplot}
\title{Dotplot for COZI Results}
\description{
Creates a dotplot where color encodes z-score and dot size encodes conditional cell ratio.
}
\usage{
cozi_dotplot(
  results,
  title = "COZI Dotplot",
  low = "#2166AC",
  mid = "white",
  high = "#B2182B",
  size_limits = c(0, 1)
)
}
\arguments{
\item{results}{Output list from \code{run_cozi()} or \code{nep_analysis()}.}
\item{title}{Plot title.}
\item{low,mid,high}{Low/mid/high colors for z-score scale.}
\item{size_limits}{Numeric length-2 vector for conditional ratio size scale limits.}
}
\value{A \code{ggplot2} plot object.}
