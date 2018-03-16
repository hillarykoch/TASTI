\name{get_triples}
\alias{get_triples}

\title{
Convert newick format file into format needed for TASTI
}
\description{
This file reads in a file of rooted triples in newick format, and converts it into a form recognized by TASTI.
}
\usage{
get_triples(file, branchLengths = TRUE)
}
\arguments{
  \item{file}{Direct path to newick tree file.}
  \item{branchLengths}{Logical. Should the format returned be one for inference with branch lengths? Defaults to true, but if FALSE, returns the form needed to make inference on topologies alone is used.}
}

\value{
A matrix (if branchLengths = TRUE) or a vector (if branchLengths = FALSE) than can be passed to \code{\link{buildST}}
}
\author{\email{hbk5086@psu.edu}}
\examples{
file <- system.file("extdata", "trees.txt", package = "TASTI")
data <- get_triples(file, branchLengths = TRUE)
head(data)
}
