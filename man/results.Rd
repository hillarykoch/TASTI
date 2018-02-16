\name{results}
\alias{results}
\docType{data}
\title{
Pre-computed results to sample analysis in the vignette
}
\description{
Computing this algorithm on 1,000 gene trees can take a while. Here are pre-computed results to demonstrate more thoroughly the use of the package.
}
\usage{data("results")}
\format{
This is a list of 3 model outputs. The first comes from doing analysis with only topologies, the second with topologies and branch lengths using numerical optimization, and the third with topologies and branch lengths using a grid search.
}
\examples{
data(results)
str(results)
}
\keyword{datasets}
