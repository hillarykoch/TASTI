\name{topos}
\alias{topos}
\docType{data}
\title{
A data set of observed gene tree topologies at 1,000 loci
}
\description{
This is a vector of gene tree topologies observed at 1,000 loci for 3 sampled individuals. The data are formatted for use with \code{\link{buildST}}.
}
\usage{data("topos")}


\examples{
data(topos)
\dontrun{
fit <- buildST(topos, theta = .005, N = 25000, gridSearch = F, branchLengths = F)
fit
}
}
\keyword{datasets}
