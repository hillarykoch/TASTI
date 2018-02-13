\name{bLens}
\alias{bLens}
\docType{data}
\title{
A data set of observed gene trees at 1,000 loci
}
\description{
This is a matrix of genee trees observed at 1,000 loci for 3 sampled individuals. The data are formatted for use with \code{\link{buildST}}.
}
\usage{data("bLens")}
\format{
  A data frame with 1000 observations on the following 5 variables.
  \describe{
    \item{\samp{((AB)C)}}{1 if ((AB)C) is the observed topology, 0 otherwise}
    \item{\samp{((BC)A)}}{1 if ((BC)A) is the observed topology, 0 otherwise}
    \item{\samp{((AC)B)}}{1 if ((AC)B) is the observed topology, 0 otherwise}
    \item{\code{t1}}{First coalescence time, going backward in time}
    \item{\code{t2}}{Second coalescence time, going backward in time}
  }
}
\examples{
data(bLens)
\dontrun{
fit <- buildST(bLens, theta = .005, N = 25000)
fit
}
}
\keyword{datasets}
