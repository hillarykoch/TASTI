\name{buildST}
\alias{buildST}
\title{
Build a species tree from gene tree data
}
\description{
Construct the most likely species tree from gene trees with or without branch lengths, using either \link{optim} or a grid search.
}
\usage{
buildST(data, theta, N = NULL, mBins = 20, tauBins = 20,
                    gridSearch = TRUE, branchLengths = TRUE, pDiscord = .05,
                    method = "Nelder-Mead", startVals = NULL)
}
\arguments{
  \item{data}{Matrix of gene tree triplet information, with or without branch lengths}
  \item{theta}{Population-scaled mutation rate}
  \item{N}{Effective population size}
  \item{mBins}{Number of bins for migration rate, if using grid search}
  \item{tauBins}{Number of bins for speciation times, if using grid search}
  \item{gridSearch}{Logical. Indicates whether a grid search will be used to optimize}
  \item{branchLengths}{Logical. Indicates whether branch lengths will be input}
  \item{pDiscord}{Minimum probability of discordance allowed. Defaults to .05}
  \item{method}{Optimization scheme passed to \code{\link{optim}}. Defaults to Nelder-Mead.}
  \item{startVals}{Initial values for \code{\link{optim}}. If NULL, the function will generate starting values uniformly at random within the parameter support.}
}
\details{
pDiscord is relevant when branchLengths = FALSE. This affects the upper bound on the migration rate in the parameter space.
}
\value{
\item{Matrix containing the maximum likelihood parameter estimates and negative log likelihood under each of the 6 possible gene tree orientations.}
}
\author{\email{hbk5086@psu.edu}}
\seealso{
\code{\link{topOptimST}}, \code{\link{gridST}}, \code{\link{blOptimST}}
}

