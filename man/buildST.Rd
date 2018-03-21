\name{buildST}
\alias{buildST}
\title{
Build a species tree from gene tree data
}
\description{
Construct the most likely species tree from gene trees with or without branch lengths, using either \link{optim} or a grid search.
}
\usage{
buildST(data, theta, N, mBins = 20, tauBins = 20,
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
\examples{
# Consider you have 3 gene trees, one with these topologies and branch lengths:
# ((A:0.01, B:0.01):0.002, C:0.012)
# ((B:0.001, C:0.001):0.0025, A:0.0035)
# ((A:0.02, C:0.02):0.01, B:0.03)
# These can be encoded in the following matrix
data <- matrix(c(1,0,0,.01,.012,
                 0,1,0,.001,.0035,
                 0,0,1,0.02,0.03), ncol = 5, byrow = T)
colnames(data) <- c("((A,B),C);", "((B,C),A);", "((A,C),B);", "t1", "t2")

# Let's say you also have an effective population size N and population-scaled
# mutation rate theta:
theta <- .005
N <- 10000

# Builds using a grid search with specified bin numbers
mle1 <- buildST(data, theta, N, mBins = 4, tauBins = 4)

# Builds using numerical optimization
mle2 <- buildST(data, theta, N, gridSearch = FALSE)

# Builds using only topologies, no branch lengths
topo_data <- colSums(data)[1:3]
\dontrun{
mle3 <- topOptimST(data, theta, N, gridSearch = FALSE, branchLengths = FALSE)
}

mle1
mle2
mle3
}

