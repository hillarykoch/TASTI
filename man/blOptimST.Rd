\name{blOptimST}
\alias{blOptimST}
\title{
Construct the species tree using gene tree topologies and branch lengths with numerical optimization
}
\description{
This function takes a matrix of counts of observed gene trees with branch lengths and constructs a species tree from them, while inferring the speciation times and the migration rate between the 2 subpopulations. The user must specify the population-scaled mutation rate theta, and the effective population size N.
}
\usage{
blOptimST(data, N, theta, method, startVals)
}
\arguments{
  \item{data}{A matrix with 5 columns and number of rows equal to the number of gene trees. Each row of the matrix corresponds to a given gene tree. The first three entries in the rows are either 0 or 1 -- 1 corresponds to the given gene tree observed. The column labelings here are the three possible labelings for a rooted gene tree triplet. These rooted triples are ((AB)C), ((BC)A), and ((AC)B). The 4th entry in the row is t1, the most recent coalescent event among the 3 taxa, going back in time. The 5th entry is t2, the second coalescent event. See Example.}
  \item{N}{Effective population size.}
  \item{theta}{Population-scaled mutation rate.}
  \item{method}{Optimization scheme passed to \link{optim}. Defaults to Nelder-Mead.}
  \item{startVals}{Numerical vector of length 2 of starting values for the optimization scheme. The first number in the vector is for the internal branch length, and the second is the migration rate.}
}
\value{
Matrix containing the maximum likelihood parameter estimates and negative log likelihood under each of the 6 possible gene tree orientations.
}
\author{\email{hbk5086@psu.edu}}
\note{
If the migration rate is significantly large, the populations are effectively unstructured.
}

\seealso{\code{\link{buildST}}}
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

mle <- blOptimST(data=data, N=N, theta=theta, method="Nelder-Mead", startVals = NULL)
mle
}
