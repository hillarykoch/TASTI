\name{gridST}
\alias{gridST}
\title{
Construct the species tree using gene tree topologies and branch lengths with a grid search
}
\description{
This function takes a matrix of counts of observed gene trees with branch lengths and constructs a species tree from them, while inferring the speciation times and the migration rate between the 2 subpopulations. The user must also specify the population-scaled mutation rate theta, and the precision of the grid search with mBins and tauBins.
}
\usage{
gridST(data, N, theta, mBins, tauBins)
}
\arguments{
  \item{data}{A matrix with 5 columns and number of rows equal to the number of gene trees. Each row of the matrix corresponds to a given gene tree. The first three entries in the rows are either 0 or 1 -- 1 corresponds to the given gene tree observed. The column labelings here are the three possible labelings for a rooted gene tree triplet. These rooted triples are ((AB)C), ((BC)A), and ((AC)B). The 4th entry in the row is t1, the most recent coalescent event among the 3 taxa, going back in time. The 5th entry is t2, the second coalescent event. See Example.}
  \item{N}{Effective population size.}
  \item{theta}{Population-scaled mutation rate.}
  \item{mBins}{Number of bins to divide the parameter space for the migration rate into, for grid search. The larger the value, the better the precision. However, larger values come with greater computational burden.}
  \item{tauBins}{Number of bins to divide the parameter space for the speciation times into, for grid search. The larger the value, the better the precision. However, larger values come with greater computational burden.}
  }
\value{
Matrix containing the maximum likelihood parameter estimates and negative log likelihood under each of the 6 possible gene tree orientations.
}
\author{\email{hbk5086@psu.edu}}
\seealso{\code{\link{buildST}}}
\examples{
# Consider you have 3 gene trees, one with these topologies and branch lengths:
# ((A:0.01, B:0.01):0.002, C:0.012)
# ((B:0.001, C:0.001):0.0025, A:0.0035)
# ((A:0.02, C:0.02):0.01, B:0.03)
# These can be encoded in the following matrix
data <- matrix(c(1,0,0,.01,.012,
                 0,1,0,.001,.0035,
                 0,0,1,0.02,0.03), ncol = 5)

# Let's say you also have an effective population size N and population-scaled
# mutation rate theta:
theta <- .005
N <- 10000

# For computational efficiency, we will choose a very small number of bins for
# the grid search (an overly coarse grid). In reality, one may want to choose
# ~20 bins or more each for migration rate and divergence times.

mle <- gridST(data=data, N=N, theta=theta, mBins=4, tauBins=4)
mle
}
