\name{topOptimST}
\alias{topOptimST}
\title{Construct the species tree using gene tree topologies and numerical optimization}
\usage{topOptimST(data, theta, method, startVals, N, pDiscord)}
\description{
This function takes a vector of counts of observed gene trees and constructs a species tree from them, while inferring the length of the internal branch as well as the migration rate between the 2 subpopulations. The user must also specify the population-scaled mutation rate theta, and the effective population size N.
}
\arguments{
    \item{data}{A numerical vector of length 3. Each vector corresponds to one of the three possible labelings for a rooted gene tree triplet. These rooted triples are ((AB)C), ((BC)A), and ((AC)B).}
    \item{theta}{Positive numerical value corresponding to the population-scaled mutation rate, $\theta=4N\mu$.}
    \item{method}{Optimization scheme passed to \link{optim}. Defaults to Nelder-Mead.}
    \item{startVals}{Numerical vector of length 2 of starting values for the optimization scheme. The first number in the vector is for the internal branch length, and the second is the migration rate.}
    \item{N}{Effective population size.}
    \item{pDiscord}{The minimum probability of gene tree discordance. Defaults to .05.}
}
\value{
    \item{Matrix containing the maximum likelihood parameter estimates and negative log likelihood under each of the 6 possible gene tree orientations.}
}
\note{If the migration rate is significantly large, the populations are effectively unstructured.}
\author{\email{hbk5086@psu.edu}}
\seealso{
\code{\link{buildST}}
}
\examples{
# Say you've observed each labeling of the gene tree rooted triple
# 50, 100, and 56 times. Then, you data vector would look like this:
data <- c(50,100,56)
names(data) <- c("((A,B),C);", "((B,C),A);", "((A,C),B);")

# Previous analyses may have told you that your taxa have an effective
# population size of N=10,000, and theta=.01
N <- 10000
theta <- .01

# Then, you can infer the maximum likelihood species tree
\dontrun{
out <- topOptimST(data=data, theta=theta, method="Nelder-Mead",
                    startVals=NULL, N=N, pDiscord=.05)
out
}
}
