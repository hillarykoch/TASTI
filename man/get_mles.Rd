\name{get_mles}
\alias{get_mles}
\title{
Get maximum likelihood parameter estimates over a list of model fits
}
\description{
If several model fits are stored in a list, use get_mles to get the maximum likelihood parameter estimates from each fit. It also will count how many times each species tree topology and configuration was the maximum likelihood species tree given the data.
}
\usage{
get_mles(fits)
}
\arguments{
  \item{fits}{A list of model fits produces using \code{\link{buildST}}}
}
\value{
A list of maximum likelihood parameter estimates corresponding to each element in the list of model fits. Additionally, one vector summarizing how often each species tree topology and configuration was the maximum likelihood specis tree, given the data, across all input model fits.
}
\author{\email{hbk5086@psu.edu}}
\seealso{\code{\link{buildST}}}
\examples{
# Suppose you have 3 data sets, each stored as an element in a list called "fits"
theta <- .005
N <- 10000
data <- list("rep1" = matrix(c(1,0,0,.01,.012,
                                0,1,0,.001,.0035,
                                0,0,1,0.02,0.03), ncol = 5, byrow = TRUE),
              "rep2" = matrix(c(1,0,0,.1,.2,
                                1,0,0,.001,.0048,
                                0,0,1,0.02,0.1,
                                0,1,0,.01,.01), ncol = 5, byrow = TRUE),
              "rep3" = matrix(c(0,1,0,.01,.15,
                                0,1,0,.11,.3,
                                0,0,1,0.02,0.04), ncol = 5, byrow = TRUE))
colnames(data) <- c("((A,B),C);", "((B,C),A);", "((A,C),B);", "t1", "t2")

# Fit a model to each replicate. Let's consider using branch lengths with a grid search
fits <- lapply(data, function(X) gridST(X, N, theta, mBins = 5, tauBins = 5))

# Extract the info!
mles <- get_mles(fits)
mles
}
