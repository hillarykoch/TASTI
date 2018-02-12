# # This is how you save a newick tree as a variable without writing to a file
# nwk <- write.tree(phylip@phylo, file = "")
# phyloApe <- phylip@phylo
#
# # Get all rooted triples
# labs <- phylip@phylo$tip.label
# triples <- combn(labs, 3)
#
# # Unrooted triple will have 1 node
# # Rooted triple will have 2 nodes
# outTip <- setdiff(labs, triples[,100])
# rooted <- drop.tip(phyloApe, outTip)
# unr <- drop.tip(phyloApe, setdiff(labs, triples[,1]))

buildST <- function(data, m, theta, N = NULL, mBins = 20, tauBins = 20,
                    gridSearch = TRUE, branchLengths = TRUE, pDiscord = .05,
                    method = "Nelder-Mead", startVals = NULL){
    if(gridSearch & !branchLengths){
        stop("grid search optimization is only available when using
             branch lengths.")
    }

    if(gridSearch & !is.null(startVals)){
        message("start values are not used in grid search.")
    }

    if(!gridSearch & !branchLengths & is.null(N)){
        stop("Effective population size (N) must be specified for inference
             with topologies")
    }

    if(gridSearch){
        gridST(data, m, theta, mBins, tauBins)
    } else if(branchLengths){
        blOptimST(data, m, theta, method, startVals)
    } else{
        topOptimST(data, theta, method, startVals, N, pDiscord)
    }
}

gridST <- function(data, N, theta, mBins, tauBins){
    cABC <- 2/theta
    upperM <- log(4*N/theta)

    tau1bound <- min(data[,4])
    tau2bound <- min(data[,5])

    r <- c(0, tau2bound)
    tbins <- seq(r[1], r[2], length.out = tauBins+1)
    tBinLen <- diff(tbins[1:2])
    tpairs <- combn(tail(tbins,tauBins),2) %>%
        cbind(matrix(rep(tail(tbins,tauBins), 2), nrow = 2, byrow = T))
    mbins <- 10^(seq(-7,upperM,length.out=mBins+1))

    lookup <- matrix(rep(NA, ncol(tpairs)*mBins*3), ncol = 3)
    lookup[,1] <- rep(tpairs[1,], mBins)-tBinLen/2
    lookup[,2] <- rep(tpairs[2,], mBins)-tBinLen/2
    lookup[,3] <- rep(tail(mbins, mBins)-.5*diff(mbins), each = ncol(tpairs))

    colnames(lookup) <- c("tau1bin", "tau2bin", "mbin")
    lookup <- data.frame(lookup)

    nll <- matrix(rep(NA, 6*nrow(lookup)), ncol = 6)
    colnames(nll) <- c("ABC", "BCA", "ACB", "BAC", "CBA", "CAB")
    for(l in 1:nrow(lookup)){
        nll[l,] <- coarseFunc(par=lookup[l,1:3], data=data, cABC)
    }
    mlEsts <- data.frame(lookup, nll)
    mle <- extractML(mlEsts)
    mle
}

topOptimST <- function(data, theta, method, startVals, N, pDiscord){
    cABC <- 2/theta
    tauBound <- -theta/2*log(3/2*pDiscord)

    if(is.null(startVals)){
        startVals <- c(runif(1, min=-7,max=log(tauBound)),
                       runif(1, min=-7, max=log(4*N/theta)))
    }

    mle <- matrix(rep(NA, 6*3), ncol = 3)
    colnames(mle) <- c("tau", "m", "nll")
    rownames(mle) <- c("ABC", "BCA", "ACB", "BAC", "CBA", "CAB")

    for(l in 1:6){
        opt <- optim(par = startVals, fn = topFunc,
                     method = method, data = data, theta = theta,
                     N = N, tauBound = tauBound, cABC = cABC, l=l)
        mle[l,] <- c(exp(opt$par), opt$value)
    }
}
