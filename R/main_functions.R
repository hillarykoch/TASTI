 buildST <- function(data, theta, N, mBins = 20, tauBins = 20,
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
        gridST(data, N, theta, mBins, tauBins)
    } else if(branchLengths){
        blOptimST(data, N, theta, method, startVals)
    } else{
        topOptimST(data, theta, method, startVals, N, pDiscord)
    }
}

gridST <- function(data, N, theta, mBins, tauBins){
    colnames(data)[1] %>%
        gsub(pattern = "\\(", replacement = "") %>%
        gsub(pattern = "\\)", replacement = "") -> labs
    spl <- str_split(labs, "")[[1]]
    commaIdx <- which(spl == ",")
    A <- paste0(spl[1:(commaIdx[1]-1)], collapse = "")
    B <- paste0(spl[(commaIdx[1]+1):(commaIdx[2]-1)], collapse = "")
    C <- paste0(spl[(commaIdx[2]+1):(length(spl)-1)], collapse = "")

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

    colnames(nll) <- c(paste0("((", A, ",", B, "),", C, ");"),
                       paste0("((", B, ",", C, "),", A, ");"),
                       paste0("((", A, ",", C, "),", B, ");"),
                       paste0("((", B, ",", A, "),", C, ");"),
                       paste0("((", C, ",", B, "),", A, ");"),
                       paste0("((", C, ",", A, "),", B, ");"))

    for(l in 1:nrow(lookup)){
        nll[l,] <- coarseFunc(par=lookup[l,1:3], data=data, cABC)
    }
    mlEsts <- data.frame(lookup, nll)
    mle <- extractML(mlEsts, rowNames = colnames(nll))
    mle
}

topOptimST <- function(data, theta, method, startVals, N, pDiscord){
    names(data)[1] %>%
        gsub(pattern = "\\(", replacement = "") %>%
        gsub(pattern = "\\)", replacement = "") -> labs
    spl <- str_split(labs, "")[[1]]
    commaIdx <- which(spl == ",")
    A <- paste0(spl[1:(commaIdx[1]-1)], collapse = "")
    B <- paste0(spl[(commaIdx[1]+1):(commaIdx[2]-1)], collapse = "")
    C <- paste0(spl[(commaIdx[2]+1):(length(spl)-1)], collapse = "")

    cABC <- 2/theta
    tauBound <- -theta/2*log(3/2*pDiscord)

    if(is.null(startVals)){
        startVals <- c(runif(1, min=-10,max=log(tauBound)),
                       runif(1, min=-10, max=log(4*N/theta)))
    }

    mle <- matrix(rep(NA, 6*3), ncol = 3)
    colnames(mle) <- c("tau", "m", "nll")
    rownames(mle) <- c(paste0("((", A, ",", B, "),", C, ");"),
                        paste0("((", B, ",", C, "),", A, ");"),
                        paste0("((", A, ",", C, "),", B, ");"),
                        paste0("((", B, ",", A, "),", C, ");"),
                        paste0("((", C, ",", B, "),", A, ");"),
                        paste0("((", C, ",", A, "),", B, ");"))

    for(l in 1:6){
        opt <- optim(par = startVals, fn = topFunc,
                     method = method, data = data, theta = theta,
                     N = N, tauBound = tauBound, cABC = cABC, l=l)
        mle[l,] <- c(exp(opt$par), opt$value)
    }
    mle
}

blOptimST <- function(data, N, theta, method, startVals){
    colnames(data)[1] %>%
        gsub(pattern = "\\(", replacement = "") %>%
        gsub(pattern = "\\)", replacement = "") -> labs
    spl <- str_split(labs, "")[[1]]
    commaIdx <- which(spl == ",")
    A <- paste0(spl[1:(commaIdx[1]-1)], collapse = "")
    B <- paste0(spl[(commaIdx[1]+1):(commaIdx[2]-1)], collapse = "")
    C <- paste0(spl[(commaIdx[2]+1):(length(spl)-1)], collapse = "")

    cABC <- 2/theta
    tau1bound <- min(data[,4])
    tau2bound <- min(data[,5])

    mle <- matrix(rep(NA, 6*4), ncol = 4)
    colnames(mle) <- c("tau1", "tau2", "m", "nll")
    rownames(mle) <- c(paste0("((", A, ",", B, "),", C, ");"),
                       paste0("((", B, ",", C, "),", A, ");"),
                       paste0("((", A, ",", C, "),", B, ");"),
                       paste0("((", B, ",", A, "),", C, ");"),
                       paste0("((", C, ",", B, "),", A, ");"),
                       paste0("((", C, ",", A, "),", B, ");"))

    for(l in 1:6){
        if(tau1bound == 0 & tau2bound == 0){
            # If there is a "star tree"
            mle[l,] <- c(0,0,NA,-Inf)
        } else if(tau1bound == 0){
            # Adjust for tau1bound=0 by fixing tau1 and optimizing over tau2 and m
            if(is.null(startVals)){
                startVals <- c(runif(1, min=-10, max = log(tau2bound)), runif(1, min=-10, max=log(4*N/theta)))
            }
            opt <- optim(par = startVals, fn = blOptimFuncAdj, data = data,
                         method = method, cABC = cABC,
                         tau2bound = tau2bound, N=N, theta=theta, l=l)
            mle[l,] <- c(0, exp(opt$par), opt$value)
        } else{
            if(is.null(startVals)){
                startVals <- c(runif(1, min = -10, max = log(tau1bound)),
                               runif(1, min = -10, max = log(tau2bound)),
                               runif(1, min=-10, max=log(4*N/theta)))
            }

            if(startVals[2] < startVals[1]){
                startVals <- startVals[c(2,1,3)]
            }

            opt <- optim(par = startVals, fn = blOptimFunc, data = data,
                         method = method, cABC=cABC,
                         tau1bound=tau1bound, tau2bound=tau2bound,
                         N=N, theta=theta, l=l)

            mle[l,] <- c(exp(opt$par), opt$value)
        }
    }
    mle
}
