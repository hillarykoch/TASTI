## These two general functions should be fine for getting all the necessary integrals
## coal1 is for the first coalescence above the root
coal1 <- function(T1, state1, state2, rate){
    rate*T1[state1,state2]
}

## Function outputs probability of coalescence along internal branch
f1Grid <- function(t1, t2, tau1, tau2, QAB, QABC, cABC) {
    if(tau2 < t1 | t2 < tau2 | t1 < tau1){
        return(0)
    } else{
        t1t1 <- expm(QAB*(t1-tau1))
        t2t1 <- expm(QAB*(tau2-t1))
        t1t2 <- expm(QABC*(t2-tau2))

        c.ab <-
            cABC*t1t1[1,3]*(t2t1[5,5]*(cABC*t1t2[10,9] + cABC*t1t2[10,12]) +
                            t2t1[5,6]*(cABC*t1t2[12,9] + cABC*t1t2[12,12])) +
            cABC*t1t1[1,4]*(t2t1[6,5]*(cABC*t1t2[10,9] + cABC*t1t2[10,12]) +
                            t2t1[6,6]*(cABC*t1t2[12,9] + cABC*t1t2[12,12]))
        if(c.ab < 0){
            c.ab <- 0
            message("negative values turned to 0.")
        }
        return(c.ab)
    }
}

## Function outputs probability of AB coalescing first, above the root
f2Grid <- function(t1, t2, tau1, tau2, QAB, QABC, cABC) {
    if(tau2 > t1 | t2 < tau2 | t1 < tau1){
        return(0)
    } else {
        val <- expm(QAB*(tau2-tau1))[1,1:4]
        norm.c <- sum(val)

        alpha1 <- val[1]/norm.c
        alpha2 <- val[2]/norm.c
        alpha3 <- val[3]/norm.c
        alpha4 <- val[4]/norm.c

        t1t2 <- expm(QABC*(t1-tau2))
        t2t1 <- expm(QABC*(t2-t1))

        c.ab <-
            alpha1*(coal1(t1t2,5,1,3*cABC)*(coal1(t2t1,9,9,cABC) +
                    coal1(t2t1,9,12,cABC)) +
                    coal1(t1t2,5,2,cABC)*(coal1(t2t1,10,9,cABC) +
                    coal1(t2t1,10,12,cABC)) +
                    coal1(t1t2,5,7,cABC)*(coal1(t2t1,11,9,cABC) +
                    coal1(t2t1,11,12,cABC)) +
                    coal1(t1t2,5,8,3*cABC)*(coal1(t2t1,12,9,cABC) +
                    coal1(t2t1,12,12,cABC))) +
            alpha2*(coal1(t1t2,6,1,3*cABC)*(coal1(t2t1,9,9,cABC) +
                    coal1(t2t1,9,12,cABC)) +
                    coal1(t1t2,6,2,cABC)*(coal1(t2t1,10,9,cABC) +
                    coal1(t2t1,10,12,cABC)) +
                    coal1(t1t2,6,7,cABC)*(coal1(t2t1,11,9,cABC) +
                    coal1(t2t1,11,12,cABC)) +
                    coal1(t1t2,6,8,3*cABC)*(coal1(t2t1,12,9,cABC) +
                    coal1(t2t1,12,12,cABC))) +
            alpha3*(coal1(t1t2,2,1,3*cABC)*(coal1(t2t1,9,9,cABC) +
                    coal1(t2t1,9,12,cABC)) +
                    coal1(t1t2,2,2,cABC)*(coal1(t2t1,10,9,cABC) +
                    coal1(t2t1,10,12,cABC)) +
                    coal1(t1t2,2,7,cABC)*(coal1(t2t1,11,9,cABC) +
                    coal1(t2t1,11,12,cABC)) +
                    coal1(t1t2,2,8,3*cABC)*(coal1(t2t1,12,9,cABC) +
                    coal1(t2t1,12,12,cABC))) +
            alpha4*(coal1(t1t2,8,1,3*cABC)*(coal1(t2t1,9,9,cABC) +
                    coal1(t2t1,9,12,cABC)) +
                    coal1(t1t2,8,2,cABC)*(coal1(t2t1,10,9,cABC) +
                    coal1(t2t1,10,12,cABC)) +
                    coal1(t1t2,8,7,cABC)*(coal1(t2t1,11,9,cABC) +
                    coal1(t2t1,11,12,cABC)) +
                    coal1(t1t2,8,8,3*cABC)*(coal1(t2t1,12,9,cABC) +
                    coal1(t2t1,12,12,cABC)))

        if(c.ab < 0){
            c.ab <- 0
            message("negative values turned to 0.")
        }

        return(c.ab)
    }
}

## Function outputs probability of BC coalescing first, above the root
f3Grid <- function(t1, t2, tau1, tau2, QAB, QABC, cABC) {
    if(tau2 > t1 | t2 < tau2 | t1 < tau1){
        return(0)
    } else {
        val <- expm(QAB*(tau2-tau1))[1,1:4]
        norm.c = sum(val)

        alpha1 <- val[1]/norm.c
        alpha2 <- val[2]/norm.c
        alpha3 <- val[3]/norm.c
        alpha4 <- val[4]/norm.c

        t1t2 <- expm(QABC*(t1-tau2))
        t2t1 <- expm(QABC*(t2-t1))

        c.bc <-
            alpha1*(coal1(t1t2,5,1,3*cABC)*(coal1(t2t1,17,17,cABC) +
                    coal1(t2t1,17,20,cABC)) +
                    coal1(t1t2,5,4,cABC)*(coal1(t2t1,18,17,cABC) +
                    coal1(t2t1,18,20,cABC)) +
                    coal1(t1t2,5,5,cABC)*(coal1(t2t1,19,17,cABC) +
                    coal1(t2t1,19,20,cABC)) +
                    coal1(t1t2,5,8,3*cABC)*(coal1(t2t1,20,17,cABC) +
                    coal1(t2t1,20,20,cABC))) +
            alpha2*(coal1(t1t2,6,1,3*cABC)*(coal1(t2t1,17,17,cABC) +
                    coal1(t2t1,17,20,cABC)) +
                    coal1(t1t2,6,4,cABC)*(coal1(t2t1,18,17,cABC) +
                    coal1(t2t1,18,20,cABC)) +
                    coal1(t1t2,6,5,cABC)*(coal1(t2t1,19,17,cABC) +
                    coal1(t2t1,19,20,cABC)) +
                    coal1(t1t2,6,8,3*cABC)*(coal1(t2t1,20,17,cABC) +
                    coal1(t2t1,20,20,cABC))) +
            alpha3*(coal1(t1t2,2,1,3*cABC)*(coal1(t2t1,17,17,cABC) +
                    coal1(t2t1,17,20,cABC)) +
                    coal1(t1t2,2,4,cABC)*(coal1(t2t1,18,17,cABC) +
                    coal1(t2t1,18,20,cABC)) +
                    coal1(t1t2,2,5,cABC)*(coal1(t2t1,19,17,cABC) +
                    coal1(t2t1,19,20,cABC)) +
                    coal1(t1t2,2,8,3*cABC)*(coal1(t2t1,20,17,cABC) +
                    coal1(t2t1,20,20,cABC))) +
            alpha4*(coal1(t1t2,8,1,3*cABC)*(coal1(t2t1,17,17,cABC) +
                    coal1(t2t1,17,20,cABC)) +
                    coal1(t1t2,8,4,cABC)*(coal1(t2t1,18,17,cABC) +
                    coal1(t2t1,18,20,cABC)) +
                    coal1(t1t2,8,5,cABC)*(coal1(t2t1,19,17,cABC) +
                    coal1(t2t1,19,20,cABC)) +
                    coal1(t1t2,8,8,3*cABC)*(coal1(t2t1,20,17,cABC) +
                    coal1(t2t1,20,20,cABC)))

        if(c.bc < 0){
            c.bc <- 0
            message("negative values turned to 0.")
        }

        return(c.bc)
    }
}

## Function outputs probability of AC coalescing first, above the root
f4Grid <- function(t1, t2, tau1, tau2, QAB, QABC, cABC) {
    if(tau2 > t1 | t2 < tau2 | t1 < tau1){
        return(0)
    } else{
        val <- expm(QAB*(tau2-tau1))[1,1:4]
        norm.c <- sum(val)

        alpha1 <- val[1]/norm.c
        alpha2 <- val[2]/norm.c
        alpha3 <- val[3]/norm.c
        alpha4 <- val[4]/norm.c

        t1t2 <- expm(QABC*(t1-tau2))
        t2t1 <- expm(QABC*(t2-t1))

        c.ac <-
            alpha1*(coal1(t1t2,5,1,3*cABC)*(coal1(t2t1,13,13,cABC) +
                    coal1(t2t1,13,16,cABC)) +
                    coal1(t1t2,5,3,cABC)*(coal1(t2t1,15,13,cABC) +
                    coal1(t2t1,15,16,cABC)) +
                    coal1(t1t2,5,6,cABC)*(coal1(t2t1,14,13,cABC) +
                    coal1(t2t1,14,16,cABC)) +
                    coal1(t1t2,5,8,3*cABC)*(coal1(t2t1,16,13,cABC) +
                    coal1(t2t1,16,16,cABC))) +
            alpha2*(coal1(t1t2,6,1,3*cABC)*(coal1(t2t1,13,13,cABC) +
                    coal1(t2t1,13,16,cABC)) +
                    coal1(t1t2,6,3,cABC)*(coal1(t2t1,15,13,cABC) +
                    coal1(t2t1,15,16,cABC)) +
                    coal1(t1t2,6,6,cABC)*(coal1(t2t1,14,13,cABC) +
                    coal1(t2t1,14,16,cABC)) +
                    coal1(t1t2,6,8,3*cABC)*(coal1(t2t1,16,13,cABC) +
                    coal1(t2t1,16,16,cABC))) +
            alpha3*(coal1(t1t2,2,1,3*cABC)*(coal1(t2t1,13,13,cABC) +
                    coal1(t2t1,13,16,cABC)) +
                    coal1(t1t2,2,3,cABC)*(coal1(t2t1,15,13,cABC) +
                    coal1(t2t1,15,16,cABC)) +
                    coal1(t1t2,2,6,cABC)*(coal1(t2t1,14,13,cABC) +
                    coal1(t2t1,14,16,cABC)) +
                    coal1(t1t2,2,8,3*cABC)*(coal1(t2t1,16,13,cABC) +
                    coal1(t2t1,16,16,cABC))) +
            alpha4*(coal1(t1t2,8,1,3*cABC)*(coal1(t2t1,13,13,cABC) +
                    coal1(t2t1,13,16,cABC)) +
                    coal1(t1t2,8,3,cABC)*(coal1(t2t1,15,13,cABC) +
                    coal1(t2t1,15,16,cABC)) +
                    coal1(t1t2,8,6,cABC)*(coal1(t2t1,14,13,cABC) +
                    coal1(t2t1,14,16,cABC)) +
                    coal1(t1t2,8,8,3*cABC)*(coal1(t2t1,16,13,cABC) +
                    coal1(t2t1,16,16,cABC)))

        if(c.ac < 0){
            c.ac <- 0
            message("negative values turned to 0.")
        }
        return(c.ac)
    }
}


# Function to compute the likelihood given a t1, t2, m
coarseFunc <- function(par, data, cABC){
    tau1 <- as.numeric(par[1])
    tau2 <- as.numeric(par[2])
    m <- as.numeric(par[3])
    t1.vec <- data[,4]
    t2.vec <- data[,5]

    ABC.idx <- data[,1] == 1
    BCA.idx <- data[,2] == 1
    ACB.idx <- data[,3] == 1


    QAB <- matrix(data = c(-2*m, 0, m, m, 0, 0,
                           0, -2*m, m, m, 0, 0,
                           m, m, -2*m-cABC, 0, cABC, 0,
                           m, m, 0, -2*m-cABC, 0, cABC,
                           0, 0, 0, 0, -m, m,
                           0, 0, 0, 0, m, -m), nrow = 6, byrow = TRUE)


    QABC <- matrix(data = c(-3*m-3*(cABC + cABC + cABC), m, m, m, 0, 0, 0,
                            0, 3*cABC, 0, 0, 0, 3*cABC, 0, 0, 0, 3*cABC, 0, 0,
                            0, 0, 0, m,-3*m-cABC,0,0,m,m,0,0,0,cABC,0,0,0,0,0,
                            0,0,0,0,0,0,0,m,0,-3*m-cABC,0,m,0,m,0,0,0,0,0,0,
                            cABC,0,0,0,0,0,0,0,0,m,0,0,-3*m-cABC,0,m,m,0,0,0,0,
                            0,0,0,0,0,0,cABC,0,0,0,0,0,m,m,0,-3*m-cABC,0,0,m,0,
                            0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,m,0,m,0,-3*m-cABC,0,
                            m,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,m,m,0,0,
                            -3*m-cABC,m,0,0,cABC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            m,m,m,-3*m-3*(cABC + cABC + cABC),0,0,0,3*cABC,0,0,
                            0,3*cABC,0,0,0,3*cABC,0,0,0,0,0,0,0,0,0,0,-2*m-cABC,
                            m,m,0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,m,
                            -2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m,0,
                            -2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m,m,
                            -2*m-cABC,0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,
                            0,0,0,-2*m-cABC,m,m,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,m,-2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,m,0,-2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,m,m,-2*m-cABC,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,-2*m-cABC,m,m,0,cABC,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,m,-2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,m,0,-2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,m,m,-2*m-cABC,0,cABC,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0,-m,m,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,m,-m), nrow = 22, byrow = TRUE)

    QAB <- matrix(unlist(QAB), nrow = 6)
    QABC <- matrix(unlist(QABC), nrow = 22)
    epsilon <- 10^(-15)
    reps <- nrow(data)

    ll.raw <- matrix(rep(NA,3*reps), ncol = 3)
    for(i in 1:reps){
        ll.raw[i,1] <- log(f1Grid(t1.vec[i], t2.vec[i],
                                  tau1, tau2, QAB, QABC, cABC) +
                               f2Grid(t1.vec[i], t2.vec[i],
                                      tau1, tau2, QAB, QABC, cABC) + epsilon)
        ll.raw[i,2] <- log(f3Grid(t1.vec[i], t2.vec[i],
                              tau1, tau2, QAB, QABC, cABC) + epsilon)
        ll.raw[i,3] <- log(f4Grid(t1.vec[i], t2.vec[i],
                              tau1, tau2, QAB, QABC, cABC) + epsilon)
    }

    ll.ABC <- sum(ll.raw[ABC.idx,1]) +
        sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,3])
    ll.BCA <- sum(ll.raw[ABC.idx,3]) +
        sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,2])
    ll.ACB <- sum(ll.raw[ABC.idx,3]) +
        sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,1])
    ll.BAC <- sum(ll.raw[ABC.idx,1]) +
        sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,2])
    ll.CBA <- sum(ll.raw[ABC.idx,2]) +
        sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,3])
    ll.CAB <- sum(ll.raw[ABC.idx,2]) +
        sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,1])

    ll <- c(ll.ABC, ll.BCA, ll.ACB, ll.BAC, ll.CBA, ll.CAB)

    return(-ll)
}

extractML <- function(mlEsts, rowNames){
    optIdx <- -mlEsts[,4:9] %>% t %>% max.col
    mle <- cbind(mlEsts[optIdx,1:3], diag(as.matrix(mlEsts[optIdx,4:9])))
    colnames(mle) <- c(colnames(mle[1:3]), "nll")
    rownames(mle) <- rowNames
    mle
}

f1Top <- function(t1, cABC, QAB){
    return(c(cABC, cABC)%*%(expm(QAB*t1)[1,3:4]))
}

f2Top <- function(t1, tau, cABC, QAB, QABC){
    val <- expm(QAB*tau)[1,1:4]
    norm.c <- sum(val)

    alpha1 <- val[1]/norm.c
    alpha2 <- val[2]/norm.c
    alpha3 <- val[3]/norm.c
    alpha4 <- val[4]/norm.c

    mat <- expm(QABC*(t1-tau))

    c.ab <- alpha1*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[5,c(1:2,7:8)] +
        alpha2*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[6,c(1:2,7:8)] +
        alpha3*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[2,c(1:2,7:8)] +
        alpha4*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[8,c(1:2,7:8)]

    return(c.ab)
}

f3Top <- function(t1, tau, cABC, QAB, QABC){
    val <- expm(QAB*tau)[1,1:4]
    norm.c <- sum(val)

    alpha1 <- val[1]/norm.c
    alpha2 <- val[2]/norm.c
    alpha3 <- val[3]/norm.c
    alpha4 <- val[4]/norm.c

    mat <- expm(QABC*(t1-tau))

    c.bc <- alpha1*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[5, c(1,4,5,8)] +
        alpha2*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[6, c(1,4,5,8)] +
        alpha3*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[2, c(1,4,5,8)] +
        alpha4*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[8, c(1,4,5,8)]

    return(c.bc)
}

f4Top <- function(t1, tau, cABC, QAB, QABC){
    val <- expm(QAB*tau)[1,1:4]
    norm.c <- sum(val)

    alpha1 <- val[1]/norm.c
    alpha2 <- val[2]/norm.c
    alpha3 <- val[3]/norm.c
    alpha4 <- val[4]/norm.c

    mat <- expm(QABC*(t1-tau))

    c.ac <- alpha1*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[5, c(1,3,6,8)] +
        alpha2*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[6, c(1,3,6,8)] +
        alpha3*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[2, c(1,3,6,8)] +
        alpha4*c(3*cABC, cABC, cABC, 3*cABC)%*%mat[8, c(1,3,6,8)]

    return(c.ac)
}

topFunc <- function(par, method, data, theta, N, tauBound, cABC, l){
    if(!(-10 <= par[1] & par[1] <= log(tauBound))){
        return(Inf)
    }

    if(!(-10 <= par[2] & par[2] <= 4*N/theta)){
        return(Inf)
    }

    tau <- exp(par[1])
    m <- exp(par[2])

    QAB <- matrix(data = c(-2*m, 0, m, m, 0, 0,
                           0, -2*m, m, m, 0, 0,
                           m, m, -2*m-cABC, 0, cABC, 0,
                           m, m, 0, -2*m-cABC, 0, cABC,
                           0, 0, 0, 0, -m, m,
                           0, 0, 0, 0, m, -m), nrow = 6, byrow = TRUE)

    QABC <- matrix(data = c(-3*m-3*(cABC + cABC + cABC), m, m, m, 0, 0, 0,
                            0, 3*cABC, 0, 0, 0, 3*cABC, 0, 0, 0, 3*cABC, 0, 0,
                            0, 0, 0, m,-3*m-cABC,0,0,m,m,0,0,0,cABC,0,0,0,0,0,
                            0,0,0,0,0,0,0,m,0,-3*m-cABC,0,m,0,m,0,0,0,0,0,0,
                            cABC,0,0,0,0,0,0,0,0,m,0,0,-3*m-cABC,0,m,m,0,0,0,0,
                            0,0,0,0,0,0,cABC,0,0,0,0,0,m,m,0,-3*m-cABC,0,0,m,0,
                            0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,m,0,m,0,-3*m-cABC,0,
                            m,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,m,m,0,0,
                            -3*m-cABC,m,0,0,cABC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            m,m,m,-3*m-3*(cABC + cABC + cABC),0,0,0,3*cABC,0,0,
                            0,3*cABC,0,0,0,3*cABC,0,0,0,0,0,0,0,0,0,0,-2*m-cABC,
                            m,m,0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,m,
                            -2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m,0,
                            -2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m,m,
                            -2*m-cABC,0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,
                            0,0,0,-2*m-cABC,m,m,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,m,-2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,m,0,-2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,m,m,-2*m-cABC,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,-2*m-cABC,m,m,0,cABC,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,m,-2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,m,0,-2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,m,m,-2*m-cABC,0,cABC,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0,-m,m,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,m,-m), nrow = 22, byrow = TRUE)

    gf1.est <- adaptIntegrate(f1Top, lowerLimit = 0, upperLimit = tau,
                                 cABC, QAB)$integral
    gf2.est <- adaptIntegrate(f2Top, lowerLimit = tau, upperLimit = 100,
                                 tau, cABC, QAB, QABC)$integral*(1-gf1.est)
    gf3.est <- adaptIntegrate(f3Top, lowerLimit = tau, upperLimit = 100,
                                 tau, cABC, QAB, QABC)$integral*(1-gf1.est)
    gf4.est <- adaptIntegrate(f4Top, lowerLimit = tau, upperLimit = 100,
                                 tau, cABC, QAB, QABC)$integral*(1-gf1.est)
    p.est <- c(gf1.est + gf2.est, gf3.est, gf4.est)

    n1 <- data[1]
    n2 <- data[2]
    n3 <- data[3]

    epsilon <- 10^(-15)
    ll.ABC <- n1*log(p.est[1]+epsilon) +
        n2*log(p.est[2]+epsilon) + n3*log(p.est[3]+epsilon)
    ll.BCA <- n1*log(p.est[3]+epsilon) +
        n2*log(p.est[1]+epsilon) + n3*log(p.est[2]+epsilon)
    ll.ACB <- n1*log(p.est[3]+epsilon) +
        n2*log(p.est[2]+epsilon) + n3*log(p.est[1]+epsilon)
    ll.BAC <- n1*log(p.est[1]+epsilon) +
        n2*log(p.est[3]+epsilon) + n3*log(p.est[2]+epsilon)
    ll.CBA <- n1*log(p.est[2]+epsilon) +
        n2*log(p.est[1]+epsilon) + n3*log(p.est[3]+epsilon)
    ll.CAB <- n1*log(p.est[2]+epsilon) +
        n2*log(p.est[3]+epsilon) + n3*log(p.est[1]+epsilon)
    ll <- c(ll.ABC, ll.BCA, ll.ACB, ll.BAC, ll.CBA, ll.CAB)

    return(-ll[l])
}

blOptimFunc <- function(par, data, method, cABC, tau1bound, tau2bound, N, theta, l){
        par <- exp(par)
        if(!(par[1] >= 0 & par[1] <= tau1bound)){
            return(Inf)
        }

        if(!(par[2] >=0 & par[2] <= tau2bound)){
            return(Inf)
        }

        if(!(par[3] >= 0 & par[3] <= 4*N/theta)){
            return(Inf)
        }

        tau1 <- min(par[1:2])
        tau2 <- max(par[1:2])
        m <- par[3]
        t1.vec <- data[,4]
        t2.vec <- data[,5]

        ABC.idx <- data[,1] == 1
        BCA.idx <- data[,2] == 1
        ACB.idx <- data[,3] == 1

        QAB <- matrix(data = c(-2*m, 0, m, m, 0, 0,
                               0, -2*m, m, m, 0, 0,
                               m, m, -2*m-cABC, 0, cABC, 0,
                               m, m, 0, -2*m-cABC, 0, cABC,
                               0, 0, 0, 0, -m, m,
                               0, 0, 0, 0, m, -m), nrow = 6, byrow = TRUE)

        QABC <- matrix(data = c(-3*m-3*(cABC + cABC + cABC), m, m, m, 0, 0, 0,
                                0, 3*cABC, 0, 0, 0, 3*cABC, 0, 0, 0, 3*cABC, 0, 0,
                                0, 0, 0, m,-3*m-cABC,0,0,m,m,0,0,0,cABC,0,0,0,0,0,
                                0,0,0,0,0,0,0,m,0,-3*m-cABC,0,m,0,m,0,0,0,0,0,0,
                                cABC,0,0,0,0,0,0,0,0,m,0,0,-3*m-cABC,0,m,m,0,0,0,0,
                                0,0,0,0,0,0,cABC,0,0,0,0,0,m,m,0,-3*m-cABC,0,0,m,0,
                                0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,m,0,m,0,-3*m-cABC,0,
                                m,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,m,m,0,0,
                                -3*m-cABC,m,0,0,cABC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                m,m,m,-3*m-3*(cABC + cABC + cABC),0,0,0,3*cABC,0,0,
                                0,3*cABC,0,0,0,3*cABC,0,0,0,0,0,0,0,0,0,0,-2*m-cABC,
                                m,m,0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,m,
                                -2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m,0,
                                -2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m,m,
                                -2*m-cABC,0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,
                                0,0,0,-2*m-cABC,m,m,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,m,-2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,m,0,-2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,m,m,-2*m-cABC,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,-2*m-cABC,m,m,0,cABC,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,m,-2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,m,0,-2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,m,m,-2*m-cABC,0,cABC,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,-m,m,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,m,-m), nrow = 22, byrow = TRUE)

        QAB <- matrix(unlist(QAB), nrow = 6)
        QABC <- matrix(unlist(QABC), nrow = 22)
        epsilon <- 10^(-15)
        reps <- nrow(data)

        ll.raw <- matrix(rep(NA,3*reps), ncol = 3)
        for(i in 1:reps){
            ll.raw[i,1] <- log(f1Grid(t1.vec[i], t2.vec[i],
                                      tau1, tau2, QAB, QABC, cABC) +
                                   f2Grid(t1.vec[i], t2.vec[i],
                                          tau1, tau2, QAB, QABC, cABC) + epsilon)
            ll.raw[i,2] <- log(f3Grid(t1.vec[i], t2.vec[i],
                                      tau1, tau2, QAB, QABC, cABC) + epsilon)
            ll.raw[i,3] <- log(f4Grid(t1.vec[i], t2.vec[i],
                                      tau1, tau2, QAB, QABC, cABC) + epsilon)
        }

        ll.ABC <- sum(ll.raw[ABC.idx,1]) +
            sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,3])
        ll.BCA <- sum(ll.raw[ABC.idx,3]) +
            sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,2])
        ll.ACB <- sum(ll.raw[ABC.idx,3]) +
            sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,1])
        ll.BAC <- sum(ll.raw[ABC.idx,1]) +
            sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,2])
        ll.CBA <- sum(ll.raw[ABC.idx,2]) +
            sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,3])
        ll.CAB <- sum(ll.raw[ABC.idx,2]) +
            sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,1])

        ll <- c(ll.ABC, ll.BCA, ll.ACB, ll.BAC, ll.CBA, ll.CAB)

        return(-ll[l])
}

blOptimFuncAdj <- function(par, data, method, cABC, tau2bound, N, theta, l){
    tau1 <- 0
    par <- exp(par)

    if(!(par[1] >= 0 & par[1] <= tau2bound)){
        return(Inf)
    }
    if(!(par[2] >= 0 & par[2] <= 4*N/theta)){
        return(Inf)
    }

    tau2 <- par[1]
    m <- par[2]
    t1.vec <- data[,4]
    t2.vec <- data[,5]

    ABC.idx <- data[,1] == 1
    BCA.idx <- data[,2] == 1
    ACB.idx <- data[,3] == 1

    QAB <- matrix(data = c(-2*m, 0, m, m, 0, 0,
                           0, -2*m, m, m, 0, 0,
                           m, m, -2*m-cABC, 0, cABC, 0,
                           m, m, 0, -2*m-cABC, 0, cABC,
                           0, 0, 0, 0, -m, m,
                           0, 0, 0, 0, m, -m), nrow = 6, byrow = TRUE)

    QABC <- matrix(data = c(-3*m-3*(cABC + cABC + cABC), m, m, m, 0, 0, 0,
                            0, 3*cABC, 0, 0, 0, 3*cABC, 0, 0, 0, 3*cABC, 0, 0,
                            0, 0, 0, m,-3*m-cABC,0,0,m,m,0,0,0,cABC,0,0,0,0,0,
                            0,0,0,0,0,0,0,m,0,-3*m-cABC,0,m,0,m,0,0,0,0,0,0,
                            cABC,0,0,0,0,0,0,0,0,m,0,0,-3*m-cABC,0,m,m,0,0,0,0,
                            0,0,0,0,0,0,cABC,0,0,0,0,0,m,m,0,-3*m-cABC,0,0,m,0,
                            0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,m,0,m,0,-3*m-cABC,0,
                            m,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,m,m,0,0,
                            -3*m-cABC,m,0,0,cABC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            m,m,m,-3*m-3*(cABC + cABC + cABC),0,0,0,3*cABC,0,0,
                            0,3*cABC,0,0,0,3*cABC,0,0,0,0,0,0,0,0,0,0,-2*m-cABC,
                            m,m,0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,m,
                            -2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m,0,
                            -2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m,m,
                            -2*m-cABC,0,0,0,0,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,
                            0,0,0,-2*m-cABC,m,m,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,m,-2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,m,0,-2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,m,m,-2*m-cABC,0,0,0,0,0,cABC,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,-2*m-cABC,m,m,0,cABC,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,m,-2*m,0,m,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,m,0,-2*m,m,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,m,m,-2*m-cABC,0,cABC,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0,-m,m,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,m,-m), nrow = 22, byrow = TRUE)

    QAB <- matrix(unlist(QAB), nrow = 6)
    QABC <- matrix(unlist(QABC), nrow = 22)
    epsilon <- 10^(-15)
    reps <- nrow(data)

    ll.raw <- matrix(rep(NA,3*reps), ncol = 3)
    for(i in 1:reps){
        ll.raw[i,1] <- log(f1Grid(t1.vec[i], t2.vec[i],
                                  tau1, tau2, QAB, QABC, cABC) +
                               f2Grid(t1.vec[i], t2.vec[i],
                                      tau1, tau2, QAB, QABC, cABC) + epsilon)
        ll.raw[i,2] <- log(f3Grid(t1.vec[i], t2.vec[i],
                                  tau1, tau2, QAB, QABC, cABC) + epsilon)
        ll.raw[i,3] <- log(f4Grid(t1.vec[i], t2.vec[i],
                                  tau1, tau2, QAB, QABC, cABC) + epsilon)
    }

    ll.ABC <- sum(ll.raw[ABC.idx,1]) +
        sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,3])
    ll.BCA <- sum(ll.raw[ABC.idx,3]) +
        sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,2])
    ll.ACB <- sum(ll.raw[ABC.idx,3]) +
        sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,1])
    ll.BAC <- sum(ll.raw[ABC.idx,1]) +
        sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,2])
    ll.CBA <- sum(ll.raw[ABC.idx,2]) +
        sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,3])
    ll.CAB <- sum(ll.raw[ABC.idx,2]) +
        sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,1])

    ll <- c(ll.ABC, ll.BCA, ll.ACB, ll.BAC, ll.CBA, ll.CAB)

    return(-ll[l])
}


