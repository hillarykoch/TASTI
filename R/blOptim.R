# #############################################
# ##
# ## Distribution of gene tree histories
# ##
# #############################################
# library(expm)
#
# t.pos <- gregexpr(pattern = 't', basename(getwd()))[[1]][1]
# reps <- substr(basename(getwd()),1,t.pos-1)
# iter <- 100
#
# ## Count the number of trees simulated, not from our super naive sampling method
# trees <- read.csv(paste0("counts_branchlens",reps,".csv"), header = F)
#
# ## Break the data into matrices corresponding to each of the 100 iterations
# reps <- as.integer(reps)
# treelist <- list()
# tau1bound <- tau2bound <- rep(NA,iter)
# for(i in 1:iter){
#     treelist[[i]] <- trees[(reps*(i-1)+1):((reps*(i-1)+1)+reps-1),]
#     tau1bound[i] <- min(treelist[[i]][,4]) # tau1 must be less than this (tau1 comes before t1 always)
#     tau2bound[i] <- min(treelist[[i]][,5]) # tau2 must be less than this (cant have 2 coalescence on internal branch)
# }
#
#
# theta <- .005
# cAB1 <- 2/theta
# cAB2 <- cAB1
# cAC1 <- cAB1
# cAC2 <- cAB1
# cBC1 <- cAB1
# cBC2 <- cAB1
# cABC1 <- cAB1
# cABC2 <- cAB1
# m1 <- NA
# m2 <- m1
#
# # Use this function when A enters in subpop1, B enters in subpop2, C enters in subpop2
# func <- function(par,k,m){
#     par <- exp(par)
#     # First, make sure everything is in bounds
#     if(!(par[1] >= 0 & par[1] <= tau1bound[k])){
#         return(Inf)
#     }
#
#     if(!(par[2] >=0 & par[2] <= tau2bound[k])){
#         return(Inf)
#     }
#
#     if(!(par[3] >= 0 & par[3] <= 10^7)){
#         return(Inf)
#     }
#
#     # Now make sure everything is correct just in case
#     tau1 <- min(par[1:2])
#     tau2 <- max(par[1:2])
#     m1 <- m2 <- par[3]
#     t1.vec <- treelist[[k]][,4]
#     t2.vec <- treelist[[k]][,5]
#
#     ## Get the indices at which each topology was observed
#     ## for this particular iteration
#     ABC.idx <- treelist[[k]][,1] == 1
#     BCA.idx <- treelist[[k]][,2] == 1
#     ACB.idx <- treelist[[k]][,3] == 1
#
#
#     QAB <- matrix(data = c(-2*m1,  0,        m1,        m1,    0,     0,
#                            0,  -2*m1,        m1,        m1,    0,     0,
#                            m1,     m1,-2*m1-cAB1,         0, cAB1,    0,
#                            m1,     m1,         0,-2*m1-cAB2,    0, cAB2,
#                            0,      0,          0,         0,  -m2,   m2,
#                            0,      0,          0,         0,   m1,  -m1), nrow = 6, byrow = TRUE)
#
#
#     QABC <- matrix(data = c(-3*m2-3*(cAB1 + cAC1 + cBC1),m2,m2,m2,0,0,0,0,3*cAB1,0,0,0,3*cAC1,0,0,0,3*cBC1,0,0,0,0,0,
#                             m2,-3*m2-cAB1,0,0,m2,m2,0,0,0,cAB1,0,0,0,0,0,0,0,0,0,0,0,0,
#                             m2,0,-3*m2-cAC1,0,m2,0,m2,0,0,0,0,0,0,cAC1,0,0,0,0,0,0,0,0,
#                             m2,0,0,-3*m2-cBC1,0,m2,m2,0,0,0,0,0,0,0,0,0,0,cBC1,0,0,0,0,
#                             0,m2,m2,0,-3*m2-cBC2,0,0,m2,0,0,0,0,0,0,0,0,0,0,cBC2,0,0,0,
#                             0,m2,0,m2,0,-3*m2-cAC2,0,m2,0,0,0,0,0,0,cAC2,0,0,0,0,0,0,0,
#                             0,0,m2,m2,0,0,-3*m2-cAB2,m2,0,0,cAB2,0,0,0,0,0,0,0,0,0,0,0,
#                             0,0,0,0,m2,m2,m2,-3*m2-3*(cAB2 + cAC2 + cBC2),0,0,0,3*cAB2,0,0,0,3*cAC2,0,0,0,3*cBC2,0,0,
#                             0,0,0,0,0,0,0,0,-2*m2-cABC1,m2,m2,0,0,0,0,0,0,0,0,0,cABC1,0,
#                             0,0,0,0,0,0,0,0,m2,-2*m2,0,m2,0,0,0,0,0,0,0,0,0,0,
#                             0,0,0,0,0,0,0,0,m2,0,-2*m2,m2,0,0,0,0,0,0,0,0,0,0,
#                             0,0,0,0,0,0,0,0,0,m2,m2,-2*m2-cABC2,0,0,0,0,0,0,0,0,0,cABC2,
#                             0,0,0,0,0,0,0,0,0,0,0,0,-2*m2-cABC1,m2,m2,0,0,0,0,0,cABC1,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,m2,-2*m2,0,m2,0,0,0,0,0,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,m2,0,-2*m2,m2,0,0,0,0,0,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,m2,m2,-2*m2-cABC2,0,0,0,0,0,cABC2,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2*m2-cABC1,m2,m2,0,cABC1,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m2,-2*m2,0,m2,0,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m2,0,-2*m2,m2,0,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m2,m2,-2*m2-cABC2,0,cABC2,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-m2,m2,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m2,-m2), nrow = 22, byrow = TRUE)
#
#     ## These two general functions should be fine for getting all the necessary integrals
#     ## coal1 is for the first coalescence above the root
#     coal1 <- function(T1, state1, state2, rate){
#         #rate*expm(QABC*(T1-lower))[state1,state2]
#         rate*T1[state1,state2]
#     }
#
#     ## Function outputs probability of coalescence along internal branch
#     f1 <- function(t1,t2) {
#         if(tau2 < t1){
#             return(0)
#         } else{
#             t1t1 <- expm(QAB*(t1-tau1))
#             t2t1 <- expm(QAB*(tau2-t1))
#             t1t2 <- expm(QABC*(t2-tau2))
#
#             c.ab <-
#                 cAB1*t1t1[1,3]*(t2t1[5,5]*(cABC1*t1t2[10,9] + cABC2*t1t2[10,12]) +
#                                     t2t1[5,6]*(cABC1*t1t2[12,9] + cABC2*t1t2[12,12])) +
#                 cAB2*t1t1[1,4]*(t2t1[6,5]*(cABC1*t1t2[10,9] + cABC2*t1t2[10,12]) +
#                                     t2t1[6,6]*(cABC1*t1t2[12,9] + cABC2*t1t2[12,12]))
#             return(c.ab)
#         }
#     }
#
#     ## Function outputs probability of AB coalescing first, above the root
#     f2 <- function(t1,t2) {
#         if(tau2 > t1){
#             return(0)
#         } else {
#             val <- expm(QAB*(tau2-tau1))[1,1:4]
#             norm.c <- sum(val)
#
#             alpha1 <- val[1]/norm.c
#             alpha2 <- val[2]/norm.c
#             alpha3 <- val[3]/norm.c
#             alpha4 <- val[4]/norm.c
#
#             t1t2 <- expm(QABC*(t1-tau2))
#             t2t1 <- expm(QABC*(t2-t1))
#
#             c.ab <-
#                 alpha1*(coal1(t1t2,5,1,3*cABC1)*(coal1(t2t1,9,9,cAB1) + coal1(t2t1,9,12,cAB2)) +
#                             coal1(t1t2,5,2,cABC1)*(coal1(t2t1,10,9,cAB1) + coal1(t2t1,10,12,cAB2)) +
#                             coal1(t1t2,5,7,cABC2)*(coal1(t2t1,11,9,cAB1) + coal1(t2t1,11,12,cAB2)) +
#                             coal1(t1t2,5,8,3*cABC2)*(coal1(t2t1,12,9,cAB1) + coal1(t2t1,12,12,cAB2))) +
#                 alpha2*(coal1(t1t2,6,1,3*cABC1)*(coal1(t2t1,9,9,cAB1) + coal1(t2t1,9,12,cAB2)) +
#                             coal1(t1t2,6,2,cABC1)*(coal1(t2t1,10,9,cAB1) + coal1(t2t1,10,12,cAB2)) +
#                             coal1(t1t2,6,7,cABC2)*(coal1(t2t1,11,9,cAB1) + coal1(t2t1,11,12,cAB2)) +
#                             coal1(t1t2,6,8,3*cABC2)*(coal1(t2t1,12,9,cAB1) + coal1(t2t1,12,12,cAB2))) +
#                 alpha3*(coal1(t1t2,2,1,3*cABC1)*(coal1(t2t1,9,9,cAB1) + coal1(t2t1,9,12,cAB2)) +
#                             coal1(t1t2,2,2,cABC1)*(coal1(t2t1,10,9,cAB1) + coal1(t2t1,10,12,cAB2)) +
#                             coal1(t1t2,2,7,cABC2)*(coal1(t2t1,11,9,cAB1) + coal1(t2t1,11,12,cAB2)) +
#                             coal1(t1t2,2,8,3*cABC2)*(coal1(t2t1,12,9,cAB1) + coal1(t2t1,12,12,cAB2))) +
#                 alpha4*(coal1(t1t2,8,1,3*cABC1)*(coal1(t2t1,9,9,cAB1) + coal1(t2t1,9,12,cAB2)) +
#                             coal1(t1t2,8,2,cABC1)*(coal1(t2t1,10,9,cAB1) + coal1(t2t1,10,12,cAB2)) +
#                             coal1(t1t2,8,7,cABC2)*(coal1(t2t1,11,9,cAB1) + coal1(t2t1,11,12,cAB2)) +
#                             coal1(t1t2,8,8,3*cABC2)*(coal1(t2t1,12,9,cAB1) + coal1(t2t1,12,12,cAB2)))
#
#             return(c.ab)
#         }
#     }
#
#     ## Function outputs probability of BC coalescing first, above the root
#     f3 <- function(t1,t2) {
#         if(tau2 > t1){
#             return(0)
#         } else {
#             val <- expm(QAB*(tau2-tau1))[1,1:4]
#             norm.c = sum(val)
#
#             alpha1 <- val[1]/norm.c
#             alpha2 <- val[2]/norm.c
#             alpha3 <- val[3]/norm.c
#             alpha4 <- val[4]/norm.c
#
#             t1t2 <- expm(QABC*(t1-tau2))
#             t2t1 <- expm(QABC*(t2-t1))
#
#             c.bc <-
#                 alpha1*(coal1(t1t2,5,1,3*cABC1)*(coal1(t2t1,17,17,cAB1) + coal1(t2t1,17,20,cAB2)) +
#                             coal1(t1t2,5,4,cABC1)*(coal1(t2t1,18,17,cAB1) + coal1(t2t1,18,20,cAB2)) +
#                             coal1(t1t2,5,5,cABC2)*(coal1(t2t1,19,17,cAB1) + coal1(t2t1,19,20,cAB2)) +
#                             coal1(t1t2,5,8,3*cABC2)*(coal1(t2t1,20,17,cAB1) + coal1(t2t1,20,20,cAB2))) +
#                 alpha2*(coal1(t1t2,6,1,3*cABC1)*(coal1(t2t1,17,17,cAB1) + coal1(t2t1,17,20,cAB2)) +
#                             coal1(t1t2,6,4,cABC1)*(coal1(t2t1,18,17,cAB1) + coal1(t2t1,18,20,cAB2)) +
#                             coal1(t1t2,6,5,cABC2)*(coal1(t2t1,19,17,cAB1) + coal1(t2t1,19,20,cAB2)) +
#                             coal1(t1t2,6,8,3*cABC2)*(coal1(t2t1,20,17,cAB1) + coal1(t2t1,20,20,cAB2))) +
#                 alpha3*(coal1(t1t2,2,1,3*cABC1)*(coal1(t2t1,17,17,cAB1) + coal1(t2t1,17,20,cAB2)) +
#                             coal1(t1t2,2,4,cABC1)*(coal1(t2t1,18,17,cAB1) + coal1(t2t1,18,20,cAB2)) +
#                             coal1(t1t2,2,5,cABC2)*(coal1(t2t1,19,17,cAB1) + coal1(t2t1,19,20,cAB2)) +
#                             coal1(t1t2,2,8,3*cABC2)*(coal1(t2t1,20,17,cAB1) + coal1(t2t1,20,20,cAB2))) +
#                 alpha4*(coal1(t1t2,8,1,3*cABC1)*(coal1(t2t1,17,17,cAB1) + coal1(t2t1,17,20,cAB2)) +
#                             coal1(t1t2,8,4,cABC1)*(coal1(t2t1,18,17,cAB1) + coal1(t2t1,18,20,cAB2)) +
#                             coal1(t1t2,8,5,cABC2)*(coal1(t2t1,19,17,cAB1) + coal1(t2t1,19,20,cAB2)) +
#                             coal1(t1t2,8,8,3*cABC2)*(coal1(t2t1,20,17,cAB1) + coal1(t2t1,20,20,cAB2)))
#
#             return(c.bc)
#         }
#     }
#
#     ## Function outputs probability of AC coalescing first, above the root
#     f4 <- function(t1,t2) {
#         if(tau2 > t1){
#             return(0)
#         } else{
#             val <- expm(QAB*(tau2-tau1))[1,1:4]
#             norm.c <- sum(val)
#
#             alpha1 <- val[1]/norm.c
#             alpha2 <- val[2]/norm.c
#             alpha3 <- val[3]/norm.c
#             alpha4 <- val[4]/norm.c
#
#             t1t2 <- expm(QABC*(t1-tau2))
#             t2t1 <- expm(QABC*(t2-t1))
#
#             c.ac <-
#                 alpha1*(coal1(t1t2,5,1,3*cABC1)*(coal1(t2t1,13,13,cAB1) + coal1(t2t1,13,16,cAB2)) +
#                             coal1(t1t2,5,3,cABC1)*(coal1(t2t1,15,13,cAB1) + coal1(t2t1,15,16,cAB2)) +
#                             coal1(t1t2,5,6,cABC2)*(coal1(t2t1,14,13,cAB1) + coal1(t2t1,14,16,cAB2)) +
#                             coal1(t1t2,5,8,3*cABC2)*(coal1(t2t1,16,13,cAB1) + coal1(t2t1,16,16,cAB2))) +
#                 alpha2*(coal1(t1t2,6,1,3*cABC1)*(coal1(t2t1,13,13,cAB1) + coal1(t2t1,13,16,cAB2)) +
#                             coal1(t1t2,6,3,cABC1)*(coal1(t2t1,15,13,cAB1) + coal1(t2t1,15,16,cAB2)) +
#                             coal1(t1t2,6,6,cABC2)*(coal1(t2t1,14,13,cAB1) + coal1(t2t1,14,16,cAB2)) +
#                             coal1(t1t2,6,8,3*cABC2)*(coal1(t2t1,16,13,cAB1) + coal1(t2t1,16,16,cAB2))) +
#                 alpha3*(coal1(t1t2,2,1,3*cABC1)*(coal1(t2t1,13,13,cAB1) + coal1(t2t1,13,16,cAB2)) +
#                             coal1(t1t2,2,3,cABC1)*(coal1(t2t1,15,13,cAB1) + coal1(t2t1,15,16,cAB2)) +
#                             coal1(t1t2,2,6,cABC2)*(coal1(t2t1,14,13,cAB1) + coal1(t2t1,14,16,cAB2)) +
#                             coal1(t1t2,2,8,3*cABC2)*(coal1(t2t1,16,13,cAB1) + coal1(t2t1,16,16,cAB2))) +
#                 alpha4*(coal1(t1t2,8,1,3*cABC1)*(coal1(t2t1,13,13,cAB1) + coal1(t2t1,13,16,cAB2)) +
#                             coal1(t1t2,8,3,cABC1)*(coal1(t2t1,15,13,cAB1) + coal1(t2t1,15,16,cAB2)) +
#                             coal1(t1t2,8,6,cABC2)*(coal1(t2t1,14,13,cAB1) + coal1(t2t1,14,16,cAB2)) +
#                             coal1(t1t2,8,8,3*cABC2)*(coal1(t2t1,16,13,cAB1) + coal1(t2t1,16,16,cAB2)))
#
#             return(c.ac)
#         }
#     }
#
#     epsilon <- 10^(-15)
#     ## Put the probabilities of each tree into one nice vector
#     ll.raw <- matrix(rep(NA,3*reps), ncol = 3)
#     for(i in 1:reps){
#         ll.raw[i,1] <- log(f1(t1.vec[i], t2.vec[i]) + f2(t1.vec[i], t2.vec[i]) + epsilon)
#         ll.raw[i,2] <- log(f3(t1.vec[i], t2.vec[i]) + epsilon)
#         ll.raw[i,3] <- log(f4(t1.vec[i], t2.vec[i]) + epsilon)
#     }
#
#     ll.ABC <- sum(ll.raw[ABC.idx,1]) + sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,3])
#     ll.BCA <- sum(ll.raw[ABC.idx,3]) + sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,2])
#     ll.ACB <- sum(ll.raw[ABC.idx,3]) + sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,1])
#     ll.BAC <- sum(ll.raw[ABC.idx,1]) + sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,2])
#     ll.CBA <- sum(ll.raw[ABC.idx,2]) + sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,3])
#     ll.CAB <- sum(ll.raw[ABC.idx,2]) + sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,1])
#
#     ll <- c(ll.ABC, ll.BCA, ll.ACB, ll.BAC, ll.CBA, ll.CAB)
#
#     return(-ll[m]) ## use negative because optim minimizes
# }
#
# func.adj <- function(par,k,m){
#     tau1 = 0
#     par <- exp(par)
#     # First, make sure everything is in bounds
#     if(!(par[1] >= 0 & par[1] <= tau2bound[k])){
#         return(Inf)
#     }
#     if(!(par[2] >= 0 & par[2] <= 10^7)){
#         return(Inf)
#     }
#
#     # Now make sure everything is correct just in case
#     tau2 <- par[1]
#     m1 <- m2 <- par[2]
#
#     t1.vec <- treelist[[k]][,4]
#     t2.vec <- treelist[[k]][,5]
#
#     ## Get the indices at which each topology was observed
#     ## for this particular iteration
#     ABC.idx <- treelist[[k]][,1] == 1
#     BCA.idx <- treelist[[k]][,2] == 1
#     ACB.idx <- treelist[[k]][,3] == 1
#
#     QAB <- matrix(data = c(-2*m1,      0,        m1,        m1,    0,    0,
#                            0,  -2*m1,        m1,        m1,    0,    0,
#                            m1,     m1,-2*m1-cAB1,         0, cAB1,    0,
#                            m1,     m1,         0,-2*m1-cAB2,    0, cAB2,
#                            0,      0,         0,         0,  -m2,   m2,
#                            0,      0,         0,         0,   m1,  -m1), nrow = 6, byrow = TRUE)
#
#
#     QABC <- matrix(data = c(-3*m2-3*(cAB1 + cAC1 + cBC1),m2,m2,m2,0,0,0,0,3*cAB1,0,0,0,3*cAC1,0,0,0,3*cBC1,0,0,0,0,0,
#                             m2,-3*m2-cAB1,0,0,m2,m2,0,0,0,cAB1,0,0,0,0,0,0,0,0,0,0,0,0,
#                             m2,0,-3*m2-cAC1,0,m2,0,m2,0,0,0,0,0,0,cAC1,0,0,0,0,0,0,0,0,
#                             m2,0,0,-3*m2-cBC1,0,m2,m2,0,0,0,0,0,0,0,0,0,0,cBC1,0,0,0,0,
#                             0,m2,m2,0,-3*m2-cBC2,0,0,m2,0,0,0,0,0,0,0,0,0,0,cBC2,0,0,0,
#                             0,m2,0,m2,0,-3*m2-cAC2,0,m2,0,0,0,0,0,0,cAC2,0,0,0,0,0,0,0,
#                             0,0,m2,m2,0,0,-3*m2-cAB2,m2,0,0,cAB2,0,0,0,0,0,0,0,0,0,0,0,
#                             0,0,0,0,m2,m2,m2,-3*m2-3*(cAB2 + cAC2 + cBC2),0,0,0,3*cAB2,0,0,0,3*cAC2,0,0,0,3*cBC2,0,0,
#                             0,0,0,0,0,0,0,0,-2*m2-cABC1,m2,m2,0,0,0,0,0,0,0,0,0,cABC1,0,
#                             0,0,0,0,0,0,0,0,m2,-2*m2,0,m2,0,0,0,0,0,0,0,0,0,0,
#                             0,0,0,0,0,0,0,0,m2,0,-2*m2,m2,0,0,0,0,0,0,0,0,0,0,
#                             0,0,0,0,0,0,0,0,0,m2,m2,-2*m2-cABC2,0,0,0,0,0,0,0,0,0,cABC2,
#                             0,0,0,0,0,0,0,0,0,0,0,0,-2*m2-cABC1,m2,m2,0,0,0,0,0,cABC1,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,m2,-2*m2,0,m2,0,0,0,0,0,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,m2,0,-2*m2,m2,0,0,0,0,0,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,m2,m2,-2*m2-cABC2,0,0,0,0,0,cABC2,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2*m2-cABC1,m2,m2,0,cABC1,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m2,-2*m2,0,m2,0,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m2,0,-2*m2,m2,0,0,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m2,m2,-2*m2-cABC2,0,cABC2,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-m2,m2,
#                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,m2,-m2), nrow = 22, byrow = TRUE)
#
#     ## These two general functions should be fine for getting all the necessary integrals
#     ## coal1 is for the first coalescence above the root
#     coal1 <- function(T1, state1, state2, rate){
#         #rate*expm(QABC*(T1-lower))[state1,state2]
#         rate*T1[state1,state2]
#     }
#
#
#     ## Function outputs probability of coalescence along internal branch
#     f1 <- function(t1,t2) {
#         if(tau2 < t1){
#             return(0)
#         } else{
#             t1t1 <- expm(QAB*(t1-tau1))
#             t2t1 <- expm(QAB*(tau2-t1))
#             t1t2 <- expm(QABC*(t2-tau2))
#
#             c.ab <-
#                 cAB1*t1t1[1,3]*(t2t1[5,5]*(cABC1*t1t2[10,9] + cABC2*t1t2[10,12]) +
#                                     t2t1[5,6]*(cABC1*t1t2[12,9] + cABC2*t1t2[12,12])) +
#                 cAB2*t1t1[1,4]*(t2t1[6,5]*(cABC1*t1t2[10,9] + cABC2*t1t2[10,12]) +
#                                     t2t1[6,6]*(cABC1*t1t2[12,9] + cABC2*t1t2[12,12]))
#             return(c.ab)
#         }
#     }
#
#     ## Function outputs probability of AB coalescing first, above the root
#     f2 <- function(t1,t2) {
#         if(tau2 > t1){
#             return(0)
#         } else {
#             val <- expm(QAB*(tau2-tau1))[1,1:4]
#             norm.c <- sum(val)
#
#             alpha1 <- val[1]/norm.c
#             alpha2 <- val[2]/norm.c
#             alpha3 <- val[3]/norm.c
#             alpha4 <- val[4]/norm.c
#
#             t1t2 <- expm(QABC*(t1-tau2))
#             t2t1 <- expm(QABC*(t2-t1))
#
#             c.ab <-
#                 alpha1*(coal1(t1t2,5,1,3*cABC1)*(coal1(t2t1,9,9,cAB1) + coal1(t2t1,9,12,cAB2)) +
#                             coal1(t1t2,5,2,cABC1)*(coal1(t2t1,10,9,cAB1) + coal1(t2t1,10,12,cAB2)) +
#                             coal1(t1t2,5,7,cABC2)*(coal1(t2t1,11,9,cAB1) + coal1(t2t1,11,12,cAB2)) +
#                             coal1(t1t2,5,8,3*cABC2)*(coal1(t2t1,12,9,cAB1) + coal1(t2t1,12,12,cAB2))) +
#                 alpha2*(coal1(t1t2,6,1,3*cABC1)*(coal1(t2t1,9,9,cAB1) + coal1(t2t1,9,12,cAB2)) +
#                             coal1(t1t2,6,2,cABC1)*(coal1(t2t1,10,9,cAB1) + coal1(t2t1,10,12,cAB2)) +
#                             coal1(t1t2,6,7,cABC2)*(coal1(t2t1,11,9,cAB1) + coal1(t2t1,11,12,cAB2)) +
#                             coal1(t1t2,6,8,3*cABC2)*(coal1(t2t1,12,9,cAB1) + coal1(t2t1,12,12,cAB2))) +
#                 alpha3*(coal1(t1t2,2,1,3*cABC1)*(coal1(t2t1,9,9,cAB1) + coal1(t2t1,9,12,cAB2)) +
#                             coal1(t1t2,2,2,cABC1)*(coal1(t2t1,10,9,cAB1) + coal1(t2t1,10,12,cAB2)) +
#                             coal1(t1t2,2,7,cABC2)*(coal1(t2t1,11,9,cAB1) + coal1(t2t1,11,12,cAB2)) +
#                             coal1(t1t2,2,8,3*cABC2)*(coal1(t2t1,12,9,cAB1) + coal1(t2t1,12,12,cAB2))) +
#                 alpha4*(coal1(t1t2,8,1,3*cABC1)*(coal1(t2t1,9,9,cAB1) + coal1(t2t1,9,12,cAB2)) +
#                             coal1(t1t2,8,2,cABC1)*(coal1(t2t1,10,9,cAB1) + coal1(t2t1,10,12,cAB2)) +
#                             coal1(t1t2,8,7,cABC2)*(coal1(t2t1,11,9,cAB1) + coal1(t2t1,11,12,cAB2)) +
#                             coal1(t1t2,8,8,3*cABC2)*(coal1(t2t1,12,9,cAB1) + coal1(t2t1,12,12,cAB2)))
#
#             return(c.ab)
#         }
#     }
#
#     ## Function outputs probability of BC coalescing first, above the root
#     f3 <- function(t1,t2) {
#         if(tau2 > t1){
#             return(0)
#         } else {
#             val <- expm(QAB*(tau2-tau1))[1,1:4]
#             norm.c = sum(val)
#
#             alpha1 <- val[1]/norm.c
#             alpha2 <- val[2]/norm.c
#             alpha3 <- val[3]/norm.c
#             alpha4 <- val[4]/norm.c
#
#             t1t2 <- expm(QABC*(t1-tau2))
#             t2t1 <- expm(QABC*(t2-t1))
#
#             c.bc <-
#                 alpha1*(coal1(t1t2,5,1,3*cABC1)*(coal1(t2t1,17,17,cAB1) + coal1(t2t1,17,20,cAB2)) +
#                             coal1(t1t2,5,4,cABC1)*(coal1(t2t1,18,17,cAB1) + coal1(t2t1,18,20,cAB2)) +
#                             coal1(t1t2,5,5,cABC2)*(coal1(t2t1,19,17,cAB1) + coal1(t2t1,19,20,cAB2)) +
#                             coal1(t1t2,5,8,3*cABC2)*(coal1(t2t1,20,17,cAB1) + coal1(t2t1,20,20,cAB2))) +
#                 alpha2*(coal1(t1t2,6,1,3*cABC1)*(coal1(t2t1,17,17,cAB1) + coal1(t2t1,17,20,cAB2)) +
#                             coal1(t1t2,6,4,cABC1)*(coal1(t2t1,18,17,cAB1) + coal1(t2t1,18,20,cAB2)) +
#                             coal1(t1t2,6,5,cABC2)*(coal1(t2t1,19,17,cAB1) + coal1(t2t1,19,20,cAB2)) +
#                             coal1(t1t2,6,8,3*cABC2)*(coal1(t2t1,20,17,cAB1) + coal1(t2t1,20,20,cAB2))) +
#                 alpha3*(coal1(t1t2,2,1,3*cABC1)*(coal1(t2t1,17,17,cAB1) + coal1(t2t1,17,20,cAB2)) +
#                             coal1(t1t2,2,4,cABC1)*(coal1(t2t1,18,17,cAB1) + coal1(t2t1,18,20,cAB2)) +
#                             coal1(t1t2,2,5,cABC2)*(coal1(t2t1,19,17,cAB1) + coal1(t2t1,19,20,cAB2)) +
#                             coal1(t1t2,2,8,3*cABC2)*(coal1(t2t1,20,17,cAB1) + coal1(t2t1,20,20,cAB2))) +
#                 alpha4*(coal1(t1t2,8,1,3*cABC1)*(coal1(t2t1,17,17,cAB1) + coal1(t2t1,17,20,cAB2)) +
#                             coal1(t1t2,8,4,cABC1)*(coal1(t2t1,18,17,cAB1) + coal1(t2t1,18,20,cAB2)) +
#                             coal1(t1t2,8,5,cABC2)*(coal1(t2t1,19,17,cAB1) + coal1(t2t1,19,20,cAB2)) +
#                             coal1(t1t2,8,8,3*cABC2)*(coal1(t2t1,20,17,cAB1) + coal1(t2t1,20,20,cAB2)))
#
#             return(c.bc)
#         }
#     }
#
#     ## Function outputs probability of AC coalescing first, above the root
#     f4 <- function(t1,t2) {
#         if(tau2 > t1){
#             return(0)
#         } else{
#             val <- expm(QAB*(tau2-tau1))[1,1:4]
#             norm.c <- sum(val)
#
#             alpha1 <- val[1]/norm.c
#             alpha2 <- val[2]/norm.c
#             alpha3 <- val[3]/norm.c
#             alpha4 <- val[4]/norm.c
#
#             t1t2 <- expm(QABC*(t1-tau2))
#             t2t1 <- expm(QABC*(t2-t1))
#
#             c.ac <-
#                 alpha1*(coal1(t1t2,5,1,3*cABC1)*(coal1(t2t1,13,13,cAB1) + coal1(t2t1,13,16,cAB2)) +
#                             coal1(t1t2,5,3,cABC1)*(coal1(t2t1,15,13,cAB1) + coal1(t2t1,15,16,cAB2)) +
#                             coal1(t1t2,5,6,cABC2)*(coal1(t2t1,14,13,cAB1) + coal1(t2t1,14,16,cAB2)) +
#                             coal1(t1t2,5,8,3*cABC2)*(coal1(t2t1,16,13,cAB1) + coal1(t2t1,16,16,cAB2))) +
#                 alpha2*(coal1(t1t2,6,1,3*cABC1)*(coal1(t2t1,13,13,cAB1) + coal1(t2t1,13,16,cAB2)) +
#                             coal1(t1t2,6,3,cABC1)*(coal1(t2t1,15,13,cAB1) + coal1(t2t1,15,16,cAB2)) +
#                             coal1(t1t2,6,6,cABC2)*(coal1(t2t1,14,13,cAB1) + coal1(t2t1,14,16,cAB2)) +
#                             coal1(t1t2,6,8,3*cABC2)*(coal1(t2t1,16,13,cAB1) + coal1(t2t1,16,16,cAB2))) +
#                 alpha3*(coal1(t1t2,2,1,3*cABC1)*(coal1(t2t1,13,13,cAB1) + coal1(t2t1,13,16,cAB2)) +
#                             coal1(t1t2,2,3,cABC1)*(coal1(t2t1,15,13,cAB1) + coal1(t2t1,15,16,cAB2)) +
#                             coal1(t1t2,2,6,cABC2)*(coal1(t2t1,14,13,cAB1) + coal1(t2t1,14,16,cAB2)) +
#                             coal1(t1t2,2,8,3*cABC2)*(coal1(t2t1,16,13,cAB1) + coal1(t2t1,16,16,cAB2))) +
#                 alpha4*(coal1(t1t2,8,1,3*cABC1)*(coal1(t2t1,13,13,cAB1) + coal1(t2t1,13,16,cAB2)) +
#                             coal1(t1t2,8,3,cABC1)*(coal1(t2t1,15,13,cAB1) + coal1(t2t1,15,16,cAB2)) +
#                             coal1(t1t2,8,6,cABC2)*(coal1(t2t1,14,13,cAB1) + coal1(t2t1,14,16,cAB2)) +
#                             coal1(t1t2,8,8,3*cABC2)*(coal1(t2t1,16,13,cAB1) + coal1(t2t1,16,16,cAB2)))
#
#             return(c.ac)
#         }
#     }
#
#     epsilon <- 10^(-15)
#     ## Put the probabilities of each tree into one nice vector
#     ll.raw <- matrix(rep(NA,3*reps), ncol = 3)
#     for(i in 1:reps){
#         ll.raw[i,1] <- log(f1(t1.vec[i], t2.vec[i]) + f2(t1.vec[i], t2.vec[i]) + epsilon)
#         ll.raw[i,2] <- log(f3(t1.vec[i], t2.vec[i]) + epsilon)
#         ll.raw[i,3] <- log(f4(t1.vec[i], t2.vec[i]) + epsilon)
#     }
#
#     ll.ABC <- sum(ll.raw[ABC.idx,1]) + sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,3])
#     ll.BCA <- sum(ll.raw[ABC.idx,3]) + sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,2])
#     ll.ACB <- sum(ll.raw[ABC.idx,3]) + sum(ll.raw[BCA.idx,2]) + sum(ll.raw[ACB.idx,1])
#     ll.BAC <- sum(ll.raw[ABC.idx,1]) + sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,2])
#     ll.CBA <- sum(ll.raw[ABC.idx,2]) + sum(ll.raw[BCA.idx,1]) + sum(ll.raw[ACB.idx,3])
#     ll.CAB <- sum(ll.raw[ABC.idx,2]) + sum(ll.raw[BCA.idx,3]) + sum(ll.raw[ACB.idx,1])
#
#     ll <- c(ll.ABC, ll.BCA, ll.ACB, ll.BAC, ll.CBA, ll.CAB)
#
#     return(-ll[m]) ## use negative because optim minimizes
# }
#
#
# library(foreach)
# library(doParallel)
# library(parallel)
#
# # setup parallel backend to use given number of cores
# cores <- 20
# cluster <- makeCluster(cores, outfile = "log.txt")
# registerDoParallel(cluster)
# clusterEvalQ(cluster, library(expm))
# clusterEvalQ(cluster, library(cubature))
#
# like.ests <- foreach(j = 1:iter, .combine = 'cbind', .packages = "foreach") %dopar% {
#     foreach(l = 1:6, .combine = 'rbind') %do% {
#         if(tau1bound[j] == tau2bound[j] & tau1bound[j] == 0){
#             # If there is a "star" tree
#             c(0,0,NA,-Inf)
#         } else if(tau1bound[j] == 0){
#             # Adjust for tau1bound=0 by fixing tau1 and optimizing over tau2 and m
#             init <- c(runif(1, min = log(10^(-10)), max = log(tau2bound[j])), runif(1, min=-7, max=7))
#             opt <- optim(par = init, fn = func.adj, method = "Nelder-Mead", control = list(trace = F), k=j, m=l)
#
#             # Get back optimal values
#             par <- exp(opt$par)
#             est1 <- 0
#             est2 <- par[1]
#             est3 <- par[2]
#             obj <- opt$value
#
#             c(est1,est2,est3,obj)
#         } else{
#             init <- c(runif(1, min = log(10^(-10)), max = log(tau1bound[j])),
#                       runif(1, min = log(10^(-10)), max = log(tau2bound[j])), runif(1, min=-7, max=7))
#
#             if(init[2] < init[1]){
#                 init <- init[c(2,1,3)]
#             }
#
#             opt <- optim(par = init, fn = func, method = "Nelder-Mead", control = list(trace = F), k=j, m=l)
#
#             # Get back optimal values
#             par <- exp(opt$par)
#             est1 <- min(par[1:2])
#             est2 <- max(par[1:2])
#             est3 <- par[3]
#
#             obj <- opt$value
#             c(est1,est2,est3,obj)
#         }
#     }
# }
#
# stopCluster(cluster)
# rownames(like.ests) <- c("(AB)C","(BC)A","(AC)B","(BA)C","(CB)A","(AC)B")
#
# write.csv(like.ests, file = paste0("branchlikeests",reps,".csv"))
#
#
