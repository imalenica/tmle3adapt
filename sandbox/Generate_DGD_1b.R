#####################################################
# Simulation 1b
# Adaptive sequential trial with multiple biomarkers
# collected at t=2 time points and 2 different assays
#####################################################
library(bindata)

gen_data <- function(n=1000, rho=0.2, p=3, s=4) {
  
  Qbar0 <- function(A, W) {
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    #When W3>x, optimal A=1
    #When W3=<x, optimal A=0
    Qbar <- data.frame(ifelse(W3 > 1, plogis(2 * A), plogis(-2*A)))
    return(Qbar)
  }
  
  #Generate baseline covariates:
  W <- matrix(rnorm(n * p), nrow = n)
  colnames(W) <- paste("W", seq_len(p), sep = "")
  
  #Generate initial A (first round of the trial):
  A <- rbinom(n, 1, 0.5)
  
  Sigma <- matrix(rho, nrow=(s+1), ncol=(s+1)) + diag((s+1))*(1-rho)
  
  #TO DO: generalize this part
  Q <- Qbar0(A,W)

  x<-t(apply(Q, 1, function(i) rmvbin(1, margprob = c(rep(i,5)), bincorr = Sigma)))
  
  S <- x[,2:5]
  colnames(S) <- paste("S", seq_len(s), sep = "")
  
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W))
  Y0 <- as.numeric(u < Qbar0(0, W))
  Y1 <- as.numeric(u < Qbar0(1, W))
  d0 <- as.numeric(Qbar0(1, W) > Qbar0(0, W))
  Yd0 <- as.numeric(u < Qbar0(d0, W))
  blip <- data.frame(Qbar0(1, W) - Qbar0(0, W))
  colnames(blip)<-"blip"
  
  return(data.frame(W, A=A, S, Y, Y0, Y1, Yd0, d0, blip))
}

#set.seed(11)
#data_full_rho0.2 <- gen_data(1000000, rho=0.2)
#data_rho0.2 <- gen_data(1000, rho=0.2)

#data_full_rho0.8 <- gen_data(1000000, rho=0.8)
#data_rho0.8 <- gen_data(1000, rho=0.8)

gen_data_adapt <- function(n=1000, rho=0.2, p=3, s=4, Gstar, W) {
  
  Qbar0 <- function(A, W) {
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    Qbar <- data.frame(ifelse(W3 > 1, plogis(2 * A), plogis(-2*A)))
    return(Qbar)
  }
  
  A <- rbinom(n, 1, Gstar)
  
  Sigma <- matrix(rho, nrow=(s+1), ncol=(s+1)) + diag((s+1))*(1-rho)
  
  #TO DO: generalize this part
  Q <- Qbar0(A,W)
  
  x<-t(apply(Q, 1, function(i) rmvbin(1, margprob = c(rep(i,5)), bincorr = Sigma)))
  
  S <- x[,2:5]
  colnames(S) <- paste("S", seq_len(s), sep = "")
  
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W))
  Y0 <- as.numeric(u < Qbar0(0, W))
  Y1 <- as.numeric(u < Qbar0(1, W))
  d0 <- as.numeric(Qbar0(1, W) > Qbar0(0, W))
  Yd0 <- as.numeric(u < Qbar0(d0, W))
  blip <- data.frame(Qbar0(1, W) - Qbar0(0, W))
  colnames(blip)<-"blip"
  
  return(data.frame(W, A=A, S, Y))
}

gen_data_adapt_truth <- function(n=1000, rho=0.2, p=3, s=4, Gstar, W) {
  
  Qbar0 <- function(A, W) {
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    Qbar <- data.frame(ifelse(W3 > 1, plogis(2 * A), plogis(-2*A)))
    return(Qbar)
  }
  
  A <- rbinom(n, 1, Gstar)
  
  Sigma <- matrix(rho, nrow=(s+1), ncol=(s+1)) + diag((s+1))*(1-rho)
  
  #TO DO: generalize this part
  Q <- Qbar0(A,W)
  
  x<-t(apply(Q, 1, function(i) rmvbin(1, margprob = c(rep(i,5)), bincorr = Sigma)))
  
  S <- x[,2:5]
  colnames(S) <- paste("S", seq_len(s), sep = "")
  
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W))
  Y0 <- as.numeric(u < Qbar0(0, W))
  Y1 <- as.numeric(u < Qbar0(1, W))
  d0 <- as.numeric(Qbar0(1, W) > Qbar0(0, W))
  Yd0 <- as.numeric(u < Qbar0(d0, W))
  blip <- data.frame(Qbar0(1, W) - Qbar0(0, W))
  colnames(blip)<-"blip"
  
  return(data.frame(gen_data_adapt <- function(n=1000, rho=0.2, p=3, s=4, Gstar, W) {
    
    Qbar0 <- function(A, W) {
      W1 <- W[, 1]
      W2 <- W[, 2]
      W3 <- W[, 3]
      Qbar <- data.frame(ifelse(W3 > 1, plogis(2 * A), plogis(-2*A)))
      return(Qbar)
    }
    
    A <- rbinom(n, 1, Gstar)
    
    Sigma <- matrix(rho, nrow=(s+1), ncol=(s+1)) + diag((s+1))*(1-rho)
    
    #TO DO: generalize this part
    Q <- Qbar0(A,W)
    
    x<-t(apply(Q, 1, function(i) rmvbin(1, margprob = c(rep(i,5)), bincorr = Sigma)))
    
    S <- x[,2:5]
    colnames(S) <- paste("S", seq_len(s), sep = "")
    
    u <- runif(n)
    Y <- as.numeric(u < Qbar0(A, W))
    Y0 <- as.numeric(u < Qbar0(0, W))
    Y1 <- as.numeric(u < Qbar0(1, W))
    d0 <- as.numeric(Qbar0(1, W) > Qbar0(0, W))
    Yd0 <- as.numeric(u < Qbar0(d0, W))
    blip <- data.frame(Qbar0(1, W) - Qbar0(0, W))
    colnames(blip)<-"blip"
    
    return(data.frame(W, A=A, S, Y, Y0, Y1, d0, Yd0))
  }))
}
