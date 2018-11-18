#####################################################
# Simulation 1a 
# Adaptive sequential trial with multiple biomarkers
# collected at t=4 time points 
#####################################################

library(devtools)
library(data.table)

gen_data = function(n = 1000) {
  
  logit <- qlogis
  expit <- plogis
  
  #Necessary functions:
  runifdisc<-function(n, min=0, max=1) {sample(min:max, n, replace=T)}
  rexpit <- function(x) {rbinom(n=length(x), size=1, prob=plogis(x))}
  
  Qbar0 <- function(A, W, S) {
    
    #Baseline Covariates:
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    
    #Surrogates:
    S1a <- S[,1]
    S2a <- S[,2]
    S3a <- S[,3]

    Qbar <- data.frame(ifelse(W2 == 1, plogis(2*A + S3a), plogis(-2*A + S3a)))
    return(Qbar)
  }
  
  W1 <- rnorm(n = n, mean = 0.5, sd = 0.2)
  W2 <- rbinom(n, 1, 0.6)
  W3 <- rnorm(n = n, mean = 1, sd = 0.4)
  
  W <- cbind.data.frame(W1,W2,W3)
  
  A <- rbinom(n, 1, 0.5)
  
  #S1a: not predictive of outcome at all, random
  S1a <- rbinom(n=n, size = 1, prob = 0.3)
  S2a <- ifelse(W2 == 1, rexpit(0.2*A + 0.2*W1), rexpit(0.2*A - 0.2*W1))
  S3a <- ifelse(W2 == 1, rexpit(1.2*A + 0.7*W1), rexpit(1.2*A - 1.2*W1))
  
  S <-cbind.data.frame(S1a,S2a,S3a)
  
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W, S))
  Y0 <- as.numeric(u < Qbar0(0, W, S))
  Y1 <- as.numeric(u < Qbar0(1, W, S))
  d0 <- as.numeric(Qbar0(1, W, S) > Qbar0(0, W, S))
  Yd0 <- as.numeric(u < Qbar0(d0, W, S))
  blip = Qbar0(1, W, S) - Qbar0(0, W, S)
  names(blip)="blip"
  data.frame(W, A, S, Y, Y0, Y1, Yd0, d0, blip = blip)
}

#set.seed(11)
#data_full <- gen_data(1000000)
#data <- gen_data(1000)
#data<-data[,1:8]
#True Psi is 0.805 (based on Y)
#psi<-mean(data_full$Yd0)
#mean(data_full$Y1)
#mean(data_full$Y0)
#rm(gen_data)

#devtools::use_data(data, data_full, psi, internal = TRUE)

gen_data_adapt_truth = function(n = 1000, Gstar, W) {
  
  logit <- qlogis
  expit <- plogis
  
  #Necessary functions:
  runifdisc<-function(n, min=0, max=1) {sample(min:max, n, replace=T)}
  rexpit <- function(x) {rbinom(n=length(x), size=1, prob=plogis(x))}
  
  Qbar0 <- function(A, W, S) {
    
    #Baseline Covariates:
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    
    #Surrogates:
    S1a <- S[,1]
    S2a <- S[,2]
    S3a <- S[,3]
    
    Qbar <- data.frame(ifelse(W2 == 1, plogis(2*A + S3a), plogis(-2*A + S3a)))
    return(Qbar)
  }
  
  A <- rbinom(n, 1, Gstar)
  
  W1 <- W[, 1]
  W2 <- W[, 2]
  W3 <- W[, 3]
  
  #S1a: not predictive of outcome at all, random
  S1a <- rbinom(n=n, size = 1, prob = 0.3)
  S2a <- ifelse(W2 == 1, rexpit(0.2*A + 0.2*W1), rexpit(0.2*A - 0.2*W1))
  S3a <- ifelse(W2 == 1, rexpit(1.2*A + 0.7*W1), rexpit(1.2*A - 1.2*W1))
  
  S <-cbind.data.frame(S1a,S2a,S3a)
  
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W, S))
  Y0 <- as.numeric(u < Qbar0(0, W, S))
  Y1 <- as.numeric(u < Qbar0(1, W, S))
  d0 <- as.numeric(Qbar0(1, W, S) > Qbar0(0, W, S))
  Yd0 <- as.numeric(u < Qbar0(d0, W, S))
  blip = Qbar0(1, W, S) - Qbar0(0, W, S)
  names(blip)="blip"
  
  data.frame(W, A, S, Y, Y0, Y1, Yd0, d0, blip = blip)
}

gen_data_adapt = function(n = 1000, Gstar, W) {
  
  logit <- qlogis
  expit <- plogis
  
  #Necessary functions:
  runifdisc<-function(n, min=0, max=1) {sample(min:max, n, replace=T)}
  rexpit <- function(x) {rbinom(n=length(x), size=1, prob=plogis(x))}
  
  Qbar0 <- function(A, W, S) {
    
    #Baseline Covariates:
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    
    #Surrogates:
    S1a <- S[,1]
    S2a <- S[,2]
    S3a <- S[,3]
    
    Qbar <- data.frame(ifelse(W2 == 1, plogis(2*A + S3a), plogis(-2*A + S3a)))
    return(Qbar)
  }
  
  A <- rbinom(n, 1, Gstar)
  
  W1 <- W[, 1]
  W2 <- W[, 2]
  W3 <- W[, 3]
  
  #S1a: not predictive of outcome at all, random
  S1a <- rbinom(n=n, size = 1, prob = 0.3)
  S2a <- ifelse(W2 == 1, rexpit(0.2*A + 0.2*W1), rexpit(0.2*A - 0.2*W1))
  S3a <- ifelse(W2 == 1, rexpit(1.2*A + 0.7*W1), rexpit(1.2*A - 1.2*W1))
  
  S <-cbind.data.frame(S1a,S2a,S3a)
  
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W, S))

  data.frame(W, A, S, Y)
}


