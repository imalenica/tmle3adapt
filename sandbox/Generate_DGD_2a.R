#####################################################
# Simulation 2a (parameter: ATE)
#
# Test ATE parameter
# 
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
    
    #Qbar <- plogis(2*A - 1.4*W2 + S3a)
    Qbar <- rexpit(2*A - S3a)
    return(Qbar)
  }
  
  W1 <- rnorm(n = n, mean = 2, sd = 0.2)
  W2 <- rbinom(n, 1, 0.7)
  W3 <- rbinom(n, 1, 0.3)
  
  W <- cbind.data.frame(W1,W2,W3)
  
  A <- rbinom(n, 1, 0.5)
  
  #S1a: not predictive of outcome at all, random
  S1a <- rbinom(n=n, size = 1, prob = 0.3)
  S2a <- rexpit(A - 0.3*W2)
  S3a <- rexpit(A - 0.8*W2 + 1.5*W3)
  
  S <-cbind.data.frame(S1a,S2a,S3a)
  
  #u <- runif(n)
  #Y <- as.numeric(u < Qbar0(A, W, S))
  #Y0 <- as.numeric(u < Qbar0(0, W, S))
  #Y1 <- as.numeric(u < Qbar0(1, W, S))
  Y <- Qbar0(A,W,S)
  Y0 <- Qbar0(0,W,S)
  Y1 <- Qbar0(1,W,S)
  
  data.frame(W, A, S, Y, Y0, Y1)
}

set.seed(11)
data_full <- gen_data(1000000)
data <- gen_data(1000)
data<-data[,1:8]
#True Psi for ATE is 0.426038 
psi<-mean(data_full$Y1) - mean(data_full$Y0)
mean(data_full$Y1)
mean(data_full$Y0)
devtools::use_data(data, psi, internal = TRUE, overwrite = TRUE)

