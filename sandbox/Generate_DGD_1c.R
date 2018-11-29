#####################################################
# Simulation 1c
# 
# Test Continuous Outcome
#
# Adaptive sequential trial with multiple biomarkers
# collected at t=3 time points 
# First surrogate is not predictive of the outcome
#####################################################

library(devtools)
library(data.table)

gen_data <- function(n = 1000) {
  
  logit <- qlogis
  expit <- plogis
  
  #Necessary functions:
  runifdisc<-function(n, min=0, max=1) {sample(min:max, n, replace=T)}
  rexpit <- function(x) {rbinom(n=length(x), size=1, prob=plogis(x))}
  
  Qbar0 <- function(A, W, S) {
    
    #Baseline Covariates:
    W1 <- W[, 1]
    W2 <- W[, 2]
    
    #Surrogates:
    S1a <- S[,1]
    S2a <- S[,2]
    S3a <- S[,3]
    
    Qbar <- ifelse(W3 == 1, 2*A + 2*S3a, -2*A + 2*S3a)
    return(Qbar)
  }
  
  #Sample age uniformly with equal probability:
  W1 <- runifdisc(n, min = 2, max = 14)
  
  #Gender:
  W2 <- rbinom(n, 1, 0.5)
  W3 <- rbinom(n, 1, 0.6)
  
  W <- cbind.data.frame(W1,W2,W3)
  
  #Randomized treatment (vaccine:placebo):
  A <- rbinom(n, 1, 0.5)
  
  S1a<-rnorm(n, mean=0, sd=0.2)
  S2a<-ifelse(W3==1, rnorm(n, mean=2*A*W2 + 0.2*log(W1, base = 10), sd=0.2),
              rnorm(n, mean=-2*A*W2 + 0.2*log(W1, base = 10), sd=0.4))
  S3a<-ifelse(W3==1, rnorm(n, mean=3*A*W2 + 0.2*log(W1, base = 10), sd=0.2),  
              rnorm(n, mean=-3*A*W2 + 0.2*log(W1, base = 10), sd=0.4))

  S <-cbind.data.frame(S1a,S2a,S3a)
  
  Y <- rnorm(n, mean=Qbar0(A, W, S), sd=0.4)
  Y0 <- rnorm(n, mean=Qbar0(0, W, S), sd=0.4)
  Y1 <- rnorm(n, mean=Qbar0(1, W, S), sd=0.4)
  d0 <- as.numeric(Qbar0(1, W, S) > Qbar0(0, W, S))
  Yd0 <- rnorm(n, mean=Qbar0(d0, W, S), sd=0.4)
  blip = Qbar0(1, W, S) - Qbar0(0, W, S)
  names(blip)<-"blip"
  data.frame(W, A, S, Y, Y0, Y1, Yd0, d0, blip = blip)
}

set.seed(11)
data_full <- gen_data(1000000)
data <- gen_data(1000)
data<-data[,1:8]
#True Psi is 1.8409 (based on Y)
psi<-mean(data_full$Yd0)
mean(data_full$Y1)
mean(data_full$Y0)
rm(gen_data)

#devtools::use_data(data, data_full, psi, internal = TRUE)

gen_data_adapt_truth <- function(n = 1000, Gstar, W) {
  
  logit <- qlogis
  expit <- plogis
  
  #Necessary functions:
  runifdisc<-function(n, min=0, max=1) {sample(min:max, n, replace=T)}
  rexpit <- function(x) {rbinom(n=length(x), size=1, prob=plogis(x))}
  
  Qbar0 <- function(A, W, S) {
    
    #Baseline Covariates:
    W1 <- W[, 1]
    W2 <- W[, 2]
    
    #Surrogates:
    S1a <- S[,1]
    S2a <- S[,2]
    S3a <- S[,3]
    
    Qbar <- ifelse(W3 == 1, 2*A + 2*S3a, -2*A + 2*S3a)
    return(Qbar)
  }
  
  #Randomized treatment (vaccine:placebo):
  A <- rbinom(n, 1, Gstar)
  
  S1a<-rnorm(n, mean=0, sd=0.2)
  S2a<-ifelse(W3==1, rnorm(n, mean=2*A*W2 + 0.2*log(W1, base = 10), sd=0.2),
              rnorm(n, mean=-2*A*W2 + 0.2*log(W1, base = 10), sd=0.4))
  S3a<-ifelse(W3==1, rnorm(n, mean=3*A*W2 + 0.2*log(W1, base = 10), sd=0.2),  
              rnorm(n, mean=-3*A*W2 + 0.2*log(W1, base = 10), sd=0.4))
  
  S <-cbind.data.frame(S1a,S2a,S3a)
  
  Y <- rnorm(n, mean=Qbar0(A, W, S), sd=0.4)
  Y0 <- rnorm(n, mean=Qbar0(0, W, S), sd=0.4)
  Y1 <- rnorm(n, mean=Qbar0(1, W, S), sd=0.4)
  d0 <- as.numeric(Qbar0(1, W, S) > Qbar0(0, W, S))
  Yd0 <- rnorm(n, mean=Qbar0(d0, W, S), sd=0.4)
  blip = Qbar0(1, W, S) - Qbar0(0, W, S)
  names(blip)<-"blip"
  data.frame(W, A, S, Y, Y0, Y1, Yd0, d0, blip = blip)
}

gen_data_adapt_truth <- function(n = 1000, Gstar, W) {
  
  logit <- qlogis
  expit <- plogis
  
  #Necessary functions:
  runifdisc<-function(n, min=0, max=1) {sample(min:max, n, replace=T)}
  rexpit <- function(x) {rbinom(n=length(x), size=1, prob=plogis(x))}
  
  Qbar0 <- function(A, W, S) {
    
    #Baseline Covariates:
    W1 <- W[, 1]
    W2 <- W[, 2]
    
    #Surrogates:
    S1a <- S[,1]
    S2a <- S[,2]
    S3a <- S[,3]
    
    Qbar <- ifelse(W3 == 1, 2*A + 2*S3a, -2*A + 2*S3a)
    return(Qbar)
  }
  
  #Randomized treatment (vaccine:placebo):
  A <- rbinom(n, 1, Gstar)
  
  S1a<-rnorm(n, mean=0, sd=0.2)
  S2a<-ifelse(W3==1, rnorm(n, mean=2*A*W2 + 0.2*log(W1, base = 10), sd=0.2),
              rnorm(n, mean=-2*A*W2 + 0.2*log(W1, base = 10), sd=0.4))
  S3a<-ifelse(W3==1, rnorm(n, mean=3*A*W2 + 0.2*log(W1, base = 10), sd=0.2),  
              rnorm(n, mean=-3*A*W2 + 0.2*log(W1, base = 10), sd=0.4))
  
  S <-cbind.data.frame(S1a,S2a,S3a)
  
  Y <- rnorm(n, mean=Qbar0(A, W, S), sd=0.4)

  data.frame(W, A, S, Y)
}
