#####################################################
# Simulation 1a 
# Adaptive sequential trial with multiple biomarkers
# collected at t=4 time points and 2 different assays
#####################################################

library(devtools)
library(data.table)

logit <- qlogis
expit <- plogis

#Necessary functions:
runifdisc<-function(n, min=0, max=1) sample(min:max, n, replace=T)
rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))

Qbar0 <- function(A, W, S) {
  
  #Baseline Covariates:
  W1 <- W[, 1]
  W2 <- W[, 2]
  
  #Surrogates:
  S1a <- S[,1]
  S2a <- S[,2]
  S3a <- S[,3]
  S4a <- S[,4]
  
  S1b <- S[,5]
  S2b <- S[,6]
  S3b <- S[,7]
  S4b <- S[,8]

  #Should there be dependence on surrogates?
  #What if for certain covariates one surrogate is better than the other?
  #Qbar <- ifelse(W2 == 0, rexpit(-2*A + log(W1, base = 10)), rexpit(-2*A + 2*log(W1, base = 10)))
  Qbar <- plogis(2*A - 1.5*log(W1, base = 10))
  return(Qbar)
}

gen_data <- function(n = 1000) {
  
  #Sample age uniformly with equal probability:
  W1 <- runifdisc(n, min = 2, max = 14)
  
  #Gender:
  W2 <- rbinom(n, 1, 0.5)
  
  W <- cbind.data.frame(W1,W2)
  
  #Randomized treatment (vaccine:placebo):
  A <- rbinom(n, 1, 0.5)
  
  #Seropositivity 
  #0 is "testing positive"
  #1 is "testing negative"
  #We want to maximize,aka as many testing negative (1)
  #More likely to test positive with age and if male (just an example)
  S1a <- rexpit(0.8*A - 0.5*log(W1, base = 10))
  S2a <- rexpit(1.2*A - 0.5*log(W1, base = 10))
  S3a <- rexpit(1.5*A - 0.5*log(W1, base = 10))
  S4a <- rexpit(1.8*A - 0.5*log(W1, base = 10))
  
  S1b <- rexpit(0.8*A - 0.7*log(W1, base = 10) - 0.5*W2)
  S2b <- rexpit(1.2*A - 0.7*log(W1, base = 10) - 0.5*W2)
  S3b <- rexpit(1.5*A - 0.7*log(W1, base = 10) - 0.5*W2)
  S4b <- rexpit(1.8*A - 0.7*log(W1, base = 10) - 0.5*W2)
  
  S <-cbind.data.frame(S1a,S2a,S3a,S4a,S1b,S2b,S3b,S4b)
  
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W, S))
  Y0 <- as.numeric(u < Qbar0(0, W, S))
  Y1 <- as.numeric(u < Qbar0(1, W, S))
  d0 <- as.numeric(Qbar0(1, W, S) > Qbar0(0, W, S))
  Yd0 <- as.numeric(u < Qbar0(d0, W, S))
  data.frame(W, A, S, Y, Y0, Y1, Yd0, d0, blip = Qbar0(1, W, S) - Qbar0(0, W, S))
}

set.seed(11)
data_full <- gen_data(1000000)
data <- gen_data(1000)
data<-data[,1:12]
#True Psi is 0.671368 (based on Y)
psi<-mean(data_full$Yd0)
mean(data_full$Y1)
mean(data_full$Y0)
rm(expit,gen_data,logit,Qbar0,rexpit,runifdisc)

devtools::use_data(data)

