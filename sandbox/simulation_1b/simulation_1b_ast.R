########################################
# Simulation 1b: True Y
########################################

library(here)
library(data.table)
library(tidyverse)
library(sl3)
library(tmle3adapt)
library(tmle3)
library(tmle3mopttx)
library(origami)
library(devtools)

load_all()

load(here("data_rho0.2.rda"))
data <- data_rho0.2[, 1:9]

#True mean under true rule:
psi <- mean(data_full_rho0.2$Yd0)

# Define sl3 library and metalearners:
xgboost_50 <- Lrnr_xgboost$new(nrounds = 50)
xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
xgboost_300 <- Lrnr_xgboost$new(nrounds = 300)
xgboost_500 <- Lrnr_xgboost$new(nrounds = 500)
glmnet_0.2 <- Lrnr_glmnet$new(alpha = 0.2, lambda = 300)
glmnet_0.6 <- Lrnr_glmnet$new(alpha = 0.6, lambda = 300)
glmnet_0.8 <- Lrnr_glmnet$new(alpha = 0.8, lambda = 300)

lrn1 <- Lrnr_mean$new()
lrn2 <- Lrnr_glm_fast$new()
lrn3 <- Lrnr_hal9001$new(degrees = 3, n_folds = 3)

Q_learner <- Lrnr_sl$new(
  learners = list(xgboost_300, glmnet_0.8,lrn2, lrn3),
  #learners = list(xgboost_300, xgboost_500, glmnet_0.8, lrn1, lrn2, lrn3),
  metalearner = Lrnr_nnls$new()
)

g_learner <- Lrnr_sl$new(
  learners = list(xgboost_50, xgboost_100, glmnet_0.2, glmnet_0.8, lrn2),
  metalearner = Lrnr_nnls$new()
)

b_learner <- Lrnr_sl$new(
  learners = list(xgboost_50, xgboost_100, glmnet_0.2, glmnet_0.8, lrn2),
  metalearner = Lrnr_nnls$new()
)

s_learner <- Lrnr_sl$new(
  learners = list(xgboost_300, glmnet_0.8,lrn2, lrn3),
  #learners = list(xgboost_50, xgboost_100, xgboost_300, xgboost_500,glmnet_0.2, glmnet_0.8, lrn2, lrn3),
  metalearner = Lrnr_nnls$new()
)

learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner, S = s_learner)
learners = learner_list
surrogate = TRUE
S = c("S1", "S2", "S3", "S4")
W = c("W1","W2","W3")
A = "A"
Y = "Y"
V = NULL
param = "opt"
training_size = 800
test_size = 20
mini_batch = 25
Gexploit = 0.1
Gexplore = 0.01
n_max = 1600
by = 200
n = 1000
rho=0.2
rule_outcome="Y"
opt_surrogate="TMLE"
source(here("sandbox/Generate_DGD_1b.R"))
source(here("sandbox/simulation_1b/tmle3_adapt_sim1b.R"))

set.seed(1234)
MC = 500
est_d0 <- list()
se <- data.frame(matrix(NA, nrow = MC, ncol = 4))
cov <- data.frame(matrix(NA, nrow = MC, ncol = 4))
colnames(cov)<-c("Initial", "Update 1","Update 2",  "Update 3")
colnames(se)<-c("Initial", "Update 1","Update 2",  "Update 3")

# Monte Carlo iterations
for(i in 1:MC) {

  res<-tmle3_sadapt_sim1b(surrogate = surrogate,
                          S = S,
                          W = W,
                          V = V,
                          A = A,
                          Y = Y,
                          learners = learners,
                          data=data,
                          param = param,
                          training_size = training_size,
                          test_size = test_size,
                          mini_batch = mini_batch,
                          Gexploit = Gexploit,
                          Gexplore = Gexplore,
                          n_max = n_max,
                          by = by,
                          n = nrow(data),
                          gen_data = gen_data,
                          gen_data_adapt = gen_data_adapt,
                          gen_data_adapt_truth = gen_data_adapt_truth,
                          rho = rho,
                          rule_outcome=rule_outcome,
                          opt_surrogate=opt_surrogate)
  
  fin<-res$summary
  
  #Check coverage at each step:
  cov[i,]<-apply(fin, 1, function(x) ifelse((x[9]<=psi & x[10]>=psi),1,0) )
  
  #Get standard error:
  se[i,] <-  t(fin$se)
  
  #Save estimates:
  est_d0<-c(est_d0,list(fin))

}

##################
# n=500

learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner, S = s_learner)
learners = learner_list
surrogate = FALSE
S = c("S1", "S2", "S3", "S4") 
W = c("W1","W2","W3")
A = "A"
Y = "Y"
V = NULL
param = "opt"
training_size = 300
test_size = 20
mini_batch = 25
Gexploit = 0.1
Gexplore = 0.01
n_max = 1100
by = 200
n = 500
rho=0.2
data<-data[1:500,]

MC = 500
est_d0_n500 <- list()
se_n500 <- data.frame(matrix(NA, nrow = MC, ncol = 5))
cov_n500 <- data.frame(matrix(NA, nrow = MC, ncol = 5))
colnames(cov_n500)<-c("Initial", "Update 1","Update 2",  "Update 3",  "Update 4")
colnames(se_n500)<-c("Initial", "Update 1","Update 2",  "Update 3",  "Update 4")

# Monte Carlo iterations
for(i in 1:MC) {
  
  res<-tmle3_sadapt_sim1b(surrogate = surrogate,
                          S = S,
                          W = W,
                          V = V,
                          A = A,
                          Y = Y,
                          learners = learners,
                          data=data,
                          param = param,
                          training_size = training_size,
                          test_size = test_size,
                          mini_batch = mini_batch,
                          Gexploit = Gexploit,
                          Gexplore = Gexplore,
                          n_max = n_max,
                          by = by,
                          n = nrow(data),
                          gen_data = gen_data,
                          gen_data_adapt = gen_data_adapt,
                          gen_data_adapt_truth = gen_data_adapt_truth,
                          rho = rho)
  
  fin<-res$summary
  
  #Check coverage at each step:
  cov_n500[i,]<-apply(fin, 1, function(x) ifelse((x[9]<=psi & x[10]>=psi),1,0) )
  
  #Get standard error:
  se_n500[i,] <-  t(fin$se)
  
  #Save estimates:
  est_d0_n500<-c(est_d0_n500,list(fin))
  
}

