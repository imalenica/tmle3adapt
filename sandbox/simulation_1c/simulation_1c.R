###############################################
# Simulation 1c: True Y and data-adaptive Y
###############################################

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

load(here("Data/data_sim1c.rda"))

#True mean under true rule:
psi <- mean(data_full$Yd0)

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
  learners = list(xgboost_300, lrn2, lrn3),
  #learners = list(xgboost_300, xgboost_500, glmnet_0.8, lrn1, lrn2, lrn3),
  metalearner = Lrnr_nnls$new()
)

g_learner <- Lrnr_sl$new(
  learners = list(xgboost_50, xgboost_100, lrn2),
  metalearner = Lrnr_nnls$new()
)

b_learner <- Lrnr_sl$new(
  learners = list(xgboost_50, xgboost_100, lrn2),
  metalearner = Lrnr_nnls$new()
)

s_learner <- Lrnr_sl$new(
  learners = list(xgboost_300,lrn2, lrn3),
  #learners = list(xgboost_50, xgboost_100, xgboost_300, xgboost_500,glmnet_0.2, glmnet_0.8, lrn2, lrn3),
  metalearner = Lrnr_nnls$new()
)

learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner, S = s_learner)
learners = learner_list
surrogate = TRUE
S = c("S1a", "S2a", "S3a")
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
rule_outcome="Y"
opt_surrogate="TMLE"
data_adaptive=FALSE
source(here("sandbox/Generate_DGD_1c.R"))
source(here("sandbox/tmle3_adapt_sim.R"))

set.seed(1234)
MC = 500
est_d0 <- list()
se <- data.frame(matrix(NA, nrow = MC, ncol = 4))
cov <- data.frame(matrix(NA, nrow = MC, ncol = 4))
colnames(cov)<-c("Initial", "Update 1","Update 2",  "Update 3")
colnames(se)<-c("Initial", "Update 1","Update 2",  "Update 3")

# Monte Carlo iterations
for(i in 1:MC) {
  
  res<-tmle3_adapt_sim(surrogate = surrogate,
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
                          rule_outcome=rule_outcome,
                          opt_surrogate=opt_surrogate,
                          data_adaptive=data_adaptive)
  
  fin<-res$summary
  
  #Check coverage at each step:
  cov[i,]<-apply(fin, 1, function(x) ifelse((x[9]<=psi & x[10]>=psi),1,0) )
  
  #Get standard error:
  se[i,] <-  t(fin$se)
  
  #Save estimates:
  est_d0<-c(est_d0,list(fin))
  
}
