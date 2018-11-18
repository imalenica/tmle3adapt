###############################################
# Simulation 1a: True Y and data-adaptive Y
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

load(here("Data/data_sim1a.rda"))

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
S = c("S1a","S2a","S3a")
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
n_max = 1400
by = 200
n = 1000
rule_outcome="Y"
opt_surrogate="TMLE"
data_adaptive=TRUE
source(here("sandbox/Generate_DGD_1a.R"))
source(here("sandbox/tmle3_adapt_sim_osc.R"))

set.seed(1234)
MC = 500

res_all<-list()
mse<-data.frame(matrix(NA, nrow = MC, ncol = 6))
colnames(mse) <- c("Initial TMLE", "Initial SL", "Update 1 TMLE", "Update 1 SL", "Update 2 TMLE", "Update 2 SL")
cov <- data.frame(matrix(NA, nrow = MC, ncol = 6))
colnames(cov)<-c("Initial TMLE","Initial SL","Update 1 TMLE", "Update 1 SL",
                 "Update 2 TMLE", "Update 2 SL")
se_tmle <- data.frame(matrix(NA, nrow = MC, ncol = 3))
se_sl <- data.frame(matrix(NA, nrow = MC, ncol = 3))
colnames(se_tmle) <- c("Initial", "Update 1", "Update 2")
colnames(se_sl) <- c("Initial", "Update 1", "Update 2")

# Monte Carlo iterations
for(i in 1:MC) {
  
  res<-tmle3_adapt_sim_osc(surrogate = surrogate,
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
  
  summary_tmle <- res$summary_tmle
  summary_sl <- res$summary_sl
  
  #Save all results:
  res_all <- c(res_all, list(res))
  
  #MSE:
  mse[i,] <- res$mse
  
  #Coverage of data adative parameter (w.r.t. true Y):
  psi_dn <- data.frame(dn=res$psi_dn)
  
  cov[i, c(1,3,5) ] <- ifelse((summary_tmle$lower_transformed<=psi_dn & 
                                 summary_tmle$upper_transformed>=psi_dn),1,0)
  cov[i, c(2,4,6) ] <- ifelse((summary_sl$lower_transformed<=psi_dn & 
                                 summary_sl$upper_transformed>=psi_dn),1,0)
  
  #Get standard error:
  se_tmle[i,] <-  t(summary_tmle$se)
  se_sl[i,] <-  t(summary_sl$se)
  
}
