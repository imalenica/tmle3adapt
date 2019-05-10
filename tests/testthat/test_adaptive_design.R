context("Test usual adaptive sequential design")

library(testthat)
library(sl3)
library(data.table)
library(tmle3adapt)
library(tmle3)
library(tmle3mopttx)
library(R6)
library(devtools)
library(here)
load_all()
source(here("sandbox/Generate_DGD_2a.R"))

set.seed(1234)
data <- gen_data(n=200)
data <- data[,1:5]

# True mean under true rule:
psi <- mean(data_full$Yd0)

### True mean under the data-adaptive rule:

# Define sl3 library and metalearners:
xgboost_50 <- Lrnr_xgboost$new(nrounds = 50)
xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
xgboost_300 <- Lrnr_xgboost$new(nrounds = 300)
glmnet_0.2 <- Lrnr_glmnet$new(alpha = 0.2, lambda = 300)
glmnet_0.6 <- Lrnr_glmnet$new(alpha = 0.6, lambda = 300)
glmnet_0.8 <- Lrnr_glmnet$new(alpha = 0.8, lambda = 300)
lrn1 <- Lrnr_mean$new()
lrn2 <- Lrnr_glm_fast$new()
lrn3 <- Lrnr_hal9001$new()

Q_learner <- Lrnr_sl$new(
  # learners = list(xgboost_50,xgboost_100,glmnet_0.2,glmnet_0.8,lrn2),
  learners = list(xgboost_100, lrn2, lrn1),
  metalearner = Lrnr_nnls$new()
)

g_learner <- Lrnr_sl$new(
  # learners = list(xgboost_50, xgboost_100, glmnet_0.2, glmnet_0.8, lrn2),
  learners = list(xgboost_100, lrn2, lrn1),
  metalearner = Lrnr_nnls$new()
)

b_learner <- Lrnr_sl$new(
  # learners = list(xgboost_50, xgboost_100, glmnet_0.2, glmnet_0.8, lrn2),
  learners = list(xgboost_100, lrn2, lrn1),
  metalearner = Lrnr_nnls$new()
)

s_learner <- Lrnr_sl$new(
  # learners = list(xgboost_50, xgboost_100, glmnet_0.2, glmnet_0.8, lrn2),
  learners = list(xgboost_100, lrn2, lrn1),
  metalearner = Lrnr_nnls$new()
)

learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner, S = s_learner)

#####################################################
# Adaptive Sequential Trial with actual outcome
#####################################################

### First round:

# Define spec:
tmle_spec_adapt <- tmle3_adapt(
  S = c("S1", "S2", "S3", "S4"),
  V = c("W1", "W2"),
  learners = learner_list,
  param = "opt",
  training_size = 100,
  test_size = 20,
  mini_batch = 25,
  grid_GG = FALSE,
  Gexploit = 0.1,
  Gexplore = 0.08
)

# Define nodes:
node_list <- list(W = c("W1", "W2", "W3"), A = "A", Y = "Y")

# Define data (note that we use the estimated optimal surrogate here!)
# CV is defined as the rolling window, not actual Super Learner
tmle_task <- tmle_spec_adapt$make_tmle_task(data, node_list)

# Define likelihood:
initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)

updater <- tmle_spec_adapt$make_updater()
targeted_likelihood <- tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)

tmle_params <- tmle_spec_adapt$make_params(tmle_task, likelihood=targeted_likelihood)
fit <- fit_tmle3(tmle_task, likelihood=targeted_likelihood, tmle_params=tmle_params, updater)
fit

#summary(initial_likelihood$get_likelihood(tmle_task, node = "Y"))
#summary(targeted_likelihood$get_likelihood(tmle_task, node = "Y"))

# Needs to be counterfactual task!
blip_task <- tmle_spec_adapt$get_blip_cf(tmle_task)
rule1 <- tmle_spec_adapt$get_Gstar(blip_task, initial_likelihood)

# Now use this Q estimate in initial_likelihood with a new task-
# but now our task is defined by the new obtained covariates
data_new <- gen_data(n = 200)
W <- cbind.data.frame(W1=data_new$W1,W2=data_new$W2,W3=data_new$W3)
tmle_task_new <- tmle_spec_adapt$make_tmle_task(data_new, node_list)
blip_task <- tmle_spec_adapt$get_blip_cf(tmle_task_new)
dn <- tmle_spec_adapt$get_Gstar(blip_task, initial_likelihood)

# New data with the collected W, but the randomization probability is changed
# based on the initial fit of Q.
data_targeted <- lapply(dn, function(x){
  gen_data_adapt(n = 200, Gstar = x, W)
})

# Combine data:
data <- lapply(data_targeted, function(x){
  rbind.data.frame(data, x)
})
  
########
# Repeat!
########
tmle_tasks <- lapply(data, function(x){
  tmle_spec_adapt$make_tmle_task(x, node_list)
})
  
#Define likelihood:
initial_likelihoods <- lapply(tmle_tasks, function(tmle_task){
  tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)
})
                             
updaters <- lapply(seq(length(tmle_tasks)), function(x){
  tmle_spec_adapt$make_updater()
})

targeted_likelihoods <- lapply(seq(length(tmle_tasks)), function(x){
  tmle_spec_adapt$make_targeted_likelihood(initial_likelihoods[[x]], updaters[[x]])
})
  
tmle_params <- lapply(seq(length(tmle_tasks)), function(x){
  tmle_spec_adapt$make_params(tmle_tasks[[x]], targeted_likelihoods[[x]])
})
  
fits <- lapply(seq(length(tmle_tasks)), function(x){
  fit_tmle3(tmle_tasks[[x]], targeted_likelihoods[[x]], tmle_params[[x]], updaters[[x]])
})
  
fits

#Pick the one that maximizes the outcome?
#Has the lowest variance?
#

summary(initial_likelihood$get_likelihood(tmle_task, node = "Y"))
summary(targeted_likelihood$get_likelihood(tmle_task, node = "Y"))

# Needs to be counterfactual task!
blip_task <- tmle_spec_adapt$get_blip_cf(tmle_task)
rule1 <- tmle_spec_adapt$get_Gstar(blip_task, initial_likelihood)
