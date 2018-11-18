context("Test finding the optimal surrogate")

library(sl3)
library(data.table)
library(tmle3adapt)
library(tmle3)
library(tmle3mopttx)

set.seed(1234)

data("data_rho0.2.rda")
data <- data_rho0.2[, 1:9]

# True mean under true rule:
psi <- mean(data_full_rho0.2$Yd0)

# True mean under the data-adaptive rule:

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
  learners = list(glmnet_0.8, lrn2, lrn1, lrn3),
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
  learners = list(xgboost_50, xgboost_100, glmnet_0.2, glmnet_0.8, lrn2),
  metalearner = Lrnr_nnls$new()
)

learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner, S = s_learner)

#####################################################
# Adaptive Sequential Trial with actual outcome
#####################################################

# Define spec:
tmle_spec_adapt <- tmle3_adapt(
  surrogate = TRUE,
  S = c("S1", "S2", "S3", "S4"),
  V = c("W1", "W2"),
  learners = learner_list,
  param = "opt",
  training_size = 200,
  test_size = 20,
  mini_batch = 25,
  Gexploit = 0.1,
  Gexplore = 0.01
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

tmle_params <- tmle_spec_adapt$make_params(tmle_task, targeted_likelihood)
fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
fit

##################################
# Learn the optimal surrogate
##################################

# Define spec:
tmle_spec <- tmle3_surrogate(
  S = c("S1", "S2", "S3", "S4"),
  V = c("W1", "W2"),
  learners = learner_list,
  param = "opt",
  training_size = 200,
  test_size = 20,
  mini_batch = 50
)

# Define nodes:
node_list <- list(W = c("W1", "W2", "W3"), A = "A", Y = "Y")

# Define data:
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# Define likelihood:
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
data <- tmle_spec$make_params(tmle_task, initial_likelihood)

#####################################################
# Adaptive Sequential Trial with surrogate outcome
#####################################################

# TO DO: Resolve problem with negative prediction :/

# Define spec:
tmle_spec_adapt <- tmle3_adapt(
  surrogate = TRUE,
  S = c("S1", "S2", "S3", "S4"),
  V = c("W1", "W2"),
  learners = learner_list,
  param = "opt",
  training_size = 200,
  test_size = 20,
  mini_batch = 25,
  Gexploit = 0.1,
  Gexplore = 0.01
)

# Define nodes:
node_list <- list(W = c("W1", "W2", "W3"), A = "A", Y = "Y")

# Define data (note that we use the estimated optimal surrogate here!)
# CV is defined as the rolling window, not actual Super Learner
tmle_task <- tmle_spec_adapt$make_tmle_task(data, node_list)

# Define likelihood:
initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)
summary(initial_likelihood$get_likelihood(tmle_task, node = "Y"))

updater <- tmle_spec_adapt$make_updater()
targeted_likelihood <- tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)

self <- tmle_spec_adapt
private <- tmle_spec_adapt$.__enclos_env__$private
likelihood <- targeted_likelihood

tmle_params <- tmle_spec_adapt$make_params(tmle_task, targeted_likelihood)
fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
fit
