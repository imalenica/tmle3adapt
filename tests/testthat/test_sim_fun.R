context("Test the simulation function")

library(testthat)
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
learners <- learner_list
surrogate <- TRUE
S <- c("S1", "S2", "S3", "S4")
W <- c("W1", "W2", "W3")
A <- "A"
Y <- "Y"
V <- NULL
param <- "opt"
training_size <- 100
test_size <- 20
mini_batch <- 25
Gexploit <- 0.1
Gexplore <- 0.01
n_max <- 1600
by <- 200
n <- 1000
source(here("Sandbox/Generate_DGD_1b.R"))

est <- tmle3_sadapt(
  surrogate = surrogate, S = S, W = W, V = NULL, A = A, Y = Y,
  learners = learner_list, data = data,
  param = param, training_size = training_size, test_size = test_size, mini_batch = mini_batch,
  Gexploit = Gexploit, Gexplore = Gexplore,
  n_max = n_max, by = by, n = n,
  gen_data = gen_data, gen_data_adapt = gen_data_adapt
)
