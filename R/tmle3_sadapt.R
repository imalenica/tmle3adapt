#' Learns the Optimal Surrogate based on the observed final outcome Y, for the
#' Mean Under the Optimal Individualized Rule and Average Treatment Effect
#' target parameters.
#'
#' O=(W,A,S,Y)
#' W=Covariates
#' A=Treatment (binary)
#' S=Potential Surrogates
#' Y=Outcome (binary or bounded continuous)
#'
#' @param surrogate If \code{TRUE}, adaptively learns the target parameter with respect to the
#' estimated optimal surrogate.
#' @param S Covariates to consider for the Optimal Surrogate estimation.
#' @param W Baseline covariates used in the trial.
#' @param V Covariates the rule depends on. Not used currently.
#' @param A Treatment used on the trial.
#' @param Y Outcome node used in the trial.
#' @param learners List of learners used for Q,g,S and B.
#' @param param Target parameter. Current implementation supports Mean under the Optimal Individualized
#' Treatment (opt) and Average Treatment Effect (are)
#' @param training_size Size of the initial training set. Necessary part of online Super Learner.
#' @param test_size Size of the test set. Necessary part of online Super Learner.
#' @param mini_batch Size of the increase in the initial training size, added per each iteration of the
#' online Super Learner.
#' @param Gexploit Sequence t_n. Default is 0.1.
#' @param Gexplore Sequence e_n. Default is 0.05.
#' @param n_max Maximum sample size after which the trial ends.
#' @param by Size of the iterative sample size increase in the sequential trial.
#' @param n Initial sample size.
#' @param gen_data Data generating distribution of the process. Used for simulations.
#' @param gen_data_adapt Data generating distribution of the process with used specifed rule. Used for simulations.
#'
#' @export
#'

# TO DO: How to generalize this?
# Right now, tailored for the simulation
tmle3_sadapt <- function(surrogate = TRUE, S = S, W = c("W1", "W2", "W3"), V = NULL, A = "A", Y = "Y",
                         learners = learner_list, data,
                         param = "opt", training_size = 100, test_size = 20, mini_batch = 25,
                         Gexploit = 0.1, Gexplore = 0.01,
                         n_max = 1600, by = 200, n = nrow(data),
                         gen_data = gen_data, gen_data_adapt = gen_data_adapt) {
  if (surrogate) {

    # Define spec:
    tmle_spec <- tmle3_surrogate(
      S = S,
      V = V,
      learners = learner_list,
      param = param,
      training_size = training_size, 
      test_size = test_size,
      mini_batch = mini_batch
    )

    # Define nodes:
    node_list <- list(W = W, A = A, Y = Y)

    # Define data:
    tmle_task <- tmle_spec$make_tmle_task(data, node_list)

    # Define likelihood:
    initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
    data <- tmle_spec$make_params(tmle_task, initial_likelihood)

    # Define spec:
    tmle_spec_adapt <- tmle3_adapt(
      S = S,
      V = V,
      learners = learner_list,
      param = param,
      training_size = training_size,
      test_size = test_size,
      mini_batch = mini_batch,
      Gexploit = Gexploit,
      Gexplore = Gexplore
    )

    tmle_task <- tmle_spec_adapt$make_tmle_task(data = data, node_list)

    initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)
    updater <- tmle_spec_adapt$make_updater()
    targeted_likelihood <- tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)
    tmle_params <- tmle_spec_adapt$make_params(tmle_task, targeted_likelihood)

    while (n <= n_max) {
      data <- tmle_spec_adapt$new_Gstar(gen_data, gen_data_adapt, by, old_data, node_list, initial_likelihood)
      tmle_task <- tmle_spec_adapt$make_tmle_task(data, node_list)

      # Define a new likelihood:
      initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)
      updater <- tmle_spec_adapt$make_updater()
      targeted_likelihood <- tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)

      tmle_params <- c(tmle_params, tmle_spec_adapt$make_params(tmle_task, targeted_likelihood))

      n <- n + by
    }

    fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
  } else {

    # Define spec:
    tmle_spec_adapt <- tmle3_adapt(
      S = S,
      V = V,
      learners = learner_list,
      param = param,
      training_size = training_size,
      test_size = test_size,
      mini_batch = mini_batch,
      Gexploit = Gexploit,
      Gexplore = Gexplore
    )

    # Define nodes and task:
    node_list <- list(W = W, A = A, Y = Y)

    tmle_task <- tmle_spec_adapt$make_tmle_task(data = data, node_list)

    initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)
    updater <- tmle_spec_adapt$make_updater()
    targeted_likelihood <- tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)
    tmle_params <- tmle_spec_adapt$make_params(tmle_task, targeted_likelihood)

    while (n <= n_max) {
      data <- tmle_spec_adapt$new_Gstar(gen_data, gen_data_adapt, by, old_data, node_list, initial_likelihood)
      tmle_task <- tmle_spec_adapt$make_tmle_task(data, node_list)

      # Define a new likelihood:
      initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)
      updater <- tmle_spec_adapt$make_updater()
      targeted_likelihood <- tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)

      tmle_params <- c(tmle_params, tmle_spec_adapt$make_params(tmle_task, targeted_likelihood))

      n <- n + by
    }

    fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
  }
}
