MC = 500
est_d0 <- data.frame(matrix(NA, nrow = MC, ncol = 40))
cov <- data.frame(matrix(NA, nrow = MC, ncol = 8))

tmle3_sadapt_sim1b <- function(surrogate = TRUE,
                               S = S,
                               W = c("W1","W2","W3"),
                               V = NULL,
                               A = "A",
                               Y = "Y",
                               learners = learner_list,
                               data,
                               param = "opt",
                               training_size = 100,
                               test_size = 20,
                               mini_batch = 25,
                               Gexploit = 0.1,
                               Gexplore = 0.01,
                               n_max = 1600,
                               by = 200,
                               n = nrow(data),
                               gen_data = gen_data,
                               gen_data_adapt = gen_data_adapt,
                               gen_data_adapt_truth = gen_data_adapt_truth,
                               rho = rho) {
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
    initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task,
                                                            learner_list)
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

    initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task,
                                                                  learner_list)
    updater <- tmle_spec_adapt$make_updater()
    targeted_likelihood <-
      tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)
    tmle_params <- tmle_spec_adapt$make_params(tmle_task, targeted_likelihood)

    fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)

    #Data-adaptive param:
    #data_new<-gen_data(n = 50000)
    #tmle_task_new <- self$make_tmle_task(data_new, node_list)
    #blip_task<-self$get_blip_cf(tmle_task_new)
    #dn<-self$get_Gstar(blip_task, initial_likelihood)
    #data_dn<-gen_data_adapt_truth(n = 50000, Gstar=dn, W=data_new[,1:3])
    #psi_dn

    #Check coverage:
    est_d0[i, 1:10] <- fit$summary
    cov[i, i] <- ifelse((fit$summary$lower_transformed <= psi &&
                         fit$summary$upper_transformed >= psi), 1, 0)
    #cov[i,i+1] <- ifelse((fit$summary$lower_transformed <= psi &&
                          #fit$summary$upper_transformed >= psi), 1, 0)

    while (n <= n_max) {
      data <- tmle_spec_adapt$new_Gstar(gen_data, gen_data_adapt, by, old_data,
                                        node_list, initial_likelihood)
      tmle_task <- tmle_spec_adapt$make_tmle_task(data, node_list)

      # Define a new likelihood:
      initial_likelihood <-
        tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)
      updater <- tmle_spec_adapt$make_updater()
      targeted_likelihood <-
        tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)

      tmle_params <- c(tmle_params,
                       tmle_spec_adapt$make_params(tmle_task,
                                                   targeted_likelihood))
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
    initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task,
                                                                  learner_list)
    updater <- tmle_spec_adapt$make_updater()
    targeted_likelihood <-
      tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)
    tmle_params <- tmle_spec_adapt$make_params(tmle_task, targeted_likelihood)

    while (n <= n_max) {
      data <- tmle_spec_adapt$new_Gstar(gen_data, gen_data_adapt, by, old_data,
                                        node_list,initial_likelihood)
      tmle_task <- tmle_spec_adapt$make_tmle_task(data, node_list)

      # Define a new likelihood:
      initial_likelihood <-
        tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)
      updater <- tmle_spec_adapt$make_updater()
      targeted_likelihood <-
        tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)

      tmle_params <- c(tmle_params,
                       tmle_spec_adapt$make_params(tmle_task,
                                                   targeted_likelihood))
      n <- n + by
    }
    fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
  }
}


