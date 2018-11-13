
tmle3_sadapt_sim1b <- function(surrogate = FALSE,
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
                               rho = rho, rule_outcome=rule_outcome,
                               opt_surrogate=opt_surrogate) {
  if (surrogate) {
    
    ## Learn the optimal surrogate:
    tmle_spec <- tmle3_surrogate(
      S = S,
      V = V,
      learners = learner_list,
      param = param, 
      opt_surrogate="SL",
      rule_outcome=rule_outcome
    )

    # Define nodes:
    node_list <- list(W = W, A = A, Y = Y)

    # Define data:
    tmle_task <- tmle_spec$make_tmle_task(data, node_list)

    # Define likelihood: (only need P(A=1|W))
    #Estimates:
    #E(Y|A,W), P(A=1|W), P(W)
    initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task,
                                                            learner_list)
    data <- tmle_spec$make_params(tmle_task, initial_likelihood)

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

    #COMPUTATIONAL: Can borrow from prior spec?
    #Estimates:
    #E(Y_s|A,W), P(A=1|W)
    initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task,
                                                                  learner_list)
    updater <- tmle_spec_adapt$make_updater()
    targeted_likelihood <-
      tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)
    tmle_params <- tmle_spec_adapt$make_params(tmle_task, targeted_likelihood)

    #fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)

    #Data-adaptive param:
    #data_new<-gen_data(n = 50000)
    #tmle_task_new <- self$make_tmle_task(data_new, node_list)
    #blip_task<-self$get_blip_cf(tmle_task_new)
    #dn<-self$get_Gstar(blip_task, initial_likelihood)
    #data_dn<-gen_data_adapt_truth(n = 50000, Gstar=dn, W=data_new[,1:3])
    #psi_dn

    while (n < n_max) {
      #Returns targeted data
      data_targ <- tmle_spec_adapt$new_Gstar(gen_data, gen_data_adapt, W=NULL, by,
                                        node_list, initial_likelihood)
      #Returns old and new (targeted) data, with outcome being the surrogate as learned previously
      data<-tmle_spec_adapt$new_data(inter=data_targ,old_data=data, 
                                     tmle_spec=tmle_spec, node_list=node_list)

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
    
   } else if(!surrogate) {
    
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
    targeted_likelihood <- tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)
    tmle_params <- tmle_spec_adapt$make_params(tmle_task, targeted_likelihood)

    while (n <= n_max) {
      data_new <- tmle_spec_adapt$new_Gstar(gen_data, gen_data_adapt, by, old_data,
                                        node_list,initial_likelihood)
      data<-rbind.data.frame(data,data_new)
      tmle_task <- tmle_spec_adapt$make_tmle_task(data, node_list)

      # Define a new likelihood:
      initial_likelihood <- tmle_spec_adapt$make_initial_likelihood(tmle_task, learner_list)
      updater <- tmle_spec_adapt$make_updater()
      targeted_likelihood <- tmle_spec_adapt$make_targeted_likelihood(initial_likelihood, updater)

      tmle_params <- c(tmle_params, tmle_spec_adapt$make_params(tmle_task,
                                                                 targeted_likelihood))
      n <- n + by
    }
    
    fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
  }
  
  return(fit)
  
}


