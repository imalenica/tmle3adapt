#Function to run simulations for the adaptive sequential design.
#Optimal surrogate comparison

tmle3_adapt_sim_osc <- function(surrogate = FALSE,
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
                            rule_outcome=rule_outcome,
                            opt_surrogate=opt_surrogate,
                            data_adaptive=FALSE) {
  if (surrogate) {
    
    psi_dn <- NULL
    mse <- NULL
    data_original_round1<-data
    
    ########################################
    ## Learn the targeted optimal surrogate: 
    ########################################
    
    tmle_spec_tmle <- tmle3_surrogate(
      S = S,
      V = V,
      learners = learner_list,
      param = param, 
      opt_surrogate="TMLE",
      rule_outcome=rule_outcome
    )
    
    # Define nodes:
    node_list <- list(W = W, A = A, Y = Y)
    
    # Define data:
    tmle_task <- tmle_spec_tmle$make_tmle_task(data, node_list)
    
    # Define likelihood: (only need P(A=1|W))
    #Estimates:
    #E(Y|A,W), P(A=1|W), P(W)
    initial_likelihood <- tmle_spec_tmle$make_initial_likelihood(tmle_task,learner_list)
    data_tmle <- tmle_spec_tmle$make_params(tmle_task, initial_likelihood)
    
    ########################################
    ## Learn the SL optimal surrogate: 
    ########################################
    
    tmle_spec_sl <- tmle3_surrogate(
      S = S,
      V = V,
      learners = learner_list,
      param = param, 
      opt_surrogate="SL",
      rule_outcome=rule_outcome
    )

    data_sl <- tmle_spec_sl$make_params(tmle_task, initial_likelihood)
    
    ########################################
    ## Adapt
    ########################################

    ######
    ##TMLE
    ######
    
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
    
    tmle_task_tmle <- tmle_spec_adapt$make_tmle_task(data = data_tmle, node_list)
    
    #Estimates:
    #E(Y_s|A,W), P(A=1|W)
    initial_likelihood_tmle <- tmle_spec_adapt$make_initial_likelihood(tmle_task_tmle,
                                                                      learner_list)
    updater_tmle <- tmle_spec_adapt$make_updater()
    targeted_likelihood_tmle <-
      tmle_spec_adapt$make_targeted_likelihood(initial_likelihood_tmle, updater_tmle)
    tmle_params_tmle <- tmle_spec_adapt$make_params(tmle_task_tmle, targeted_likelihood_tmle)
    
    fit_tmle <- fit_tmle3(tmle_task_tmle, targeted_likelihood_tmle, tmle_params_tmle, updater_tmle)
    
    #####
    ##SL
    #####
    
    tmle_task_sl <- tmle_spec_adapt$make_tmle_task(data = data_sl, node_list)
    
    #Estimates:
    #E(Y_s|A,W), P(A=1|W)
    initial_likelihood_sl <- tmle_spec_adapt$make_initial_likelihood(tmle_task_sl,
                                                                       learner_list)
    updater_sl <- tmle_spec_adapt$make_updater()
    targeted_likelihood_sl <-
      tmle_spec_adapt$make_targeted_likelihood(initial_likelihood_sl, updater_sl)
    tmle_params_sl <- tmle_spec_adapt$make_params(tmle_task_sl, targeted_likelihood_sl)
    
    fit_sl <- fit_tmle3(tmle_task_sl, targeted_likelihood_sl, tmle_params_sl, updater_sl)
    
    #Data-adaptive param:
    if(data_adaptive==TRUE){
      
      #NOTE: This is data-adaptive parameter w.r.t the TRUE Y.
      # We want to see how close are we to the true data-adaptive parameter, 
      # had we actually obsereved Y.
      
      data_adapt_full<-gen_data(n = 100000)
      tmle_task_new <- tmle_spec_adapt$make_tmle_task(data_adapt_full, node_list)
      blip_task<-tmle_spec_adapt$get_blip_cf(tmle_task_new)
      dn<-tmle_spec_adapt$get_Gstar(blip_task, initial_likelihood)
      data_dn<-gen_data_adapt_truth(n = 100000, Gstar=dn, W=data.frame(data_adapt_full[,1:3]))
      psi_dn<-mean(data_dn$Yd0)
    }
    
    #For MSE later:
    mse <- c(TMLE=(psi_dn-fit_tmle$summary$tmle_est)^2,SL=(psi_dn-fit_sl$summary$tmle_est)^2)
    
    while (n < n_max) {
      
      ########
      ## TMLE
      ########
      
      #Returns new randomization probabilities; obtained by knowing E(Y_t|W,A)
      #Fix W and true Y!
      data_targ_tmle <- tmle_spec_adapt$new_Gstar(gen_data, gen_data_adapt, W=NULL, by,
                                             node_list, initial_likelihood_tmle)
      #Returns old and new (targeted) data, with outcome being the surrogate as learned previously
      data_tmle<-tmle_spec_adapt$new_data(inter=data_targ_tmle,old_data=data_tmle, 
                                     tmle_spec=tmle_spec_tmle, node_list=node_list)
      
      #Create a data-set that contains the data with actual Ys.
      data_original<-rbind.data.frame(data_original_round1,data_targ_tmle)
        
      #COMPUTATIONAL: Can we extract non-targeted initial likelihood from fit object?
      tmle_task_data_adapt <- tmle_spec_adapt$make_tmle_task(data_original, node_list)
      initial_likelihood_new <- tmle_spec_adapt$make_initial_likelihood(tmle_task_data_adapt, learner_list)
      
      fit_tmle_new <- tmle3(tmle_spec_adapt, data_tmle, node_list, learner_list)
      fit_tmle <- c(fit_tmle,fit_tmle_new)
      
      #Data-adaptive param:
      if(data_adaptive==TRUE){
        
        data_adapt_full<-gen_data(n = 100000)
        tmle_task_new <- tmle_spec_adapt$make_tmle_task(data_adapt_full, node_list)
        blip_task<-tmle_spec_adapt$get_blip_cf(tmle_task_new)
        dn<-tmle_spec_adapt$get_Gstar(blip_task, initial_likelihood_new)
        data_dn<-gen_data_adapt_truth(n = 100000, Gstar=dn, W=data.frame(data_adapt_full[,1:3]))
        psi_dn_new<-mean(data_dn$Yd0)
        psi_dn<-c(psi_dn,psi_dn_new)
      }
      
      ########
      ## SL
      ########
      
      #Returns new randomization probabilities; obtained by knowing E(Y_s|W,A)
      #Fix W and true Y!
      data_targ_sl <- tmle_spec_adapt$new_Gstar(gen_data=NULL, gen_data_adapt=NULL, 
                                                W=tmle_spec_adapt$get_newW, by,
                                                node_list, initial_likelihood_sl)
      #Replace so that S matches:
      data_targ_sl[,5:7]<-data_targ_tmle[,5:7]
      #Returns old and new (targeted) data, with outcome being the surrogate as learned previously
      data_sl<-tmle_spec_adapt$new_data(inter=data_targ_sl,old_data=data_sl, 
                                          tmle_spec=tmle_spec_sl, node_list=node_list)
      fit_sl_new<-tmle3(tmle_spec_adapt, data_sl, node_list, learner_list)
      fit_sl <- c(fit_sl,fit_sl_new)
      
      #For MSE later:
      mse_new <- c(TMLE=(psi_dn_new-fit_tmle_new$summary$tmle_est)^2,
                   SL=(psi_dn_new-fit_sl_new$summary$tmle_est)^2)
      mse<-c(mse,mse_new)
      
      n <- n + by
    }
    
    summary_tmle <- do.call(rbind, lapply(fit_tmle, function(x) x$summary))
    summary_sl <- do.call(rbind, lapply(fit_sl, function(x) x$summary))
    
  } 
  
  return(list(summary_tmle=summary_tmle,
              summary_sl=summary_sl,
              psi_dn=psi_dn, 
              mse=mse,
              fit_tmle=fit_tmle,
              fit_sl=fit_sl))
  
}


