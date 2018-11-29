#Function to run simulations for a new trial (sampling from the same population)
#Optimal surrogate comparison

tmle3_new_trial_osc <- function(surrogate = TRUE,
                                #Column names for the surrogates
                                S = S,
                                #Column names for the baseline covariates
                                W = c("W1","W2","W3"),
                                #Covariates the rule depends on. 
                                #Not implemented at the momemnt
                                V = NULL,
                                #Column name for the treatment
                                A = "A",
                                #Column name for the outcome of interest
                                Y = "Y",
                                #List of learners for the Online SL
                                learners = learner_list,
                                #Initial data for the 1st trial (if available)
                                data = NULL,
                                #Target parameter (opt vs. ate)
                                param = "opt",
                                #Initial training size for the Online SL
                                training_size = 100,
                                #Validation size for the Online SL
                                test_size = 20,
                                #Increase size for the Online SL
                                mini_batch = 25,
                                #Sequence t for OIT
                                Gexploit = 0.1,
                                #Sequence e for OIT
                                Gexplore = 0.01,
                                #Sample size of the first trial
                                n1 = 500,
                                #Sample size of the second trial
                                n2 = 500,
                                #Function to generate data of the format
                                #O=(W,A,S,Y,Y1,Y0)
                                gen_data = gen_data,
                                #Used only for the OIT parameter
                                #Rule learned w.r.t. Y or Optimal Surrogate?
                                rule_outcome=rule_outcome,
                                #Estimator of the Optimal Surrogate (TMLE vs. SL)
                                opt_surrogate=opt_surrogate,
                                #If we know the truth, put here:
                                psi=NULL
                                ) {
  if (surrogate) {
    
    mse <- NULL
    
    #Option to play with different initial sample size.
    if(is.null(data)){
      data<-gen_data(n = n1)
    }

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
    ## Target Towards the Target Parameter
    ########################################
    
    if(param == "ate"){
      
      tmle_spec_ate <- tmle_ate(baseline = 0, contrast = 1)
      
      ####################################
      ### Actual Y:
      ####################################
      
      tmle_task_Y <- tmle_spec_ate$make_tmle_task(data, node_list)
      
      #Estimates:
      #E(Y|A,W), P(A=1|W)
      initial_likelihood_Y <- tmle_spec_ate$make_initial_likelihood(tmle_task_Y,
                                                                    learner_list)
      updater_Y <- tmle_spec_ate$make_updater()
      targeted_likelihood_Y <-
        tmle_spec_ate$make_targeted_likelihood(initial_likelihood_Y, updater_Y)
      tmle_params_Y <- tmle_spec_ate$make_params(tmle_task_Y, targeted_likelihood_Y)
      
      fit_Y <- fit_tmle3(tmle_task_Y, targeted_likelihood_Y, tmle_params_Y, updater_Y)

      ####################################
      ## Targeted optimal surrogate:
      ####################################
      
      tmle_task_tmle <- tmle_spec_ate$make_tmle_task(data_tmle, node_list)
      
      #Estimates:
      #E(Y_s^*|A,W), P(A=1|W)
      initial_likelihood_tmle <- tmle_spec_ate$make_initial_likelihood(tmle_task_tmle,
                                                                         learner_list)
      updater_tmle <- tmle_spec_ate$make_updater()
      targeted_likelihood_tmle <-
        tmle_spec_ate$make_targeted_likelihood(initial_likelihood_tmle, updater_tmle)
      tmle_params_tmle <- tmle_spec_ate$make_params(tmle_task_tmle, targeted_likelihood_tmle)
      
      fit_tmle <- fit_tmle3(tmle_task_tmle, targeted_likelihood_tmle, tmle_params_tmle, updater_tmle)
      
      ####################################
      ## Super-Learner optimal surrogate:
      ####################################
      
      tmle_task_sl <- tmle_spec_ate$make_tmle_task(data_sl, node_list)
      
      #Estimates:
      #E(Y_s|A,W), P(A=1|W)
      initial_likelihood_sl <- tmle_spec_ate$make_initial_likelihood(tmle_task_sl,
                                                                       learner_list)
      updater_sl <- tmle_spec_ate$make_updater()
      targeted_likelihood_sl <-
        tmle_spec_ate$make_targeted_likelihood(initial_likelihood_sl, updater_sl)
      tmle_params_sl <- tmle_spec_ate$make_params(tmle_task_sl, targeted_likelihood_sl)
      
      fit_sl <- fit_tmle3(tmle_task_sl, targeted_likelihood_sl, tmle_params_sl, updater_sl)
      
      #MSE in the original trial:
      mse <- c(Y=(psi-fit_Y$summary$tmle_est[[3]])^2,
              TMLE=(psi-fit_tmle$summary$tmle_est[[3]])^2,
               SL=(psi-fit_sl$summary$tmle_est[[3]])^2)
      
      summary_Y <- fit_Y$summary[3,]
      summary_tmle <- fit_tmle$summary[3,]
      summary_sl <- fit_sl$summary[3,]

    }else if(param == "opt"){
      
    }
    
    ########################################
    ## New Trial
    ########################################
    
    #Sample data from a new trial:
    data_new <- gen_data(n=n2)
    
    #Set up spec for ATE:
    tmle_spec_ate <- tmle_ate(baseline = 0, contrast = 1)
    
    ####################################
    ### Actual Y:
    ####################################
    
    tmle_task_Y <- tmle_spec_ate$make_tmle_task(data_new, node_list)
    
    #Estimates:
    #E(Y|A,W), P(A=1|W)
    initial_likelihood_Y <- tmle_spec_ate$make_initial_likelihood(tmle_task_Y,
                                                                  learner_list)
    updater_Y <- tmle_spec_ate$make_updater()
    targeted_likelihood_Y <-
      tmle_spec_ate$make_targeted_likelihood(initial_likelihood_Y, updater_Y)
    tmle_params_Y <- tmle_spec_ate$make_params(tmle_task_Y, targeted_likelihood_Y)
    
    fit_Y_new <- fit_tmle3(tmle_task_Y, targeted_likelihood_Y, tmle_params_Y, updater_Y)
    fit_Y <- c(fit_Y,fit_Y_new)
    
    ####################################
    ## Targeted optimal surrogate:
    ####################################
    
    data_tmle <- tmle_spec_tmle$get_targeted_surrogate(inter=data_new, param=param)
    
    tmle_task_tmle <- tmle_spec_ate$make_tmle_task(data_tmle, node_list)
    
    #Estimates:
    #E(Y_s^*|A,W), P(A=1|W)
    initial_likelihood_tmle <- tmle_spec_ate$make_initial_likelihood(tmle_task_tmle,
                                                                     learner_list)
    updater_tmle <- tmle_spec_ate$make_updater()
    targeted_likelihood_tmle <-
      tmle_spec_ate$make_targeted_likelihood(initial_likelihood_tmle, updater_tmle)
    tmle_params_tmle <- tmle_spec_ate$make_params(tmle_task_tmle, targeted_likelihood_tmle)
    
    fit_tmle_new <- fit_tmle3(tmle_task_tmle, targeted_likelihood_tmle, tmle_params_tmle, updater_tmle)
    fit_tmle <- c(fit_tmle,fit_tmle_new)
    
    ####################################
    ## Super-Learner optimal surrogate:
    ####################################
    
    data_sl <- tmle_spec_tmle$get_SL_surrogate(inter=data_new)
    
    tmle_task_sl <- tmle_spec_ate$make_tmle_task(data_sl, node_list)
    
    #Estimates:
    #E(Y_s|A,W), P(A=1|W)
    initial_likelihood_sl <- tmle_spec_ate$make_initial_likelihood(tmle_task_sl,
                                                                   learner_list)
    updater_sl <- tmle_spec_ate$make_updater()
    targeted_likelihood_sl <-
      tmle_spec_ate$make_targeted_likelihood(initial_likelihood_sl, updater_sl)
    tmle_params_sl <- tmle_spec_ate$make_params(tmle_task_sl, targeted_likelihood_sl)
    
    fit_sl_new <- fit_tmle3(tmle_task_sl, targeted_likelihood_sl, tmle_params_sl, updater_sl)
    fit_sl <- c(fit_sl,fit_sl_new)
    
    #MSE in the original trial:
    mse_new <- c(Y=(psi-fit_Y_new$summary$tmle_est[[3]])^2,
             TMLE=(psi-fit_tmle_new$summary$tmle_est[[3]])^2,
             SL=(psi-fit_sl_new$summary$tmle_est[[3]])^2)
    
    mse<-c(mse,mse_new)
    
    summary_Y <- rbind(summary_Y, fit_Y_new$summary[3,])
    summary_tmle <- rbind(summary_tmle, fit_tmle_new$summary[3,])
    summary_sl <- rbind(summary_sl, fit_sl_new$summary[3,])
    
    return(list(summary_Y=summary_Y,
                summary_tmle=summary_tmle,
                summary_sl=summary_sl,
                mse=mse,
                fit_Y=fit_Y,
                fit_tmle=fit_tmle,
                fit_sl=fit_sl))
  }
}



  
  
  
  
  
  