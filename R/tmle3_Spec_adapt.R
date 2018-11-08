#' Adaptively learns the Mean under the Optimal Individualized Treatment Rule or the
#' Average Treatmen Effect in a sequential trial, possibly by finding an optimal
#' surrogate outcome.
#'
#' @importFrom R6 R6Class
#'
#' @export
#

tmle3_Spec_adapt <- R6Class(
  classname = "tmle3_Spec_adapt",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(S = NULL, V = NULL, learners, param = "opt",
                              training_size, test_size, mini_batch,
                              Gexploit = 0.1, Gexplore = 0.05, ...) {
      options <- list(
        S = S, V = V, param = param, learners = learners,
        training_size = training_size, test_size = test_size, mini_batch = mini_batch,
        Gexploit = Gexploit, Gexplore = Gexplore
      )
      do.call(super$initialize, options)
    },

    ##################################
    # Function useful for simulations!
    new_Gstar = function(gen_data, gen_data_adapt, by, node_list, initial_likelihood) {
      data_new <- gen_data(n = by)
      tmle_task_new <- self$make_tmle_task(data_new, node_list, initial=TRUE)
      blip_task <- self$get_blip_cf(tmle_task_new)
      dn <- self$get_Gstar(blip_task, initial_likelihood)
      W <- data_new[, 1:3]

      data_targeted <- gen_data_adapt(n = by, Gstar = dn, W = W)

      return(data_targeted)
    },
    
    new_data = function(inter, old_data,tmle_spec){
      
      #Create the optimal surrogate:
      tmle_spec
      
    }

    bound = function(g) {
      g[g < 0.01] <- 0.01
      g[g > 0.99] <- 0.99
      return(g)
    },

    make_tmle_task = function(data, node_list, initial=FALSE, ...) {
      setDT(data)

      Y_node <- node_list$Y
      Y_vals <- unlist(data[, Y_node, with = FALSE])
      Y_variable_type <- variable_type(x = Y_vals)

      training_size <- self$get_training
      test_size <- self$get_test
      mini_batch <- self$get_batch

      if (Y_variable_type$type == "continuous") {
        min_Y <- min(Y_vals)
        max_Y <- max(Y_vals)
        range <- max_Y - min_Y
        lower <- min_Y
        upper <- max_Y
        Y_variable_type <- variable_type(
          type = "continuous",
          bounds = c(lower, upper)
        )
      }

      npsem <- list(
        define_node("W", node_list$W),
        define_node("A", node_list$A, c("W")),
        define_node("Y", node_list$Y, c("A", "W"), Y_variable_type)
      )

      if (initial) {
        tmle_task <- tmle3_Task$new(data, npsem = npsem)
      }
      else {
        #TO DO: Change this to folds_rolling_origin...
        folds <- origami::make_folds(data,
          fold_fun = folds_rolling_window,
          window_size = training_size,
          validation_size = test_size, gap = 0,
          batch = mini_batch
        )
        # TO DO: Add weights g^ref/g
        # tmle_task <- tmle3_Task$new(data, npsem = npsem, folds=folds)
        tmle_task <- tmle3_Task$new(data, npsem = npsem)
      }
      return(tmle_task)
    },
    make_updater = function() {
      updater <- tmle3_Update$new()
    },

    blik = function(A, G) {
      return(A * G + (1 - A) * (1 - G))
    },

    # Balanced design
    g_ref = function(n) {
      return(rep(0.5, n))
    },

    smoothIndicator = function(blip, Gexploit, Gexplore) {
      aa <- -(1 / 2 - Gexploit) / (2 * Gexplore^3)
      bb <- (1 / 2 - Gexploit) / (2 * Gexplore / 3)
      cc <- 1 / 2

      out <- rep(Gexploit, length(blip))
      pos <- (blip > Gexplore)
      btwn <- (-Gexplore <= blip & blip <= Gexplore)
      out[pos] <- 1 - Gexploit
      out[btwn] <- aa * blip[btwn]^3 + bb * blip[btwn] + cc

      return(out)
    },

    Q = function(task, likelihood) {
      pred <- likelihood$get_likelihood(tmle_task = task, node = "Y")
      return(pred)
    },

    get_Gstar = function(task, likelihood) {
      Gexploit <- self$get_Gexploit
      Gexplore <- self$get_Gexplore

      # TO DO: Limited to only binary at this point
      blip <- self$get_blip(task, likelihood)

      Gstar <- self$smoothIndicator(blip, Gexploit, Gexplore)

      return(Gstar)
    },

    get_weight = function(G_ref, GstarW, A) {

      # P(A=a|W)
      GA <- self$blik(A, G_ref)
      GAstarW <- self$blik(A, GstarW)

      return(GAstarW / GA)
    },

    get_blip = function(task, likelihood) {
      return(self$Q(task[[2]], likelihood) - self$Q(task[[1]], likelihood))
    },

    get_rule = function(task, likelihood) {
      return(as.numeric(self$Q(task[[2]], likelihood) - self$Q(task[[1]], likelihood) > 0))
    },

    get_blip_cf = function(tmle_task) {
      A_vals <- tmle_task$npsem$A$variable_type$levels

      # Generate counterfactual tasks for each value of A:
      cf_tasks <- lapply(A_vals, function(A_val) {
        newdata <- data.table(A = A_val)
        cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = newdata)
        return(cf_task)
      })
      return(cf_tasks)
    },

    make_params = function(tmle_task, likelihood) {
      n <- nrow(tmle_task$data)
      data <- tmle_task$get_data()
      A <- data$A
      Y <- data$Y

      A_vals <- tmle_task$npsem$A$variable_type$levels

      # Generate counterfactual tasks for each value of A:
      cf_tasks <- lapply(A_vals, function(A_val) {
        newdata <- data.table(A = A_val)
        cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = newdata)
        return(cf_task)
      })

      # TO DO: Implement a better way to learn the rule
      # TO DO: For now, the only reference design is a balanced design
      GstarW <- self$get_Gstar(task = cf_tasks, likelihood)
      G_ref <- self$g_ref(n)

      # Get weights
      weight <- self$get_weight(G_ref, GstarW, A)

      # Learn the rule:
      rA <- self$get_rule(task = cf_tasks, likelihood)

      # How to incorporate weights?
      lf_rule <- define_lf(LF_rule, "A", rule_fun = rA)
      intervens <- Param_TSM$new(observed_likelihood = likelihood, intervention_list = lf_rule)
      # intervens <- Param_TSM_weight$new(observed_likelihood=likelihood,
      #                                 weight=weight,
      #                                intervention_list=lf_rule)

      return(intervens)
    }
  ),

  active = list(
    get_training = function() {
      t <- private$.options$training_size
      return(t)
    },

    get_test = function() {
      t <- private$.options$test_size
      return(t)
    },

    get_batch = function() {
      t <- private$.options$mini_batch
      return(t)
    },

    get_Gexploit = function() {
      t <- private$.options$Gexploit
      return(t)
    },

    get_Gexplore = function() {
      t <- private$.options$Gexplore
      return(t)
    },

    get_S = function() {
      S <- private$.options$S
      return(S)
    },

    get_V = function() {
      V <- private$.options$V
      return(V)
    },

    get_S_learner = function() {
      S_lrn <- private$.options$learners$S
      return(S_lrn)
    },

    get_B_learner = function() {
      B_lrn <- private$.options$learners$B
      return(B_lrn)
    },

    get_learners = function() {
      lrn <- private$.options$learners
      return(lrn)
    },

    get_param = function() {
      param <- private$.options$param
    }
  ),
  private = list()
)

#' Adaptively learns the Mean under the Optimal Individualized Treatment Rule or the
#' Average Treatmen Effect in a sequential trial, possibly by finding an optimal
#' surrogate outcome.
#'
#' O=(W,A,S,Y)
#' W=Covariates
#' A=Treatment (binary)
#' S=Potential Surrogates
#' Y=Outcome (binary or bounded continuous)
#'
#' @param S Covariates to consider for the Optimal Surrogate estimation. Leave empty if no surrage is used.
#' @param learners List of learners used for Q,g,S and B.
#' @param param Target parameter. Current implementation supports Mean under the Optimal Individualized
#' Treatment (opt) and Average Treatment Effect (are)
#' @param V Covariates the rule depends on. Not used at the moment.
#' @param training_size Size of the initial training set. Necessary part of online Super Learner.
#' @param test_size Size of the test set. Necessary part of online Super Learner.
#' @param mini_batch Size of the increase in the initial training size, added per each iteration of the
#' online Super Learner.
#' @param Gexploit Sequence t_n. Default is 0.1.
#' @param Gexplore Sequence e_n. Default is 0.05.
#'
#' @export
#'

tmle3_adapt <- function(S = NULL, V = NULL, learners, param = "opt",
                        training_size, test_size, mini_batch,
                        Gexploit = 0.1, Gexplore = 0.05) {
  tmle3_Spec_adapt$new(
    surrogate = surrogate, S = S, V = V, learners = learners, param = param,
    training_size = training_size, test_size = test_size, mini_batch = mini_batch,
    Gexploit = Gexploit, Gexplore = Gexplore
  )
}
