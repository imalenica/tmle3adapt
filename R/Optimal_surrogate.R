#' Learns the Optimal Surrogate based on the observed final outcome Y, for the
#' Mean Under the Optimal Individualized Rule and Average Treatment Effect
#' target parameters.
#'
#' @importFrom R6 R6Class
#' @importFrom data.table data.table
#'
#' @export
#
Optimal_Surrogate <- R6Class(
  classname = "Optimal_Surrogate",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  lock_objects = FALSE,
  public = list(
    initialize = function(S, V = NULL, learners, param = "opt", tmle_task, likelihood) {
      private$.S <- S
      private$.V <- V
      private$.learners <- learners
      private$.param <- param
      private$.tmle_task <- tmle_task
      private$.likelihood <- likelihood
    },

    bound = function(g) {
      g[g < 0.01] <- 0.01
      g[g > 0.99] <- 0.99
      return(g)
    },

    Q = function(task) {
      Q_sl <- self$get_Q_pred
      return(Q_sl$predict(task))
    },

    get_rule = function(task) {
      return(as.numeric(self$Q(task[[2]]) - self$Q(task[[1]]) > 0))
    },

    # SuperLearner of the surrogate:
    surrogate_SL = function() {
      tmle_task <- self$tmle_task

      data <- tmle_task$get_data()
      S_learner <- self$get_S_learner
      S <- self$get_S
      # folds<-tmle_task$folds

      covariates <- c(S, names(tmle_task$data)[-length(tmle_task$data)])

      sur_tmle_task <- make_sl3_Task(data, covariates = covariates, outcome = "Y")
      sur_sl <- S_learner$train(sur_tmle_task)

      S_pred <- sur_sl$predict()
      private$sur_sl <- sur_sl

      return(S_pred)
    },

    # Targeted SuperLearner of the surrogate:
    surrogate_TSL = function(S_pred, tmle_task) {
      param <- self$get_param
      tmle_task <- self$tmle_task
      initial_likelihood <- self$likelihood

      if (param == "opt") {

        # Simple rule learning:
        data <- tmle_task$get_data()
        data$Y <- S_pred

        A <- data$A
        Y <- data$Y
        S <- self$get_S

        folds <- tmle_task$folds

        covariates <- c(S, names(tmle_task$data)[-length(tmle_task$data)])
        Q_tmle_task <- make_sl3_Task(data, covariates = covariates, outcome = "Y", folds = folds)

        B_learner <- self$get_B_learner
        Q_sl <- B_learner$train(Q_tmle_task)
        private$Q_sl <- Q_sl

        A_vals <- tmle_task$npsem$A$variable_type$levels

        # Generate counterfactual tasks for each value of A:
        cf_tasks <- lapply(A_vals, function(A_val) {
          newdata <- data
          newdata$A <- A_val
          cf_task <- make_sl3_Task(newdata,
            covariates = covariates,
            outcome = "Y", folds = folds
          )
          return(cf_task)
        })

        # Learn the rule:
        dn <- self$get_rule(cf_tasks)

        # Get the estimated g:
        g.ests <- initial_likelihood$get_likelihood(tmle_task = tmle_task, node = "A")
        g.ests[A == 0] <- 1 - g.ests[A == 0]
        g.ests <- self$bound(g.ests)

        # Clever covariate and fluctuation:
        HA <- as.numeric(A == dn) / g.ests
        # Q.ests vs. Y...
        eps <- coef(glm(Y ~ -1 + HA, offset = qlogis(Y), family = "quasibinomial"))

        # Update:
        Q.star <- plogis(qlogis(Y) + HA * eps)

        # Use tmle3mopptx to get the rule:
        # learner_list <- self$get_learners
        # B_learner <- self$get_B_learner
        # V <- self$get_V

        # data <- tmle_task$get_data()
        # data$Y <- S_pred
        # A<-data$A
        # Y<-data$Y

        # tmle_spec_opt <- tmle3_mopttx_blip(V = V, type = "blip1",b_learner = B_learner, maximize = TRUE, complex = FALSE)

        # node_list <- list(W = c(tmle_task$npsem$W$variables),
        #                  A = tmle_task$npsem$A$variables,
        #                  Y = tmle_task$npsem$Y$variables)

        # TO DO: Right now we are "refitting" Q, completly unnecessary
        # tmle_task_opt <- tmle_spec_opt$make_tmle_task(data, node_list)
        # initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task_opt, learner_list)
        # opt_rule <- Optimal_Rule$new(tmle_task_opt, initial_likelihood, "split-specific",
        #                            V = V, blip_type = "blip1",
        #                             blip_library = B_learner, maximize = TRUE)
        # opt_rule$fit_blip()

        # Get the rule:
        # dn <- opt_rule$rule(tmle_task_opt)

        # Get the estimated g:
        # g.ests <- initial_likelihood$get_likelihood(tmle_task = tmle_task_opt, node = "A")
        # g.ests[A==0] <- 1-g.ests[A==0]
        # g.ests <- self$bound(g.ests)

        # Clever covariate and fluctuation:
        # HA <- as.numeric(A == dn)/g.ests
        # Q.ests vs. Y...
        # eps<-coef(glm(Y ~ -1 + HA, offset=qlogis(Y), family="quasibinomial"))

        # Update:
        # Q.star<-plogis(qlogis(Y) + HA * eps)

        return(Q.star)

        # fit <- tmle3(tmle_spec_opt, data, node_list, learner_list=private$.options$learners)
        # Q<-initial_likelihood$get_likelihood(tmle_task_opt, "Y", -1)
        # e<-fit$updater$epsilons[[1]]$Y
        # H<-fit$tmle_params[[1]]$clever_covariates()
      } else if (param == "ate") {

        # TO DO
      } else {
        stop("The specified parameter is not implemented. The package currently only
             supports Mean under the Optimal Individualized Treatment (opt) and 
             Average Treatment Effect (ate)")
      }
    }
  ),
  active = list(
    tmle_task = function() {
      return(private$.tmle_task)
    },
    likelihood = function() {
      return(private$.likelihood)
    },
    V = function() {
      return(private$.V)
    },
    get_S = function() {
      S <- private$.S
      return(S)
    },
    get_V = function() {
      V <- private$.V
      return(V)
    },
    get_S_learner = function() {
      S_lrn <- private$.learners$S
      return(S_lrn)
    },
    get_B_learner = function() {
      B_lrn <- private$.learners$B
      return(B_lrn)
    },
    get_learners = function() {
      lrn <- private$.learners
      return(lrn)
    },
    get_param = function() {
      param <- private$.param
    },
    get_Q_pred = function() {
      param <- private$Q_sl
    },
    get_S_pred = function() {
      param <- private$sur_sl
    }
  ),
  private = list(
    sur_sl = list(),
    Q_sl = list(),
    .S = NULL,
    .V = NULL,
    .param = NULL,
    .learners = NULL,
    .tmle_task = list(),
    .likelihood = list()
  )
)
