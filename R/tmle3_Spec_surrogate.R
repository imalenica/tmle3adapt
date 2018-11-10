#' Learns the Optimal Surrogate based on the observed final outcome Y, for the
#' Mean Under the Optimal Individualized Rule and Average Treatment Effect
#' target parameters.
#'
#' @importFrom R6 R6Class
#'
#' @export
#

tmle3_Spec_surrogate <- R6Class(
  classname = "tmle3_Spec_surrogate",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    # TO DO: Should we learn the initial estimates using Online Learning?
    # Theoretically, no need
    initialize = function(S, V = NULL, learners, param = "opt",
                              training_size = NULL, test_size = NULL, mini_batch = NULL, ...) {
      options <- list(
        S = S, V = V, param = param, learners = learners,
        training_size = training_size, test_size = test_size, mini_batch = mini_batch
      )
      do.call(super$initialize, options)
    },

    make_tmle_task = function(data, node_list, ...) {
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

      if (!is.null(node_list$id)) {
        tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, ...)
      }
      else {
        # folds <- origami::make_folds(data,
        #  fold_fun = origami::folds_rolling_window,
        #  window_size = training_size,
        #  validation_size = test_size, gap = 0,
        #  batch = mini_batch
        # )
        # TO DO: Add weights g^ref/g
        tmle_task <- tmle3_Task$new(data, npsem = npsem)
      }
      return(tmle_task)
    },

    get_surrogate = function(inter) {
      # Get the fit we learned in the 1st part of the trial:
      osl <- self$get_sur_sl
      sur_sl <- osl$get_sur_sl

      # TO DO: Y must be last for this to work!
      covariates <- c(names(inter)[-length(inter)])

      sur_tmle_task <- make_sl3_Task(inter, covariates = covariates, outcome = "Y")
      S_pred <- sur_sl$predict(sur_tmle_task)
    },

    make_params = function(tmle_task, likelihood) {
      S <- self$get_S
      V <- self$get_V
      learners <- self$get_learners
      param <- self$get_param

      opt <- Optimal_Surrogate$new(
        S = S, V = V, learners = learners, param = param,
        tmle_task = tmle_task, likelihood = likelihood
      )

      ### Learn the SL optimal surrogate, and save the fit:
      S_pred <- opt$surrogate_SL()
      private$sur_sl <- opt$get_sur_sl
      private$opt <- opt

      ### Target towards the parameter of interest:
      Starg_pred <- opt$surrogate_TSL(S_pred = S_pred)

      data <- tmle_task$get_data()
      data$Y <- Starg_pred

      return(data)
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
    },
    get_sur_sl = function() {
      return(private$sur_sl)
    },
    get_opt = function() {
      return(private$opt)
    }
  ),
  private = list(
    sur_sl = list(),
    opt = list()
  )
)

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
#' @param S Covariates to consider for the Optimal Surrogate estimation.
#' @param learners List of learners used for Q,g,S and B.
#' @param param Target parameter. Current implementation supports Mean under the Optimal Individualized
#' Treatment (opt) and Average Treatment Effect (are)
#' @param V Covariates the rule depends on.
#' @param training_size Size of the initial training set. Necessary part of online Super Learner.
#' @param test_size Size of the test set. Necessary part of online Super Learner.
#' @param mini_batch Size of the increase in the initial training size, added per each iteration of the
#' online Super Learner.
#'
#'
#' @export
#'

tmle3_surrogate <- function(S, V = NULL, learners, param = "opt", training_size = NULL,
                            test_size = NULL, mini_batch = NULL) {
  tmle3_Spec_surrogate$new(
    S = S, V = V, learners = learners, param = param,
    training_size = training_size, test_size = test_size, mini_batch = mini_batch
  )
}
