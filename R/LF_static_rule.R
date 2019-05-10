#' Dynamic Likelihood Factor with static rule
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_static, name, type, value, ...)}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name 
#'     in the nodes specified by tmle3_Task.}
#'     
#'     \item{\code{type}}{character, either 'density', for conditional density or, 'mean' for conditional mean
#'     }
#'     \item{\code{value}}{the static value
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{value}}{the static value.}
#'     }

#'
#' @export

LF_static_rule <- R6Class(
  classname = "LF_static_rule",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, type = "density", value, ...) {
      super$initialize(name, ..., type = type)
      private$.value <- value
      private$.variable_type <- variable_type("constant", value)
    },
    get_mean = function(tmle_task, fold_number=-1) {
      return(self$value)
    },
    get_density = function(tmle_task, fold_number=-1) {
      observed <- tmle_task$get_tmle_node(self$name)
      likelihood <- as.numeric(self$value == observed)
      
      return(likelihood)
    },
    cf_values = function(tmle_task) {
      cf_values <- self$value
      return(cf_values)
    }
  ),
  active = list(
    value = function() {
      return(private$.value)
    }
  ),
  private = list(
    .name = NULL,
    .value = NULL,
    .is_degenerate = TRUE
  )
)