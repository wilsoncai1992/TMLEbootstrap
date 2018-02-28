library(R6)

#' @export
blipVarianceTMLE <- R6Class("blipVarianceTMLE",
  public = list(
    data = NULL,
    Psi = NULL,
    EIC = NULL,
    epsilon_step = 1e-3,
    tol = NULL,
    CI = NULL,
    verbose = FALSE,
    max_iter = 1e2,
    initialize = function(data, epsilon_step = NULL, verbose = NULL) {
      self$data <- data
      self$tol <- 1/length(x)
      if (!is.null(epsilon_step)) self$epsilon_step <- epsilon_step
      if (!is.null(verbose)) self$verbose <- verbose
    }
))


