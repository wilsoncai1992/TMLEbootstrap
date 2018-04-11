#' @export
avgDensity_LambdaGrid <- R6Class("avgDensity_LambdaGrid",
  public = list(
    data = NULL,
    bin_width = NULL,
    epsilon_step = NULL,
    REPEAT_BOOTSTRAP = NULL,
    inflate_lambda = NULL,

    dict_boot = dict::dict(), # dictionary of comprehensiveBootstrap
    initialize = function(data,
                          bin_width,
                          epsilon_step,
                          REPEAT_BOOTSTRAP = 2e2,
                          inflate_lambda = 1) {
      library(dict)
      self$data <- data
      self$bin_width <- bin_width
      self$epsilon_step <- epsilon_step
      self$REPEAT_BOOTSTRAP <- REPEAT_BOOTSTRAP
      self$inflate_lambda <- inflate_lambda
    },
    add_lambda = function(lambda_grid = NULL) {
      for (lambda in lambda_grid) {
        boot_here <- comprehensiveBootstrap$new(parameter = avgDensityBootstrap,
                                                x = self$data$x,
                                                bin_width = self$bin_width,
                                                epsilon_step = self$epsilon_step)
        boot_here$bootstrap(REPEAT_BOOTSTRAP = self$REPEAT_BOOTSTRAP,
                            inflate_lambda = self$inflate_lambda)
        boot_here$compute_width()
        self$dict_boot[[lambda]] <- boot_here
        message(paste(lambda, 'is added'))
      }
    }
  )
)