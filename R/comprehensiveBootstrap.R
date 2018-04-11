#' @export
comprehensiveBootstrap <- R6Class("comprehensiveBootstrap",
  public = list(
    bootOut = NULL, # regular boot
    bootOutExact = NULL, # secOrd boot
    Psi = NULL,
    CI_all = NULL, # list of length = 2; the first element is wald CI, the second is bootstrap CI
    initialize = function(parameter = NULL, ...) {
      # input
      # data etc...
      if(!as.character(parameter$inherit) == 'generalBootstrap') stop('please input generalBootstrap class!')
      # create two boot objects
      self$bootOut <- ateBootstrap$new(...)
      self$bootOutExact <- bootOut$clone(deep = TRUE) # deep copy point tmle, less repeat
    },
    bootstrap = function(...) {
      # input:
      # REPEAT_BOOTSTRAP
      self$bootOut$bootstrap(...)
      self$bootOutExact$exact_bootstrap(...)
      self$Psi <- self$bootOut # populate Psi_n
    },
    all_CI = function(){
      regularCI <- self$bootOut$all_boot_CI()
      taylorCI <- self$bootOutExact$all_boot_CI()
      self$CI_all <- list(wald = regularCI$wald,
                          reg = regularCI$boot,
                          reg_pen = regularCI$penalized,
                          reg_pen_half = regularCI$penalized_half,
                          reg_scale = regularCI$scale,
                          reg_scale_pen = regularCI$scale_penalized,
                          reg_scale_pen_half = regularCI$scale_penalized_half,

                          reg_ctr = regularCI$ctr,
                          reg_pen_ctr = regularCI$penalized_ctr,
                          reg_pen_half_ctr = regularCI$penalized_half_ctr,
                          reg_scale_ctr = regularCI$scale_ctr,
                          reg_scale_pen_ctr = regularCI$scale_penalized_ctr,
                          reg_scale_pen_half_ctr = regularCI$scale_penalized_half_ctr,

                          taylor = taylorCI$boot,
                          taylor_pen = taylorCI$penalized,
                          taylor_pen_half = taylorCI$penalized_half,
                          taylor_scale = taylorCI$scale,
                          taylor_scale_pen = taylorCI$scale_penalized,
                          taylor_scale_pen_half = taylorCI$scale_penalized_half,

                          taylor_ctr = taylorCI$ctr,
                          taylor_pen_ctr = taylorCI$penalized_ctr,
                          taylor_pen_half_ctr = taylorCI$penalized_half_ctr,
                          taylor_scale_ctr = taylorCI$scale_ctr,
                          taylor_scale_pen_ctr = taylorCI$scale_penalized_ctr,
                          taylor_scale_pen_half_ctr = taylorCI$scale_penalized_half_ctr
                          )
      # return(self$CI_all)
    }
  )
)