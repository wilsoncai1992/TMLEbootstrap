#' @export
comprehensiveBootstrap <- R6Class("comprehensiveBootstrap",
  public = list(
    bootOut = NULL, # regular boot
    bootOutExact = NULL, # secOrd boot
    Psi = NULL,
    CI_all = NULL, # list of all CI
    width_all = NULL, # the CI width in CI_all
    initialize = function(parameter = NULL, ...) {
      # input
      # data etc...
      if(!as.character(parameter$inherit) == 'generalBootstrap') stop('please input generalBootstrap class!')
      # create two boot objects
      self$bootOut <- parameter$new(...)
      self$bootOutExact <- self$bootOut$clone(deep = TRUE) # deep copy point tmle, less repeat
    },
    bootstrap = function(...) {
      # input:
      # REPEAT_BOOTSTRAP
      self$bootOut$bootstrap(...)
      self$bootOutExact$exact_bootstrap(...)
      self$Psi <- self$bootOut$Psi # populate Psi_n
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

                          secOrd = taylorCI$boot,
                          secOrd_pen = taylorCI$penalized,
                          secOrd_pen_half = taylorCI$penalized_half,
                          secOrd_scale = taylorCI$scale,
                          secOrd_scale_pen = taylorCI$scale_penalized,
                          secOrd_scale_pen_half = taylorCI$scale_penalized_half,

                          secOrd_ctr = taylorCI$ctr,
                          secOrd_pen_ctr = taylorCI$penalized_ctr,
                          secOrd_pen_half_ctr = taylorCI$penalized_half_ctr,
                          secOrd_scale_ctr = taylorCI$scale_ctr,
                          secOrd_scale_pen_ctr = taylorCI$scale_penalized_ctr,
                          secOrd_scale_pen_half_ctr = taylorCI$scale_penalized_half_ctr
                          )
      # return(self$CI_all)
    },
    compute_width = function(){
      get_width <- function(list) vapply(list, diff, FUN.VALUE = numeric(1)) # loop over list, take diff of the CI bounds
      self$width_all <- as.list(get_width(self$CI_all))
      names(self$width_all) <- names(self$CI_all)
    }
  )
)