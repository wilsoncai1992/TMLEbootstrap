# library(R6)
#' @export
generalBootstrap <- R6Class("generalBootstrap",
  public = list(
    Psi = NULL,
    CI_all = NULL, # list of length = 2; the first element is wald CI, the second is bootstrap CI
    initialize = function() {
    },
    center_CI = function(bootCI = NULL){
      # if user don't provide bootCI, use existing bootCI;
      if(is.null(bootCI)) bootCI <- self$CI_all[[2]]
      # centered bootstrap (shift 1 time)
      new_CI <- bootCI - mean(bootCI) + self$Psi
      # only shift positively
      # new_CI <- bootCI + max(0, - mean(bootCI) + self$Psi)
      return(new_CI)
    },
    scale_adjust_CI = function(bootCI = NULL){
      waldCI <- self$CI_all[[1]]
      # if user don't provide bootCI, use existing bootCI;
      if(is.null(bootCI)) bootCI <- self$CI_all[[2]]
      bootCenter <- mean(bootCI)
      r <- 1
      if(diff(bootCI) == 0) r <- 1 # catch when bootstrap Psi# are all identical
      if(diff(bootCI) < diff(waldCI)) r <- diff(bootCI)/diff(waldCI)
      # keep center the same, increase the width of the bootCI
      return((bootCI - bootCenter)/r + bootCenter)
    },
    penalized_CI = function(bootCI = NULL){
      # bias penalized bootstrap
      # if user don't provide bootCI, use existing bootCI;
      # if user input scale_adjust CI, this will output scale + penalized bootCI
      if(is.null(bootCI)) bootCI <- self$CI_all[[2]]
      delta <- mean(bootCI) - self$Psi
      bootCI[2] <- bootCI[2] + abs(delta)
      bootCI[1] <- bootCI[1] - abs(delta)
      return(bootCI)
    },
    penalized_CI_half = function(bootCI = NULL){
      # bias penalized bootstrap
      # if user don't provide bootCI, use existing bootCI;
      # if user input scale_adjust CI, this will output scale + penalized bootCI
      if(is.null(bootCI)) bootCI <- self$CI_all[[2]]
      delta <- mean(bootCI) - self$Psi
      bootCI[2] <- bootCI[2] + abs(delta)/2
      bootCI[1] <- bootCI[1] - abs(delta)/2
      return(bootCI)
    },
    shift2 = function(bootCI = NULL){
      # bias-corrected bootstrap (shift 2 times)
      new_CI <- bootCI + 2*(- mean(bootCI) + self$Psi)
      # only shift positively
      # new_CI <- bootCI + max(0, 2*(- mean(bootCI) + self$Psi))
      return(new_CI)
    },
    bias_scale = function(bootCI = NULL, n = 1e2){
      # if user don't provide bootCI, use existing bootCI;
      if(is.null(bootCI)) bootCI <- self$CI_all[[2]]
      bootCenter <- mean(bootCI)

      mse <- mean((self$bootstrap_estimates - self$Psi)^2)
      sigma_star <- sqrt(mse)
      sigma <- diff(bootCI)/1.96/2*sqrt(n) # the spread of original boot is 2*1.96*sd/sqrt(n)

      r <- 1
      if(sigma_star == 0) r <- 1 # catch when bootstrap Psi# are all identical
      if(sigma_star != 0) r <- sigma_star/sigma # always use sigma#
      # if(sigma_star > sigma) r <- sigma_star/sigma # only use sigma# to widen the CI
      # keep center the same, increase the width of the bootCI
      return((bootCI - bootCenter)*r + bootCenter)
    },
    all_boot_CI = function(){
      penalized <- self$penalized_CI()
      penalized_half <- self$penalized_CI_half()
      scale <- self$scale_adjust_CI()
      scale_penalized <- self$penalized_CI(bootCI = scale)
      scale_penalized_half <- self$penalized_CI_half(bootCI = scale)

      shift2 <- self$shift2()
      # centered versions
      center <- self$center_CI()
      penalized_center <- self$center_CI(bootCI = penalized)
      penalized_half_center <- self$center_CI(bootCI = penalized_half)
      scale_center <- self$center_CI(bootCI = scale)
      scale_penalized_center <- self$center_CI(bootCI = scale_penalized)
      scale_penalized_half_center <- self$center_CI(bootCI = scale_penalized_half)

      # bias scale
      bias_scale <- self$bias_scale()
      bias_scale_ctr <- self$center_CI(bootCI = bias_scale)
      return(list(wald = self$CI_all[[1]],

                  boot = self$CI_all[[2]],
                  penalized = penalized,
                  penalized_half = penalized_half,
                  scale = scale,
                  scale_penalized = scale_penalized,
                  scale_penalized_half = scale_penalized_half,

                  ctr = center,
                  penalized_ctr = penalized_center,
                  penalized_half_ctr = penalized_half_center,
                  scale_ctr = scale_center,
                  scale_penalized_ctr = scale_penalized_center,
                  scale_penalized_half_ctr = scale_penalized_half_center,

                  bias_scale = bias_scale,
                  bias_scale_ctr = bias_scale_ctr,
                  shift2 = shift2))
    }
    )
)