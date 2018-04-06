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
    all_boot_CI = function(){
      penalized <- self$penalized_CI()
      penalized_half <- self$penalized_CI_half()
      scale <- self$scale_adjust_CI()
      scale_penalized <- self$penalized_CI(bootCI = scale)
      scale_penalized_half <- self$penalized_CI_half(bootCI = scale)

      shift2 <- self$shift2()
      # centered versions
      center <- self$center_CI()
      center_penalized <- self$penalized_CI(bootCI = center)
      center_penalized_half <- self$penalized_CI_half(bootCI = center)
      center_scale <- self$scale_adjust_CI(bootCI = center)
      center_scale_penalized <- self$penalized_CI(self$scale_adjust_CI(bootCI = center))
      center_scale_penalized_half <- self$penalized_CI_half(self$scale_adjust_CI(bootCI = center))
      return(list(wald = self$CI_all[[1]],

                  boot = self$CI_all[[2]],
                  penalized = penalized,
                  penalized_half = penalized_half,
                  scale = scale,
                  scale_penalized = scale_penalized,
                  scale_penalized_half = scale_penalized_half,

                  ctr = center,
                  ctr_penalized = center_penalized,
                  ctr_penalized_half = center_penalized_half,
                  ctr_scale = center_scale,
                  ctr_scale_penalized = center_scale_penalized,
                  ctr_scale_penalized_half = center_scale_penalized_half,
                  shift2 = shift2))
    }
    )
)