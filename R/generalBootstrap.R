# library(R6)
#' @export
generalBootstrap <- R6Class("generalBootstrap",
  public = list(
    Psi = NULL,
    CI_all = NULL, # list of length = 2; the first element is wald CI, the second is bootstrap CI
    initialize = function() {
    },
    scale_adjust_boot_CI = function(bootCI = NULL){
      waldCI <- self$CI_all[[1]]
      # if user don't provide bootCI, use existing bootCI;
      if(is.null(bootCI)) bootCI <- self$CI_all[[2]]
      bootCenter <- mean(bootCI)
      r <- 1
      if(diff(bootCI) == 0) r <- 1 # catch when bootstrap Psi# are all identical
      if(diff(bootCI) < diff(waldCI)) r <- diff(bootCI)/diff(waldCI)
      # keep center the same, increase the width of the bootCI
      newCI <- (bootCI - bootCenter)/r + bootCenter
      return(newCI)
    },
    penalized_boot_CI = function(bootCI = NULL){
      # bias penalized bootstrap
      # if user don't provide bootCI, use existing bootCI;
      # if user input scale_adjust CI, this will output scale + penalized bootCI
      if(is.null(bootCI)) bootCI <- self$CI_all[[2]]
      delta <- mean(bootCI) - self$Psi
      bootCI[2] <- bootCI[2] + abs(delta)
      bootCI[1] <- bootCI[1] - abs(delta)
      new_CI <- bootCI
      return(new_CI)
    },
    bias_corrected_boot_CI_shift1 = function(){
      # bias-corrected bootstrap (shift 1 time)
      new_CI <- self$CI_all[[2]] - mean(self$CI_all[[2]]) + self$Psi
      # only shift positively
      # new_CI <- self$CI_all[[2]] + max(0, - mean(self$CI_all[[2]]) + self$Psi)
      return(new_CI)
    },
    bias_corrected_boot_CI_shift2 = function(){
      # bias-corrected bootstrap (shift 2 times)
      new_CI <- self$CI_all[[2]] + 2*(- mean(self$CI_all[[2]]) + self$Psi)
      # only shift positively
      # new_CI <- self$CI_all[[2]] + max(0, 2*(- mean(self$CI_all[[2]]) + self$Psi))
      return(new_CI)
    },
    all_boot_CI = function(){
      # output, wald, bootstrap, bias-penalized, shift1, shift2 CI
      penalized <- self$penalized_boot_CI()
      scale <- self$scale_adjust_boot_CI()
      scale_penalized <- self$penalized_boot_CI(bootCI = scale)
      shift1 <- self$bias_corrected_boot_CI_shift1()
      shift2 <- self$bias_corrected_boot_CI_shift2()
      return(list(wald = self$CI_all[[1]],
                  boot = self$CI_all[[2]],
                  penalized = penalized,
                  scale = scale,
                  scale_penalized = scale_penalized,
                  shift1 = shift1,
                  shift2 = shift2))
    }

    )
)