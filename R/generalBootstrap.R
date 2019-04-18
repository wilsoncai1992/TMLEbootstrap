#' @export
generalBootstrap <- R6Class("generalBootstrap",
  # bootstrap modifications
  public = list(
    Psi = NULL,
    # list of length = 2; the first element is wald CI, the second is bootstrap CI
    CI_all = NULL,
    initialize = function() {
    },
    bootstrap = function(n_bootstrap, alpha = 0.05, ...) {
      self$run_bootstrap(
        n_bootstrap = n_bootstrap, alpha = alpha, kind = "reg", ...
      )
    },
    exact_bootstrap = function(n_bootstrap, alpha = 0.05, ...) {
      self$run_bootstrap(
        n_bootstrap = n_bootstrap, alpha = alpha, kind = "sec_ord", ...
      )
    },
    exact_bootstrap_paper = function(n_bootstrap, alpha = 0.05, ...) {
      self$run_bootstrap(
        n_bootstrap = n_bootstrap, alpha = alpha, kind = "sec_ord_paper", ...
      )
    },
    center_CI = function(bootCI = NULL) {
      # if user don't provide bootCI, use existing bootCI;
      if (is.null(bootCI)) bootCI <- self$CI_all[[2]]
      # centered bootstrap (shift 1 time)
      new_CI <- bootCI - mean(bootCI) + self$Psi
      return(new_CI)
    },
    scale_adjust_CI = function(bootCI = NULL) {
      waldCI <- self$CI_all[[1]]
      # if user don't provide bootCI, use existing bootCI;
      if (is.null(bootCI)) bootCI <- self$CI_all[[2]]
      bootCenter <- mean(bootCI)
      r <- 1
      if (diff(bootCI) == 0) r <- 1 # catch when bootstrap Psi# are all identical
      if (diff(bootCI) < diff(waldCI)) r <- diff(bootCI) / diff(waldCI)
      # keep center the same, increase the width of the bootCI
      return((bootCI - bootCenter) / r + bootCenter)
    },
    penalized_CI = function(bootCI = NULL) {
      # bias penalized bootstrap
      # if user input scale_adjust CI, this will output scale + penalized bootCI

      # if user don't provide bootCI, use existing bootCI;
      if (is.null(bootCI)) bootCI <- self$CI_all[[2]]
      delta <- mean(bootCI) - self$Psi
      bootCI[2] <- bootCI[2] + abs(delta)
      bootCI[1] <- bootCI[1] - abs(delta)
      return(bootCI)
    },
    penalized_CI_half = function(bootCI = NULL) {
      # bias penalized bootstrap
      # if user input scale_adjust CI, this will output scale + penalized bootCI

      # if user don't provide bootCI, use existing bootCI;
      if (is.null(bootCI)) bootCI <- self$CI_all[[2]]
      delta <- mean(bootCI) - self$Psi
      bootCI[2] <- bootCI[2] + abs(delta) / 2
      bootCI[1] <- bootCI[1] - abs(delta) / 2
      return(bootCI)
    },
    shift2 = function(bootCI = NULL) {
      # bias-corrected bootstrap (shift 2 times)
      new_CI <- bootCI + 2 * (-mean(bootCI) + self$Psi)
      return(new_CI)
    },
    sigma_mse = function(bootCI = NULL) {
      # if user don't provide bootCI, use existing bootCI;
      if (is.null(bootCI)) bootCI <- self$CI_all[[2]]
      bootCenter <- mean(bootCI)
      mse <- mean(self$psi_bootstrap[, "reg"]^2)
      n <- nrow(self$psi_bootstrap)
      sigma_star <- sqrt(mse)

      Z_std <- scale(self$psi_bootstrap[, "reg"])
      q_z <- quantile(Z_std, probs = c(.025, .975))
      psi_n <- mean(self$CI_all[[1]])
      ci_out <- c(psi_n - q_z[2] * sigma_star, psi_n - q_z[1] * sigma_star)
      return(ci_out)
    },
    spread = function(bootCI = NULL) {
      # if user don't provide bootCI, use existing bootCI;
      if (is.null(bootCI)) bootCI <- self$CI_all[[2]]
      waldCI <- self$CI_all[[1]]
      psi_n <- mean(waldCI)
      n <- nrow(self$psi_bootstrap)
      sigma_n <- diff(waldCI) / 2 / 1.96
      bootCenter <- mean(bootCI)

      Z <- self$psi_bootstrap[, "reg"]
      Z_std <- Z / sd(Z)
      q_z <- quantile(Z_std, probs = c(.025, .975))
      ci_out <- c(psi_n - q_z[2] * sigma_n, psi_n - q_z[1] * sigma_n)
      return(ci_out)
    },
    all_boot_CI = function() {
      # return a list of all kinds of modifications
      penalized <- self$penalized_CI()
      penalized_half <- self$penalized_CI_half()
      scale <- self$scale_adjust_CI()
      scale_penalized <- self$penalized_CI(bootCI = scale)
      scale_penalized_half <- self$penalized_CI_half(bootCI = scale)

      shift2 <- self$shift2(bootCI = scale)
      # centered versions
      center <- self$center_CI()
      penalized_center <- self$center_CI(bootCI = penalized)
      penalized_half_center <- self$center_CI(bootCI = penalized_half)
      scale_center <- self$center_CI(bootCI = scale)
      scale_penalized_center <- self$center_CI(bootCI = scale_penalized)
      scale_penalized_half_center <- self$center_CI(bootCI = scale_penalized_half)

      # bias scale
      sigma_mse <- self$sigma_mse()
      sigma_mse_ctr <- self$center_CI(bootCI = sigma_mse)

      spread <- self$spread()
      return(list(
        wald = self$CI_all[[1]],

        boot = self$CI_all[[2]],
        penalized = penalized, # reg + pen
        penalized_half = penalized_half, # reg + 0.5pen
        scale = scale, # reg + scale
        scale_penalized = scale_penalized, # reg + scale + pen
        scale_penalized_half = scale_penalized_half, # reg + scale + 0.5pen

        ctr = center, # reg + center at Psi
        penalized_ctr = penalized_center, # reg + pen + center at Psi
        penalized_half_ctr = penalized_half_center, # reg + 0.5pen + center at Psi
        scale_ctr = scale_center, # reg + scale + center at Psi
        scale_penalized_ctr = scale_penalized_center, # reg + pen + scale + center at Psi
        scale_penalized_half_ctr = scale_penalized_half_center, # reg + 0.5pen + scale + center at Psi

        sigma_mse = sigma_mse, # make the width 2*1.96*sigma_star; sigma_star == sqrt(mse(Psi# - Psi_n))
        sigma_mse_ctr = sigma_mse_ctr, # sigma_mse + center
        shift2 = shift2, # compensate the bias twice; bias == |mean(Psi#) - Psi_n|

        spread = spread
      ))
    }
  )
)
