#' @export
avgDensityBootstrap <- R6Class("avgDensityBootstrap",
  public = list(
    x = NULL,
    pointTMLE = NULL,
    Psi = NULL,
    epsilon_step = 1e-2,
    CI_all = NULL,
    bootstrap_estimates = NULL,
    initialize = function(x, epsilon_step = NULL, bin_width = .3, lambda_grid = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) {
      # bootstrap average density parameter; first do a pointTMLE;
      self$x <- x
      if(!is.null(epsilon_step)) self$epsilon_step <- epsilon_step
      onestepFit <- avgDensityTMLE$new(x = self$x, epsilon_step = self$epsilon_step, verbose  = TRUE)
      onestepFit$fit_density(bin_width = bin_width, lambda_grid = lambda_grid)
      onestepFit$calc_Psi()
      onestepFit$calc_EIC()
      onestepFit$onestepTarget()
      onestepFit$inference()

      self$pointTMLE <- onestepFit
      self$Psi <- onestepFit$Psi
    },
    bootstrap = function(REPEAT_BOOTSTRAP = 2e2){
      # regular bootstrap
      SAMPLE_PER_BOOTSTRAP <- length(self$x)
      betfun <- function(data, epsilon_step = self$epsilon_step){
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(data), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d = data[indices]

        bootstrapOnestepFit <- avgDensityTMLE$new(x = d, epsilon_step = epsilon_step)
        # fit new density
        longDFOut_new <- self$pointTMLE$longDataOut$generate_df_compress(x = d)
        HAL_boot <- fit_fixed_HAL(Y = longDFOut_new$Y,
          X = longDFOut_new[,'box'],
          weights = longDFOut_new$Freq, # for df_compress only
          hal9001_object = self$pointTMLE$HAL_tuned,
          family = stats::binomial())
        yhat_boot <- predict.fixed_HAL(HAL_boot, new_data = d)

        yhat_boot[yhat_boot > 2*quantile(yhat_boot, probs = .75)] <- 0 # temporarily fix hal9001 extrapolation error
        density_boot <- empiricalDensity$new(p_density = yhat_boot, x = d)
        bootstrapOnestepFit$p_hat <- density_boot$normalize()
        # target new fit
        bootstrapOnestepFit$calc_Psi()
        bootstrapOnestepFit$calc_EIC()
        bootstrapOnestepFit$onestepTarget()

        return(c(bootstrapOnestepFit$Psi))
      }
      library(foreach)
      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP), .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6', 'SuperLearner'),
                                         # .errorhandling = 'remove',
                                         .errorhandling = 'pass',
                                         .export = c('self'),
                                         .verbose = F) %do% {
                                         # .verbose = T) %dopar% {
        if(it2 %% 10 == 0) print(it2)
        betfun(self$x, self$epsilon_step)
      }
      # save(all_bootstrap_estimates, file = 'all_bootstrap_estimates.rda')
      ALPHA <- 0.05
      # remove errors
      if( !all(sapply(all_bootstrap_estimates, class) == 'numeric') ) message(paste('Error happens.', sum(sapply(all_bootstrap_estimates, class) == 'numeric'), 'bootstraps are correct'))
      all_bootstrap_estimates <- as.numeric(all_bootstrap_estimates[sapply(all_bootstrap_estimates, class) == 'numeric'])
      self$bootstrap_estimates <- all_bootstrap_estimates

      boot1_CI <- quantile(all_bootstrap_estimates, probs = c(ALPHA/2, 1 - ALPHA/2))
      normal_CI <- self$pointTMLE$CI
      self$CI_all <- list(normal_CI, boot1_CI)
    },
    exact_bootstrap = function(REPEAT_BOOTSTRAP = 2e2){
      # exact second order expansion bootstrap
      SAMPLE_PER_BOOTSTRAP <- length(self$x)
      betfun <- function(data,
                         epsilon_step = self$epsilon_step,
                         population_x = self$x,
                         population_tmle = self$pointTMLE){
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(data), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d = data[indices]

        bootstrapOnestepFit <- avgDensityTMLE$new(x = d, epsilon_step = epsilon_step)
        # fit new density
        longDFOut_new <- self$pointTMLE$longDataOut$generate_df_compress(x = d)
        HAL_boot <- fit_fixed_HAL(Y = longDFOut_new$Y,
          X = longDFOut_new[,'box'],
          weights = longDFOut_new$Freq, # for df_compress only
          hal9001_object = self$pointTMLE$HAL_tuned,
          family = stats::binomial())
        yhat_boot <- predict.fixed_HAL(HAL_boot, new_data = d)

        # browser()
        yhat_boot[yhat_boot > 2*quantile(yhat_boot, probs = .75)] <- 0 # temporarily fix hal9001 extrapolation error
        density_boot <- empiricalDensity$new(p_density = yhat_boot, x = d)
        bootstrapOnestepFit$p_hat <- density_boot$normalize()
        # target new fit
        bootstrapOnestepFit$calc_Psi()
        bootstrapOnestepFit$calc_EIC()
        bootstrapOnestepFit$onestepTarget()

        # get R2 term
        yhat_population <- predict.fixed_HAL(HAL_boot, new_data = population_x)
        yhat_population[yhat_population > 2*quantile(yhat_population, probs = .75)] <- 0 # temporarily fix hal9001 extrapolation error
        density_population <- empiricalDensity$new(p_density = yhat_population, x = population_x)
        density_population$normalize()
        ## integration
        dummy_df <- data.frame(id = 1:length(population_x), x = population_x, p_pound = density_population$p_density, p_n = population_tmle$p_hat$p_density)
        dummy_df <- dummy_df[order(dummy_df$x),]
        dx <- c(0,diff(dummy_df$x))
        R_2 <- -sum((dummy_df$p_pound - dummy_df$p_n)^2 * dx)

        return(c(bootstrapOnestepFit$Psi - R_2))
      }

      library(foreach)
      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP), .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6', 'SuperLearner'),
                                         # .errorhandling = 'remove',
                                         .errorhandling = 'pass',
                                         .export = c('self'),
                                         .verbose = F) %do% {
                                         # .verbose = T) %dopar% {
        if(it2 %% 10 == 0) print(it2)
        betfun(self$x,
               epsilon_step = self$epsilon_step,
               population_x = self$x,
               population_tmle = self$pointTMLE)
      }
      ALPHA <- 0.05
      # remove errors
      if( !all(sapply(all_bootstrap_estimates, class) == 'numeric') ) message(paste('Error happens.', sum(sapply(all_bootstrap_estimates, class) == 'numeric'), 'bootstraps are correct'))
      all_bootstrap_estimates <- as.numeric(all_bootstrap_estimates[sapply(all_bootstrap_estimates, class) == 'numeric'])
      self$bootstrap_estimates <- all_bootstrap_estimates

      boot1_CI <- quantile(all_bootstrap_estimates, probs = c(ALPHA/2, 1 - ALPHA/2))
      normal_CI <- self$pointTMLE$CI
      self$CI_all <- list(normal_CI, boot1_CI)
    },
    exact_bootstrap_widthadjusted = function(REPEAT_BOOTSTRAP = 2e2){
      # make width of bootstrap CI at least as wide as wald CI; bootstrap CI center don't shift
      self$exact_bootstrap(REPEAT_BOOTSTRAP = REPEAT_BOOTSTRAP)
      waldCI <- self$CI_all[[1]]
      bootCI <- self$CI_all[[2]]
      bootCenter <- mean(bootCI)
      r <- 1
      if(diff(bootCI) == 0) r <- 1 # catch when bootstrap Psi# are all identical
      if(diff(bootCI) < diff(waldCI)) r <- diff(bootCI)/diff(waldCI)
      # keep center the same, increase the width of the bootCI
      bootCI <- (bootCI - bootCenter)/r + bootCenter

      self$CI_all <- list(waldCI, bootCI)
    },
    penalized_boot_CI = function(){
      # bias-penalized bootstrap
      bootCI <- self$CI_all[[2]]
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
      # output wald, bootstrap, penalized, shift1, shift2 intervals
      penalized <- self$penalized_boot_CI()
      shift1 <- self$bias_corrected_boot_CI_shift1()
      shift2 <- self$bias_corrected_boot_CI_shift2()
      return(list(normal = self$CI_all[[1]],
                  boot = self$CI_all[[2]],
                  penalized = penalized,
                  shift1 = shift1,
                  shift2 = shift2))
    }
  )
)
