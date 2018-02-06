#' @export
avgDensityBootstrap <- R6Class("avgDensityBootstrap",
  public = list(
    x = NULL,
    pointTMLE = NULL,
    Psi = NULL,
    # EIC = NULL,
    epsilon_step = 1e-2,
    # tol = 1e-3,
    CI_all = NULL,
    bootstrap_estimates = NULL,
    initialize = function(x, epsilon_step = NULL, bin_width = .3) {
      self$x <- x
      if(!is.null(epsilon_step)) self$epsilon_step <- epsilon_step
      onestepFit <- avgDensityTMLE$new(x = self$x, epsilon_step = self$epsilon_step, verbose  = TRUE)
      onestepFit$fit_density(bin_width = bin_width)
      onestepFit$calc_Psi()
      onestepFit$calc_EIC()
      onestepFit$onestepTarget()
      onestepFit$inference()

      self$pointTMLE <- onestepFit
      self$Psi <- onestepFit$Psi
    },
    bootstrap = function(REPEAT_BOOTSTRAP = 2e2){
      SAMPLE_PER_BOOTSTRAP <- length(self$x)
      betfun <- function(data, epsilon_step = self$epsilon_step){
        # browser()
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(data), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d = data[indices]

        bootstrapOnestepFit <- avgDensityTMLE$new(x = d, epsilon_step = epsilon_step)
        # fit new density
        # longDFOut_new <- self$pointTMLE$longDataOut$generate_df(x = d)
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
      # browser()

      library(foreach)

      # library(doSNOW)
      # library(tcltk)
      # nw <- parallel:::detectCores()  # number of workers
      # cl <- makeSOCKcluster(nw)
      # registerDoSNOW(cl)

      # library(Rmpi)
      # library(doMPI)
      # cl = startMPIcluster()
      # registerDoMPI(cl)
      # clusterSize(cl) # just to check

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
      # closeCluster(cl)
      # stopCluster(cl)

      ALPHA <- 0.05
      # remove errors
      if( !all(sapply(all_bootstrap_estimates, class) == 'numeric') ) message(paste('Error happens.', sum(sapply(all_bootstrap_estimates, class) == 'numeric'), 'bootstraps are correct'))
      all_bootstrap_estimates <- as.numeric(all_bootstrap_estimates[sapply(all_bootstrap_estimates, class) == 'numeric'])
      self$bootstrap_estimates <- all_bootstrap_estimates

      boot1_CI <- quantile(all_bootstrap_estimates, probs = c(ALPHA/2, 1 - ALPHA/2))
      normal_CI <- self$pointTMLE$CI
      self$CI_all <- list(normal_CI, boot1_CI)
    },
    bootstrap_2 = function(REPEAT_BOOTSTRAP = 2e2){
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
        # longDFOut_new <- self$pointTMLE$longDataOut$generate_df(x = d)
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

        # get R2 term
        # browser()
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
        # return(c(bootstrapOnestepFit$Psi - 2*R_2))
      }

      library(foreach)
      # library(doSNOW)
      # library(tcltk)
      # nw <- parallel:::detectCores()  # number of workers
      # cl <- makeSOCKcluster(nw)
      # registerDoSNOW(cl)
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
    bootstrap_2_compensate = function(REPEAT_BOOTSTRAP = 2e2){
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
        # longDFOut_new <- self$pointTMLE$longDataOut$generate_df(x = d)
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

        # get R2 term
        # browser()
        yhat_population <- predict.fixed_HAL(HAL_boot, new_data = population_x)
        yhat_population[yhat_population > 2*quantile(yhat_population, probs = .75)] <- 0 # temporarily fix hal9001 extrapolation error
        density_population <- empiricalDensity$new(p_density = yhat_population, x = population_x)
        density_population$normalize()
        ## integration
        dummy_df <- data.frame(id = 1:length(population_x), x = population_x, p_pound = density_population$p_density, p_n = population_tmle$p_hat$p_density)
        dummy_df <- dummy_df[order(dummy_df$x),]
        dx <- c(0,diff(dummy_df$x))
        R_2 <- -sum((dummy_df$p_pound - dummy_df$p_n)^2 * dx)

        # return(c(bootstrapOnestepFit$Psi - R_2))
        return(c(bootstrapOnestepFit$Psi - 2*R_2))
      }

      library(foreach)
      # library(doSNOW)
      # library(tcltk)
      # nw <- parallel:::detectCores()  # number of workers
      # cl <- makeSOCKcluster(nw)
      # registerDoSNOW(cl)
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
    wider_boot_CI = function(){
      bootCI <- self$CI_all[[2]]
      delta <- mean(bootCI) - self$Psi
      bootCI[2] <- bootCI[2] + abs(delta)
      bootCI[1] <- bootCI[1] - abs(delta)
      new_CI <- bootCI
      return(list(self$CI_all[[1]], self$CI_all[[2]], new_CI))
    },
    center_boot_CI = function(){
      new_CI <- self$CI_all[[2]] - mean(self$CI_all[[2]]) + self$Psi
      return(list(self$CI_all[[1]], self$CI_all[[2]], new_CI))
    },
    compensate_boot_CI = function(){
      new_CI <- self$CI_all[[2]] - 2*(mean(self$CI_all[[2]]) - self$Psi)
      return(list(self$CI_all[[1]], self$CI_all[[2]], new_CI))
    }
  )
)
