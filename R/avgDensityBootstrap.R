#' @export
avgDensityBootstrap <- R6Class("avgDensityBootstrap",
  inherit = generalBootstrap,
  public = list(
    x = NULL,
    pointTMLE = NULL,

    epsilon_step = 1e-2,
    bootstrap_estimates = NULL,
    targeting = NULL,
    initialize = function(x,
                          epsilon_step = NULL,
                          bin_width = .3,
                          lambda_grid = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
                          targeting = TRUE) {
      # bootstrap average density parameter; first do a pointTMLE;
      self$x <- x
      self$targeting <- targeting
      if(!is.null(epsilon_step)) self$epsilon_step <- epsilon_step
      onestepFit <- avgDensityTMLE$new(x = self$x, epsilon_step = self$epsilon_step, verbose  = TRUE)
      onestepFit$fit_density(bin_width = bin_width, lambda_grid = lambda_grid)
      onestepFit$calc_Psi()
      onestepFit$calc_EIC()
      if (self$targeting) onestepFit$onestepTarget()
      onestepFit$inference()

      self$pointTMLE <- onestepFit
      self$Psi <- onestepFit$Psi
    },
    bootstrap = function(REPEAT_BOOTSTRAP = 2e2, inflate_lambda = 1){
      # regular bootstrap
      SAMPLE_PER_BOOTSTRAP <- length(self$x)
      betfun <- function(data,
                        epsilon_step = self$epsilon_step,
                        inflate_lambda){
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
                                  family = stats::binomial(),
                                  inflate_lambda = inflate_lambda)
        yhat_boot <- predict.fixed_HAL(HAL_boot, new_data = d)

        yhat_boot[yhat_boot > 2*quantile(yhat_boot, probs = .75)] <- 0 # temporarily fix hal9001 extrapolation error
        density_boot <- empiricalDensity$new(p_density = yhat_boot, x = d)
        bootstrapOnestepFit$p_hat <- density_boot$normalize()
        # target new fit
        bootstrapOnestepFit$calc_Psi()
        bootstrapOnestepFit$calc_EIC()
        if (self$targeting) bootstrapOnestepFit$onestepTarget()

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
        betfun(self$x, self$epsilon_step, inflate_lambda = inflate_lambda)
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
    exact_bootstrap = function(REPEAT_BOOTSTRAP = 2e2, inflate_lambda = 1){
      # exact second order expansion bootstrap
      SAMPLE_PER_BOOTSTRAP <- length(self$x)
      betfun <- function(data,
                         epsilon_step = self$epsilon_step,
                         population_x = self$x,
                         population_tmle = self$pointTMLE,
                         inflate_lambda){
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
                                  family = stats::binomial(),
                                  inflate_lambda = inflate_lambda)
        yhat_boot <- predict.fixed_HAL(HAL_boot, new_data = d)

        # browser()
        yhat_boot[yhat_boot > 2*quantile(yhat_boot, probs = .75)] <- 0 # temporarily fix hal9001 extrapolation error
        density_boot <- empiricalDensity$new(p_density = yhat_boot, x = d)
        bootstrapOnestepFit$p_hat <- density_boot$normalize()
        # target new fit
        bootstrapOnestepFit$calc_Psi()
        bootstrapOnestepFit$calc_EIC()
        if (self$targeting) bootstrapOnestepFit$onestepTarget()

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
               population_tmle = self$pointTMLE,
               inflate_lambda = inflate_lambda)
      }
      ALPHA <- 0.05
      # remove errors
      if( !all(sapply(all_bootstrap_estimates, class) == 'numeric') ) message(paste('Error happens.', sum(sapply(all_bootstrap_estimates, class) == 'numeric'), 'bootstraps are correct'))
      all_bootstrap_estimates <- as.numeric(all_bootstrap_estimates[sapply(all_bootstrap_estimates, class) == 'numeric'])
      self$bootstrap_estimates <- all_bootstrap_estimates

      boot1_CI <- quantile(all_bootstrap_estimates, probs = c(ALPHA/2, 1 - ALPHA/2))
      normal_CI <- self$pointTMLE$CI
      self$CI_all <- list(normal_CI, boot1_CI)
    }
    # convex_bootstrap = function(REPEAT_BOOTSTRAP = 2e2, inflate_lambda = 1, alpha = 0.4){
    #   # exact second order expansion bootstrap
    #   SAMPLE_PER_BOOTSTRAP <- length(self$x)
    #   betfun <- function(data,
    #                      epsilon_step = self$epsilon_step,
    #                      population_x = self$x,
    #                      population_tmle = self$pointTMLE,
    #                      alpha = alpha,
    #                      inflate_lambda){
    #     # indices is the random indexes for the bootstrap sample
    #     indices <- sample(1:length(data), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE)
    #     indices2 <- sample(1:length(data), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # for convex combination
    #     d = (1-alpha)*data[indices] + alpha*data[indices2]

    #     bootstrapOnestepFit <- avgDensityTMLE$new(x = d, epsilon_step = epsilon_step)
    #     # fit new density
    #     longDFOut_new <- self$pointTMLE$longDataOut$generate_df_compress(x = d)
    #     HAL_boot <- fit_fixed_HAL(Y = longDFOut_new$Y,
    #                               X = longDFOut_new[,'box'],
    #                               weights = longDFOut_new$Freq, # for df_compress only
    #                               hal9001_object = self$pointTMLE$HAL_tuned,
    #                               family = stats::binomial(),
    #                               inflate_lambda = inflate_lambda)
    #     yhat_boot <- predict.fixed_HAL(HAL_boot, new_data = d)

    #     yhat_boot[yhat_boot > 2*quantile(yhat_boot, probs = .75)] <- 0 # temporarily fix hal9001 extrapolation error
    #     density_boot <- empiricalDensity$new(p_density = yhat_boot, x = d)
    #     bootstrapOnestepFit$p_hat <- density_boot$normalize()
    #     # target new fit
    #     bootstrapOnestepFit$calc_Psi()
    #     bootstrapOnestepFit$calc_EIC()
    #     if (self$targeting)  bootstrapOnestepFit$onestepTarget()

    #     return(c(bootstrapOnestepFit$Psi))
    #   }
    #   library(foreach)
    #   all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP), .combine = c,
    #                                      .inorder = FALSE,
    #                                      .packages = c('R6', 'SuperLearner'),
    #                                      # .errorhandling = 'remove',
    #                                      .errorhandling = 'pass',
    #                                      .export = c('self'),
    #                                      .verbose = F) %do% {
    #                                      # .verbose = T) %dopar% {
    #     if(it2 %% 10 == 0) print(it2)
    #     betfun(self$x, self$epsilon_step, inflate_lambda = inflate_lambda, alpha = alpha)
    #   }
    #   # save(all_bootstrap_estimates, file = 'all_bootstrap_estimates.rda')
    #   ALPHA <- 0.05
    #   # remove errors
    #   if( !all(sapply(all_bootstrap_estimates, class) == 'numeric') ) message(paste('Error happens.', sum(sapply(all_bootstrap_estimates, class) == 'numeric'), 'bootstraps are correct'))
    #   all_bootstrap_estimates <- as.numeric(all_bootstrap_estimates[sapply(all_bootstrap_estimates, class) == 'numeric'])
    #   self$bootstrap_estimates <- all_bootstrap_estimates

    #   boot1_CI <- quantile(all_bootstrap_estimates, probs = c(ALPHA/2, 1 - ALPHA/2))
    #   normal_CI <- self$pointTMLE$CI
    #   self$CI_all <- list(normal_CI, boot1_CI)
    # }
  )
)
