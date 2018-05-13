library(R6)
library(hal9001)
library(foreach)
#' @export
blipVarianceBootstrap <- R6Class("blipVarianceBootstrap",
  inherit = generalBootstrap,
  public = list(
    data = NULL,
    pointTMLE = NULL,

    bootstrap_estimates = NULL,
    initialize = function(data, verbose = NULL) {
      # bootstraping blip variance TMLE (binary Y)
      self$data <- data
      if(class(data$W) != 'data.frame') message('W not data.frame')
      if (!is.null(verbose)) self$verbose <- verbose
      self$pointTMLE <- blipVarianceTMLE_gentmle$new(data = data)
      self$pointTMLE$initial_fit()
      self$pointTMLE$target()

      self$Psi <- self$pointTMLE$Psi
    },
    bootstrap = function(REPEAT_BOOTSTRAP = 2e2){
      # regular bootstrap
      SAMPLE_PER_BOOTSTRAP <- length(self$data$A)
      betfun <- function(data){
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(self$data$A), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d <- list(Y = data$Y[indices],
                  A = data$A[indices],
                  W = data.frame(data$W[indices,]))

        bootstrapTmleFit <- blipVarianceTMLE_gentmle$new(data = d)
        # fit new Q, g
        # Q fit
        Q_HAL_boot <- fit_fixed_HAL(Y = d$Y,
                                    X = data.frame(d$A, d$W),
                                    hal9001_object = self$pointTMLE$Q_fit,
                                    family = stats::binomial())
        Q_AW_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(d$A, d$W))
        Q_1W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(1, d$W))
        Q_0W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(0, d$W))
        # g fit
        g_HAL_boot <- fit_fixed_HAL(Y = d$A,
                                    X = d$W,
                                    hal9001_object = self$pointTMLE$g_fit,
                                    family = stats::binomial())
        g_1W_boot <- predict.fixed_HAL(g_HAL_boot, new_data = data.frame(d$W))
        # plug into tmle
        bootstrapTmleFit$g_1W <- g_1W_boot
        bootstrapTmleFit$Q_AW <- Q_AW_boot
        bootstrapTmleFit$Q_1W <- Q_1W_boot
        bootstrapTmleFit$Q_0W <- Q_0W_boot
        bootstrapTmleFit$Q_fit <- Q_HAL_boot
        bootstrapTmleFit$g_fit <- g_HAL_boot
        bootstrapTmleFit$target()

        return(c(bootstrapTmleFit$Psi))
      }
      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP),
                                         .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6', 'hal9001', 'fixedHAL'),
                                         # .errorhandling = 'remove',
                                         .errorhandling = 'pass',
                                         .export = c('self'),
                                         .verbose = F) %do% {
                                           # .verbose = T) %dopar% {
                                           if(it2 %% 10 == 0) print(it2)
                                           betfun(self$data)
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
      SAMPLE_PER_BOOTSTRAP <- length(self$data$A)
      betfun <- function(data, population_tmle){
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(self$data$A), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d <- list(Y = data$Y[indices],
                  A = data$A[indices],
                  W = data.frame(data$W[indices,]))

        bootstrapTmleFit <- blipVarianceTMLE_gentmle$new(data = d)
        # fit new Q, g
        # Q fit
        Q_HAL_boot <- fit_fixed_HAL(Y = d$Y,
                                    X = data.frame(d$A, d$W),
                                    hal9001_object = self$pointTMLE$Q_fit,
                                    family = stats::binomial())
        Q_AW_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(d$A, d$W))
        Q_1W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(1, d$W))
        Q_0W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(0, d$W))
        # g fit
        g_HAL_boot <- fit_fixed_HAL(Y = d$A,
                                    X = d$W,
                                    hal9001_object = self$pointTMLE$g_fit,
                                    family = stats::binomial())
        g_1W_boot <- predict.fixed_HAL(g_HAL_boot, new_data = data.frame(d$W))
        # plug into tmle
        bootstrapTmleFit$g_1W <- g_1W_boot
        bootstrapTmleFit$Q_AW <- Q_AW_boot
        bootstrapTmleFit$Q_1W <- Q_1W_boot
        bootstrapTmleFit$Q_0W <- Q_0W_boot
        bootstrapTmleFit$Q_fit <- Q_HAL_boot
        bootstrapTmleFit$g_fit <- g_HAL_boot
        bootstrapTmleFit$target()

        # predict Q#, g# on population data
        g_pound_1 <- predict.fixed_HAL(bootstrapTmleFit$g_fit, new_data = data.frame(data$W))
        Q_pound_1 <- predict.fixed_HAL(bootstrapTmleFit$Q_fit, new_data = data.frame(1, data$W))
        Q_pound_0 <- predict.fixed_HAL(bootstrapTmleFit$Q_fit, new_data = data.frame(0, data$W))
        B_pound <- Q_pound_1 - Q_pound_0
        # get R2 term
        term1 <- (population_tmle$Psi - bootstrapTmleFit$Psi)^2
        # population_tmle$gentmle_object$tmledata$Q1k
        term2 <- mean((population_tmle$Q_1W - population_tmle$Q_0W - B_pound)^2)
        cross_prod <- (population_tmle$g_1W - g_pound_1)/g_pound_1*(population_tmle$Q_1W - Q_pound_1) -
          (1-population_tmle$g_1W - (1-g_pound_1))/(1-g_pound_1)*(population_tmle$Q_0W - Q_pound_0)
        term3 <- mean(2*(B_pound - bootstrapTmleFit$Psi) * cross_prod)
        # compute R2
        R2 <- term1 - term2 + term3
        return(bootstrapTmleFit$Psi - R2)
      }
      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP),
                                         .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6', 'hal9001', 'fixedHAL'),
                                         .errorhandling = 'pass',
                                         .export = c('self'),
                                         .verbose = F) %do% {
                                         # .verbose = F) %dopar% {
                                           if(it2 %% 10 == 0) print(it2)
                                           betfun(data = self$data, population_tmle = self$pointTMLE)
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
    }
))



#' @export
blipVarianceBootstrap_contY <- R6Class("blipVarianceBootstrap_contY",
  inherit = blipVarianceBootstrap,
  public = list(
    Q_0W_rescale = NULL,
    initialize = function(data, lambda1 = NULL, lambda2 = NULL, verbose = NULL) {
      # subclass of `blipVarianceBootstrap` for continuous Y;
      # pointTMLE replaced by `blipVarianceTMLE_gentmle_contY` class
      # bootstrap method replaced by continuous hal fit; feed into `blipVarianceTMLE_gentmle_contY` class
      self$data <- data
      if(class(data$W) != 'data.frame') message('W not data.frame')
      if (!is.null(verbose)) self$verbose <- verbose
      self$pointTMLE <- blipVarianceTMLE_gentmle_contY$new(data = data)
      self$pointTMLE$scaleY()
      self$pointTMLE$initial_fit(lambda1 = lambda1, lambda2 = lambda2)
      self$pointTMLE$target()
      self$pointTMLE$scaleBack_afterTMLE()

      self$Psi <- self$pointTMLE$Psi
    },
    bootstrap = function(REPEAT_BOOTSTRAP = 2e2){
      # regular bootstrap; use continuous HAL bootstrap fit
      SAMPLE_PER_BOOTSTRAP <- length(self$data$A)
      betfun <- function(data){
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(self$data$A), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d <- list(Y = data$Y[indices],
                  A = data$A[indices],
                  W = data.frame(data$W[indices,]))

        bootstrapTmleFit <- blipVarianceTMLE_gentmle_contY$new(data = d)
        bootstrapTmleFit$scaleY() # not rigorous; use population scale_Y
        # fit new Q, g
        # Q fit
        Q_HAL_boot <- fit_fixed_HAL(Y = d$Y,
                                    X = data.frame(d$A, d$W),
                                    hal9001_object = self$pointTMLE$Q_fit,
                                    family = stats::gaussian())
        Q_AW_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(d$A, d$W))
        Q_1W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(1, d$W))
        Q_0W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(0, d$W))
        # g fit
        g_HAL_boot <- fit_fixed_HAL(Y = d$A,
                                    X = d$W,
                                    hal9001_object = self$pointTMLE$g_fit,
                                    family = stats::binomial())
        g_1W_boot <- predict.fixed_HAL(g_HAL_boot, new_data = data.frame(d$W))

        # plug into tmle
        bootstrapTmleFit$g_1W <- g_1W_boot
        bootstrapTmleFit$Q_AW_rescale <- self$pointTMLE$scale_Q$scale01(newX = Q_AW_boot) # scale to (0,1) because continuous
        bootstrapTmleFit$Q_1W_rescale <- self$pointTMLE$scale_Q$scale01(newX = Q_1W_boot)
        bootstrapTmleFit$Q_0W_rescale <- self$pointTMLE$scale_Q$scale01(newX = Q_0W_boot)
        bootstrapTmleFit$Q_fit <- Q_HAL_boot
        bootstrapTmleFit$g_fit <- g_HAL_boot
        bootstrapTmleFit$target()
        bootstrapTmleFit$scaleBack_afterTMLE() # cont Y
        return(bootstrapTmleFit$Psi)
      }
      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP),
                                         .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6', 'hal9001', 'fixedHAL'),
                                         .errorhandling = 'pass',
                                         .export = c('self'),
                                         .verbose = F) %do% {
                                         # .verbose = F) %dopar% {
                                           if(it2 %% 10 == 0) print(it2)
                                           betfun(data = self$data)
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
      # exact second order expansion of the bootstrap; use continuous HAL bootstrap fit
      SAMPLE_PER_BOOTSTRAP <- length(self$data$A)
      betfun <- function(data, population_tmle){
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(self$data$A), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d <- list(Y = data$Y[indices],
                  A = data$A[indices],
                  W = data.frame(data$W[indices,]))

        bootstrapTmleFit <- blipVarianceTMLE_gentmle_contY$new(data = d)
        bootstrapTmleFit$scaleY() # not rigorous; use population scale_Y
        # fit new Q, g
        # Q fit
        Q_HAL_boot <- fit_fixed_HAL(Y = d$Y,
                                    X = data.frame(d$A, d$W),
                                    hal9001_object = self$pointTMLE$Q_fit,
                                    family = stats::gaussian())
        Q_AW_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(d$A, d$W))
        Q_1W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(1, d$W))
        Q_0W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(0, d$W))
        # g fit
        g_HAL_boot <- fit_fixed_HAL(Y = d$A,
                                    X = d$W,
                                    hal9001_object = self$pointTMLE$g_fit,
                                    family = stats::binomial())
        g_1W_boot <- predict.fixed_HAL(g_HAL_boot, new_data = data.frame(d$W))

        # plug into tmle
        bootstrapTmleFit$g_1W <- g_1W_boot
        bootstrapTmleFit$Q_AW_rescale <- self$pointTMLE$scale_Q$scale01(newX = Q_AW_boot) # scale to (0,1) because continuous
        bootstrapTmleFit$Q_1W_rescale <- self$pointTMLE$scale_Q$scale01(newX = Q_1W_boot)
        bootstrapTmleFit$Q_0W_rescale <- self$pointTMLE$scale_Q$scale01(newX = Q_0W_boot)
        bootstrapTmleFit$Q_fit <- Q_HAL_boot
        bootstrapTmleFit$g_fit <- g_HAL_boot
        bootstrapTmleFit$target()
        bootstrapTmleFit$scaleBack_afterTMLE() # cont Y
        # predict Q#, g# on population data
        g_pound_1 <- predict.fixed_HAL(bootstrapTmleFit$g_fit, new_data = data.frame(data$W))
        Q_pound_1 <- predict.fixed_HAL(bootstrapTmleFit$Q_fit, new_data = data.frame(1, data$W))
        Q_pound_0 <- predict.fixed_HAL(bootstrapTmleFit$Q_fit, new_data = data.frame(0, data$W))
        Q_pound_1 <- self$pointTMLE$scale_Q$scale01(newX = Q_pound_1)
        Q_pound_0 <- self$pointTMLE$scale_Q$scale01(newX = Q_pound_0)
        B_pound <- Q_pound_1 - Q_pound_0
        # get R2 term
        term1 <- (population_tmle$Psi - bootstrapTmleFit$Psi)^2
        # population_tmle$gentmle_object$tmledata$Q1k
        term2 <- mean((population_tmle$Q_1W_rescale - population_tmle$Q_0W_rescale - B_pound)^2)
        cross_prod <- (population_tmle$g_1W - g_pound_1)/g_pound_1*(population_tmle$Q_1W_rescale - Q_pound_1) -
          (1-population_tmle$g_1W - (1-g_pound_1))/(1-g_pound_1)*(population_tmle$Q_0W_rescale - Q_pound_0)
        term3 <- mean(2*(B_pound - bootstrapTmleFit$Psi) * cross_prod)
        # compute R2
        R2 <- term1 - term2 + term3
        return(bootstrapTmleFit$Psi - R2)
      }
      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP),
                                         .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6', 'hal9001', 'fixedHAL'),
                                         .errorhandling = 'pass',
                                         .export = c('self'),
                                         .verbose = F) %do% {
                                         # .verbose = F) %dopar% {
                                           if(it2 %% 10 == 0) print(it2)
                                           betfun(data = self$data, population_tmle = self$pointTMLE)
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
    }
))
