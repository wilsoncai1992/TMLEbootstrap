library(R6)
#' @export
blipVarianceBootstrap <- R6Class("blipVarianceBootstrap",
  public = list(
    data = NULL,
    pointTMLE = NULL,
    Psi = NULL,

    CI_all = NULL,
    bootstrap_estimates = NULL,
    initialize = function(data, verbose = NULL) {
      self$data <- data
      if(class(data$W) != 'data.frame') message('W not data.frame')
      if (!is.null(verbose)) self$verbose <- verbose
      self$pointTMLE <- blipVarianceTMLE_gentmle$new(data = data)
      self$pointTMLE$initial_fit()
      self$pointTMLE$target()

      self$Psi <- self$pointTMLE$Psi
    },
    bootstrap = function(REPEAT_BOOTSTRAP = 2e2){
      SAMPLE_PER_BOOTSTRAP <- length(self$data$A)
      betfun <- function(data){
        # browser()
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(self$data$A), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d <- list(Y = data$Y[indices],
                  A = data$A[indices],
                  W = data.frame(data$W[indices,]))

        bootstrapTmleFit <- blipVarianceTMLE_gentmle$new(data = d)
        # fit new Q, g
        library(hal9001)
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
        bootstrapTmleFit$g_AW <- g_1W_boot
        bootstrapTmleFit$Q_AW <- Q_AW_boot
        bootstrapTmleFit$Q_1W <- Q_1W_boot
        bootstrapTmleFit$Q_0W <- Q_0W_boot
        # bootstrapTmleFit$Q_fit <- Q_HAL_boot
        # bootstrapTmleFit$g_fit <- g_HAL_boot
        bootstrapTmleFit$target()

        return(c(bootstrapTmleFit$Psi))
      }
      library(foreach)
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
    bootstrap_exact = function(REPEAT_BOOTSTRAP = 2e2){
      SAMPLE_PER_BOOTSTRAP <- length(self$data$A)
      betfun <- function(data, population_tmle){
        # browser()
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(self$data$A), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d <- list(Y = data$Y[indices],
                  A = data$A[indices],
                  W = data.frame(data$W[indices,]))

        bootstrapTmleFit <- blipVarianceTMLE_gentmle$new(data = d)
        # fit new Q, g
        library(hal9001)
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
        bootstrapTmleFit$g_AW <- g_1W_boot
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
      library(foreach)
      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP),
                                         .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6', 'hal9001', 'fixedHAL'),
                                         .errorhandling = 'pass',
                                         .export = c('self'),
                                         .verbose = F) %do% {
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


