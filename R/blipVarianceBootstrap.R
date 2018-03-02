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
      self$pointTMLE <- blipVarianceTMLE_gentmle$new(data = df)
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
    }
))


