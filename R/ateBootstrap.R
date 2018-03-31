#' @export
ateBootstrap <- R6Class("ateBootstrap",
  public = list(
    data = NULL,
    pointTMLE = NULL,
    Psi = NULL,

    bootstrap_estimates = NULL,
    CI_all = NULL,
    # verbose = FALSE,
    initialize = function(data) {
      self$data <- data
      if(class(data$W) != 'data.frame') message('W not data.frame')
      tmleOut <- ateTMLE$new(data = self$data)
      tmleOut$initial_fit()
      tmleOut$target()

      self$pointTMLE <- tmleOut
      self$Psi <- self$pointTMLE$Psi
    },
    bootstrap = function(REPEAT_BOOTSTRAP = 2e2){
      # regular bootstrap
      SAMPLE_PER_BOOTSTRAP <- length(self$data$Y)
      betfun <- function(data){
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(data$Y), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        Y = data$Y[indices]
        A = data$A[indices]
        W = data$W[indices,]; W <- data.frame(W)
        d <- list(Y = Y, A = A, W = W)

        bootstrapTMLEFit <- ateTMLE$new(data = d)
        # fit new Q
        Q_boot <- fit_fixed_HAL(Y = d$Y,
                                X = data.frame(d$A, d$W),
                                hal9001_object = self$pointTMLE$Q_fit,
                                family = stats::gaussian())
        Q_1W_boot <- predict.fixed_HAL(Q_boot, new_data = data.frame(1, d$W))
        Q_0W_boot <- predict.fixed_HAL(Q_boot, new_data = data.frame(0, d$W))
        # fit new g
        g_boot <- fit_fixed_HAL(Y = d$A,
                                X = d$W,
                                hal9001_object = self$pointTMLE$g_fit,
                                family = stats::binomial())
        g1_W_boot <- predict.fixed_HAL(g_boot, new_data = data.frame(d$W)) # fixedHAL binomial does not need `plogis` transform

        bootstrapTMLEFit$Q_1W <- Q_1W_boot
        bootstrapTMLEFit$Q_0W <- Q_0W_boot
        bootstrapTMLEFit$g1_W <- g1_W_boot
        # target new fit
        bootstrapTMLEFit$target()
        # browser()
        return(c(bootstrapTMLEFit$Psi))
      }
      library(foreach)
      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP), .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6'),
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
      # exact bootstrap
      SAMPLE_PER_BOOTSTRAP <- length(self$data$Y)
      betfun <- function(data){
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(data$Y), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        Y = data$Y[indices]
        A = data$A[indices]
        W = data$W[indices,]; W <- data.frame(W)
        d <- list(Y = Y, A = A, W = W)

        bootstrapTMLEFit <- ateTMLE$new(data = d)
        # fit new Q
        Q_boot <- fit_fixed_HAL(Y = d$Y,
                                X = data.frame(d$A, d$W),
                                hal9001_object = self$pointTMLE$Q_fit,
                                family = stats::gaussian())
        Q_1W_boot <- predict.fixed_HAL(Q_boot, new_data = data.frame(1, d$W))
        Q_0W_boot <- predict.fixed_HAL(Q_boot, new_data = data.frame(0, d$W))
        # fit new g
        g_boot <- fit_fixed_HAL(Y = d$A,
                                X = d$W,
                                hal9001_object = self$pointTMLE$g_fit,
                                family = stats::binomial())
        g1_W_boot <- predict.fixed_HAL(g_boot, new_data = data.frame(d$W)) # fixedHAL binomial does not need `plogis` transform
        # plug into tmle object
        bootstrapTMLEFit$Q_1W <- Q_1W_boot
        bootstrapTMLEFit$Q_0W <- Q_0W_boot
        bootstrapTMLEFit$g1_W <- g1_W_boot
        # target new fit
        bootstrapTMLEFit$target()
        # browser()

        # get R2 term
        # predict Q#, g# on population data
        g_pound_1 <- predict.fixed_HAL(g_boot, new_data = data.frame(data$W))
        g_pound_0 <- 1-g_pound_1
        Q_pound_1 <- predict.fixed_HAL(Q_boot, new_data = data.frame(1, data$W))
        Q_pound_0 <- predict.fixed_HAL(Q_boot, new_data = data.frame(0, data$W))
        # evaluate R_2
        part1 <- (g_pound_1 - self$pointTMLE$g1_W)/g_pound_1 * (Q_pound_1 - self$pointTMLE$Q_1W)
        part0 <- (g_pound_0 - (1-self$pointTMLE$g1_W))/g_pound_0 * (Q_pound_0 - self$pointTMLE$Q_0W)
        R2 <- mean(part1 - part0)
        return(c(bootstrapTMLEFit$Psi - R2))
      }
      library(foreach)
      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP), .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6'),
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
      # bias penalized bootstrap
      bootCI <- self$CI_all[[2]]
      delta <- mean(bootCI) - self$Psi
      bootCI[2] <- bootCI[2] + abs(delta)
      bootCI[1] <- bootCI[1] - abs(delta)
      new_CI <- bootCI
      return(new_CI)
    },
    bias_corrected_boot_CI_shift1 = function(){
      # bias-corrected bootstrap (shift 1)
      new_CI <- self$CI_all[[2]] - mean(self$CI_all[[2]]) + self$Psi
      # only shift positively
      # new_CI <- self$CI_all[[2]] + max(0, - mean(self$CI_all[[2]]) + self$Psi)
      return(new_CI)
    },
    bias_corrected_boot_CI_shift2 = function(){
      # bias-corrected bootstrap (shift 2)
      new_CI <- self$CI_all[[2]] + 2*(- mean(self$CI_all[[2]]) + self$Psi)
      # only shift positively
      # new_CI <- self$CI_all[[2]] + max(0, 2*(- mean(self$CI_all[[2]]) + self$Psi))
      return(new_CI)
    },
    all_boot_CI = function(){
      # output, wald, bootstrap, bias-penalized, shift1, shift2 CI
      penalized <- self$penalized_boot_CI()
      shift1 <- self$bias_corrected_boot_CI_shift1()
      shift2 <- self$bias_corrected_boot_CI_shift2()
      return(list(normal = self$CI_all[[1]],
                  boot = self$CI_all[[2]],
                  penalized = penalized,
                  shift1 = shift1,
                  shift2 = shift2))
    }
))
