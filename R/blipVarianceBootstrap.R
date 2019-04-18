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
    targeting = NULL,
    initialize = function(data, verbose = NULL, targeting = TRUE, ...) {
      # bootstraping blip variance TMLE (binary Y)
      self$data <- data
      self$targeting <- targeting
      if (class(data$W) != "data.frame") message("W not data.frame")
      if (!is.null(verbose)) self$verbose <- verbose
      self$pointTMLE <- blipVarianceTMLE$new(data = data)
      self$pointTMLE$initial_fit(...)
      if (self$targeting) {
        self$pointTMLE$target()
      } else {
        self$pointTMLE$inference_without_target()
      }
      self$Psi <- self$pointTMLE$Psi
    },
    bootstrap_once = function(self, data, population_tmle) {
      SAMPLE_PER_BOOTSTRAP <- length(self$data$A)
      # indices is the random indexes for the bootstrap sample
      indices <- sample(
        1:length(self$data$A),
        size = SAMPLE_PER_BOOTSTRAP, replace = TRUE
      ) # user specify sample size
      d <- list(
        Y = data$Y[indices],
        A = data$A[indices],
        W = data.frame(data$W[indices, ])
      )

      bootstrapTmleFit <- blipVarianceTMLE$new(data = d)
      # fit new Q, g
      # Q fit
      Q_HAL_boot <- fit_fixed_HAL(
        Y = d$Y,
        X = data.frame(d$A, d$W),
        hal9001_object = self$pointTMLE$Q_fit,
        family = stats::binomial()
      )
      Q_AW_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(d$A, d$W))
      Q_1W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(1, d$W))
      Q_0W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(0, d$W))
      # g fit
      g_HAL_boot <- fit_fixed_HAL(
        Y = d$A,
        X = d$W,
        hal9001_object = self$pointTMLE$g_fit,
        family = stats::binomial()
      )
      g_1W_boot <- predict.fixed_HAL(g_HAL_boot, new_data = data.frame(d$W))
      # plug into tmle
      bootstrapTmleFit$g_1W <- g_1W_boot
      bootstrapTmleFit$Q_AW <- Q_AW_boot
      bootstrapTmleFit$Q_1W <- Q_1W_boot
      bootstrapTmleFit$Q_0W <- Q_0W_boot
      bootstrapTmleFit$Q_fit <- Q_HAL_boot
      bootstrapTmleFit$g_fit <- g_HAL_boot
      if (self$targeting) {
        bootstrapTmleFit$target()
      } else {
        bootstrapTmleFit$inference_without_target()
      }

      # predict Q#, g# on population data
      g_pound_1 <- predict.fixed_HAL(bootstrapTmleFit$g_fit, new_data = data.frame(data$W))
      Q_pound_1 <- predict.fixed_HAL(bootstrapTmleFit$Q_fit, new_data = data.frame(1, data$W))
      Q_pound_0 <- predict.fixed_HAL(bootstrapTmleFit$Q_fit, new_data = data.frame(0, data$W))
      B_pound <- Q_pound_1 - Q_pound_0
      # get R2 term
      term1 <- (population_tmle$Psi - bootstrapTmleFit$Psi)^2
      # population_tmle$gentmle_object$tmledata$Q1k
      term2 <- mean((population_tmle$Q_1W - population_tmle$Q_0W - B_pound)^2)
      cross_prod <- (population_tmle$g_1W - g_pound_1) / g_pound_1 * (population_tmle$Q_1W - Q_pound_1) -
        (1 - population_tmle$g_1W - (1 - g_pound_1)) / (1 - g_pound_1) * (population_tmle$Q_0W - Q_pound_0)
      term3 <- mean(2 * (B_pound - bootstrapTmleFit$Psi) * cross_prod)
      # compute R2
      R2 <- term1 - term2 + term3
      return(data.frame(
        reg = bootstrapTmleFit$Psi - self$pointTMLE$Psi,
        sec_ord = bootstrapTmleFit$Psi - R2 - self$pointTMLE$Psi
      ))
    },
    run_bootstrap = function(n_bootstrap = 2e2, alpha = 0.05, kind = NULL) {
      # all bootstrap
      if (is.null(self$bootstrap_estimates)) {
        library(foreach)
        all_bootstrap_estimates <- foreach(
          it2 = 1:n_bootstrap,
          .combine = "rbind",
          .inorder = FALSE,
          .errorhandling = "remove",
          .export = c("self")
        ) %do% {
          self$bootstrap_once(
            self = self, data = self$data, population_tmle = self$pointTMLE
          )
        }
        self$bootstrap_estimates <- all_bootstrap_estimates
        dim(self$bootstrap_estimates)
      } else {
        message("invoke cached bootstrap results")
      }
      Z_quantile <- quantile(
        self$bootstrap_estimates[, kind],
        probs = c(alpha / 2, 1 - alpha / 2)
      )
      normal_CI <- self$pointTMLE$CI
      boot1_CI <- c(
        self$pointTMLE$Psi - Z_quantile[2], self$pointTMLE$Psi - Z_quantile[1]
      )
      self$CI_all <- list(normal_CI, boot1_CI)
    }
  )
)



#' @export
blipVarianceBootstrapContinuousY <- R6Class("blipVarianceBootstrapContinuousY",
  inherit = blipVarianceBootstrap,
  public = list(
    Q_0W_rescale = NULL,
    initialize = function(
                              data,
                              lambda1 = NULL,
                              lambda2 = NULL,
                              M1 = NULL,
                              M2 = NULL,
                              verbose = NULL,
                              targeting = TRUE) {
      # subclass of `blipVarianceBootstrap` for continuous Y;
      # pointTMLE replaced by `blipVarianceTMLEContinuousY` class
      # bootstrap method replaced by continuous hal fit;
      # feed into `blipVarianceTMLEContinuousY` class
      self$data <- data
      self$targeting <- targeting
      if (class(data$W) != "data.frame") message("W not data.frame")
      if (!is.null(verbose)) self$verbose <- verbose
      self$pointTMLE <- blipVarianceTMLEContinuousY$new(data = data)
      self$pointTMLE$scale_outcome()
      self$pointTMLE$initial_fit(lambda1 = lambda1, lambda2 = lambda2, M1 = M1, M2 = M2)
      if (self$targeting) {
        self$pointTMLE$target()
      } else {
        self$pointTMLE$inference_without_target()
      }
      self$pointTMLE$scale_back_after_tmle()

      self$Psi <- self$pointTMLE$Psi
    },
    bootstrap_once = function(self, data, population_tmle) {
      SAMPLE_PER_BOOTSTRAP <- length(self$data$A)
      # indices is the random indexes for the bootstrap sample
      indices <- sample(1:SAMPLE_PER_BOOTSTRAP, size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
      d <- list(
        Y = data$Y[indices],
        A = data$A[indices],
        W = data.frame(data$W[indices, ])
      )
      bootstrapTmleFit <- blipVarianceTMLEContinuousY$new(data = d)
      # use population scale_Y
      bootstrapTmleFit$scale_Y <- self$pointTMLE$scale_Y
      bootstrapTmleFit$Y_rescale <- bootstrapTmleFit$scale_Y$scale01(newX = bootstrapTmleFit$data$Y)
      # fit new Q, g
      # Q fit
      Q_HAL_boot <- fit_fixed_HAL(
        Y = d$Y,
        X = data.frame(d$A, d$W),
        hal9001_object = self$pointTMLE$Q_fit,
        family = stats::gaussian()
      )
      Q_AW_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(d$A, d$W))
      Q_1W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(1, d$W))
      Q_0W_boot <- predict.fixed_HAL(Q_HAL_boot, new_data = data.frame(0, d$W))
      # g fit
      g_HAL_boot <- fit_fixed_HAL(
        Y = d$A,
        X = d$W,
        hal9001_object = self$pointTMLE$g_fit,
        family = stats::binomial()
      )
      g_1W_boot <- predict.fixed_HAL(g_HAL_boot, new_data = data.frame(d$W))

      # plug into tmle
      bootstrapTmleFit$g_1W <- g_1W_boot
      # scale to (0,1) because continuous
      bootstrapTmleFit$Q_AW_rescale <- self$pointTMLE$scale_Y$scale01(
        newX = Q_AW_boot
      )
      bootstrapTmleFit$Q_1W_rescale <- self$pointTMLE$scale_Y$scale01(
        newX = Q_1W_boot
      )
      bootstrapTmleFit$Q_0W_rescale <- self$pointTMLE$scale_Y$scale01(
        newX = Q_0W_boot
      )
      bootstrapTmleFit$Q_fit <- Q_HAL_boot
      bootstrapTmleFit$g_fit <- g_HAL_boot
      if (self$targeting) {
        bootstrapTmleFit$target()
      } else {
        bootstrapTmleFit$inference_without_target()
      }
      bootstrapTmleFit$scale_back_after_tmle() # cont Y
      # PnD*
      # the nuisance parameters need to be for continuous Y \in R
      PnDstar <- mean(self$pointTMLE$compute_EIC(
        A = d$A,
        gk = g_1W_boot,
        Y = d$Y,
        Qk = Q_AW_boot,
        Q1k = Q_1W_boot,
        Q0k = Q_0W_boot,
        psi = var(Q_1W_boot - Q_0W_boot)
      ))
      # predict Q#, g# on population data
      g_pound_1 <- predict.fixed_HAL(bootstrapTmleFit$g_fit, new_data = data.frame(data$W))
      Q_pound_1 <- predict.fixed_HAL(bootstrapTmleFit$Q_fit, new_data = data.frame(1, data$W))
      Q_pound_0 <- predict.fixed_HAL(bootstrapTmleFit$Q_fit, new_data = data.frame(0, data$W))
      # Q_pound_1 <- self$pointTMLE$scale_Y$scale01(newX = Q_pound_1)
      # Q_pound_0 <- self$pointTMLE$scale_Y$scale01(newX = Q_pound_0)
      Q_pound_A <- data$A * Q_pound_1 + (1 - data$A) * Q_pound_0
      B_pound <- Q_pound_1 - Q_pound_0
      # P0D*
      P0Dstar <- mean(self$pointTMLE$compute_EIC(
        A = data$A,
        gk = g_pound_1,
        Y = data$Y,
        Qk = Q_pound_A,
        Q1k = Q_pound_1,
        Q0k = Q_pound_0,
        psi = var(Q_pound_1 - Q_pound_0)
      ))
      # get R2 term
      term1 <- (population_tmle$Psi - bootstrapTmleFit$Psi)^2
      term2 <- mean(
        (population_tmle$Q_1W - population_tmle$Q_0W - B_pound)^2
      )
      cross_prod <- (population_tmle$g_1W - g_pound_1) / g_pound_1 *
        (population_tmle$Q_1W - Q_pound_1) -
        (1 - population_tmle$g_1W - (1 - g_pound_1)) / (1 - g_pound_1) *
          (population_tmle$Q_0W - Q_pound_0)
      term3 <- mean(2 * (B_pound - bootstrapTmleFit$Psi) * cross_prod)
      # compute R2
      R2 <- term1 - term2 + term3
      return(data.frame(
        reg = bootstrapTmleFit$Psi - self$pointTMLE$Psi,
        sec_ord = bootstrapTmleFit$Psi - R2 - self$pointTMLE$Psi,
        sec_ord_paper = PnDstar - P0Dstar + R2
      ))
    },
    run_bootstrap = function(n_bootstrap = 2e2, alpha = 0.05, kind = NULL) {
      # all bootstrap
      if (is.null(self$bootstrap_estimates)) {
        library(foreach)
        all_bootstrap_estimates <- foreach(
          it2 = 1:n_bootstrap,
          .combine = "rbind",
          .inorder = FALSE,
          .errorhandling = "remove",
          .export = c("self")
        ) %do% {
          self$bootstrap_once(
            self = self, data = self$data, population_tmle = self$pointTMLE
          )
        }
        self$bootstrap_estimates <- all_bootstrap_estimates
        dim(self$bootstrap_estimates)
      } else {
        message("invoke cached bootstrap results")
      }
      Z_quantile <- quantile(
        self$bootstrap_estimates[, kind], probs = c(alpha / 2, 1 - alpha / 2)
      )
      normal_CI <- self$pointTMLE$CI
      boot1_CI <- c(
        self$pointTMLE$Psi - Z_quantile[2], self$pointTMLE$Psi - Z_quantile[1]
      )
      self$CI_all <- list(normal_CI, boot1_CI)
    }
  )
)
