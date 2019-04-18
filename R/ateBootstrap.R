#' @export
ateBootstrap <- R6Class("ateBootstrap",
  inherit = generalBootstrap,
  public = list(
    data = NULL,
    lambda1 = NULL,
    M1 = NULL,
    pointTMLE = NULL,

    bootstrap_estimates = NULL,
    targeting = NULL,
    initialize = function(data,
                              lambda1 = NULL,
                              lambda2 = NULL,
                              M1 = NULL,
                              M2 = NULL,
                              targeting = TRUE) {
      # data is in list
      # lambda1 is a grid of lambda for Q
      # lambda2 is a grid of lambda for g
      self$data <- data
      self$lambda1 <- lambda1
      self$targeting <- targeting
      if (class(data$W) != "data.frame") message("W not data.frame")
      tmleOut <- ateTMLE$new(data = self$data)
      tmleOut$initial_fit(lambda1 = lambda1, lambda2 = lambda2, M1 = M1, M2 = M2)
      if (self$targeting) {
        tmleOut$target()
      } else {
        tmleOut$inference_without_target(
          data = self$data, tmleOut$Q_fit, tmleOut$g_fit, to_return = FALSE
        )
      }

      self$pointTMLE <- tmleOut
      self$Psi <- self$pointTMLE$Psi
    },
    bootstrap_once = function(self, data) {
      SAMPLE_PER_BOOTSTRAP <- length(self$data$Y)
      # indices is the random indexes for the bootstrap sample
      indices <- sample(
        1:SAMPLE_PER_BOOTSTRAP,
        size = SAMPLE_PER_BOOTSTRAP, replace = TRUE
      )
      Y <- data$Y[indices]
      A <- data$A[indices]
      W <- data$W[indices, ]
      W <- data.frame(W)
      d <- list(Y = Y, A = A, W = W)

      bootstrapTMLEFit <- ateTMLE$new(data = d)
      # fit new Q
      Q_boot <- fit_fixed_HAL(
        Y = d$Y,
        X = data.frame(d$A, d$W),
        hal9001_object = self$pointTMLE$Q_fit,
        family = stats::gaussian()
      )
      Q_1W_boot <- predict.fixed_HAL(Q_boot, new_data = data.frame(1, d$W))
      Q_0W_boot <- predict.fixed_HAL(Q_boot, new_data = data.frame(0, d$W))
      Q_AW_boot <- predict.fixed_HAL(Q_boot, new_data = data.frame(d$A, d$W))
      # fit new g
      g_boot <- fit_fixed_HAL(
        Y = d$A,
        X = d$W,
        hal9001_object = self$pointTMLE$g_fit,
        family = stats::binomial()
      )
      g1_W_boot <- predict.fixed_HAL(g_boot, new_data = data.frame(d$W))
      # plug into tmle object
      bootstrapTMLEFit$Q_1W <- Q_1W_boot
      bootstrapTMLEFit$Q_0W <- Q_0W_boot
      bootstrapTMLEFit$g1_W <- g1_W_boot
      # target new fit
      if (self$targeting) {
        bootstrapTMLEFit$target()
      } else {
        bootstrapTMLEFit$inference_without_target(data = d, NULL, NULL, to_return = FALSE)
      }
      # PnD*
      PnDstar <- mean(self$pointTMLE$compute_EIC(
        A = d$A,
        gk = g1_W_boot,
        Y = d$Y,
        Qk = Q_AW_boot,
        Q1k = Q_1W_boot,
        Q0k = Q_0W_boot,
        psi = mean(Q_1W_boot - Q_0W_boot)
      ))
      # predict Q#, g# on population data
      g_pound_1 <- predict.fixed_HAL(g_boot, new_data = data.frame(data$W))
      g_pound_0 <- 1 - g_pound_1
      Q_pound_1 <- predict.fixed_HAL(Q_boot, new_data = data.frame(1, data$W))
      Q_pound_0 <- predict.fixed_HAL(Q_boot, new_data = data.frame(0, data$W))
      Q_pound_A <- predict.fixed_HAL(Q_boot, new_data = data.frame(data$A, data$W))
      # P0D*
      P0Dstar <- mean(self$pointTMLE$compute_EIC(
        A = data$A,
        gk = g_pound_1,
        Y = data$Y,
        Qk = Q_pound_A,
        Q1k = Q_pound_1,
        Q0k = Q_pound_0,
        psi = mean(Q_pound_1 - Q_pound_0)
      ))
      # get R2 term
      part1 <- (g_pound_1 - self$pointTMLE$g1_W) / g_pound_1 *
        (Q_pound_1 - self$pointTMLE$Q_1W)
      part0 <- (g_pound_0 - (1 - self$pointTMLE$g1_W)) / g_pound_0 *
        (Q_pound_0 - self$pointTMLE$Q_0W)
      R2 <- mean(part1 - part0)
      return(data.frame(
        reg = bootstrapTMLEFit$Psi - self$pointTMLE$Psi,
        sec_ord = bootstrapTMLEFit$Psi - R2 - self$pointTMLE$Psi,
        sec_ord_paper = PnDstar - P0Dstar + R2
      ))
    },
    run_bootstrap = function(n_bootstrap = 2e2, alpha = 0.05, kind = NULL) {
      # all bootstrap
      library(foreach)
      if (is.null(self$bootstrap_estimates)) {
        all_bootstrap_estimates <- foreach(
          it2 = 1:n_bootstrap,
          .combine = "rbind",
          .inorder = FALSE,
          .errorhandling = "remove",
          .export = c("self")
        ) %do% {
          self$bootstrap_once(self = self, data = self$data)
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
