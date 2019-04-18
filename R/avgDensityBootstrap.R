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
                              lambda_grid = NULL,
                              M = NULL,
                              targeting = TRUE) {
      # bootstrap average density parameter; first do a pointTMLE;
      self$x <- x
      self$targeting <- targeting
      if (!is.null(epsilon_step)) self$epsilon_step <- epsilon_step
      onestepFit <- avgDensityTMLE$new(
        x = self$x, epsilon_step = self$epsilon_step, verbose = TRUE
      )
      onestepFit$fit_density(
        bin_width = bin_width, lambda_grid = lambda_grid, M = M, n_fold = 3
      )
      onestepFit$compute_Psi(p_hat = onestepFit$p_hat, to_return = FALSE)
      onestepFit$compute_EIC(
        p_hat = onestepFit$p_hat, Psi = onestepFit$Psi, to_return = FALSE
      )
      if (self$targeting) onestepFit$onestepTarget()
      onestepFit$inference()

      self$pointTMLE <- onestepFit
      self$Psi <- onestepFit$Psi
    },
    bootstrap_once = function(self,
                                  data,
                                  epsilon_step = self$epsilon_step,
                                  population_x = self$x,
                                  population_tmle = self$pointTMLE,
                                  inflate_lambda = 1) {
      SAMPLE_PER_BOOTSTRAP <- length(self$x)
      # indices is the random indexes for the bootstrap sample
      indices <- sample(1:length(data), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE)
      d <- data[indices]

      bootstrapOnestepFit <- avgDensityTMLE$new(
        x = d, epsilon_step = epsilon_step, verbose = FALSE
      )
      # fit new density
      longDFOut_new <- self$pointTMLE$longDataOut$generate_df_compress(x = d)
      HAL_boot <- fit_fixed_HAL(
        Y = longDFOut_new$Y,
        X = longDFOut_new[, "box"],
        weights = longDFOut_new$Freq, # for df_compress only
        hal9001_object = self$pointTMLE$HAL_tuned,
        family = stats::binomial(),
        inflate_lambda = inflate_lambda
      )
      yhat_boot <- predict.fixed_HAL(HAL_boot, new_data = d)

      # temporarily fix hal9001 extrapolation error
      yhat_boot[yhat_boot > 2 * quantile(yhat_boot, probs = .75)] <- 0
      density_boot <- empiricalDensity$new(p_density = yhat_boot, x = d)
      bootstrapOnestepFit$p_hat <- density_boot$normalize()
      # target new fit
      bootstrapOnestepFit$compute_Psi(p_hat = bootstrapOnestepFit$p_hat, FALSE)
      bootstrapOnestepFit$compute_EIC(p_hat = bootstrapOnestepFit$p_hat, Psi = bootstrapOnestepFit$Psi, FALSE)
      if (self$targeting) bootstrapOnestepFit$onestepTarget()

      # get R2 term
      yhat_population <- predict.fixed_HAL(HAL_boot, new_data = population_x)
      yhat_population[yhat_population > 2 * quantile(yhat_population, probs = .75)] <- 0 # temporarily fix hal9001 extrapolation error
      density_population <- empiricalDensity$new(p_density = yhat_population, x = population_x)
      density_population$normalize()
      ## integration
      dummy_df <- data.frame(
        id = 1:length(population_x),
        x = population_x,
        p_pound = density_population$p_density,
        p_n = population_tmle$p_hat$p_density
      )
      dummy_df <- dummy_df[order(dummy_df$x), ]
      dx <- c(0, diff(dummy_df$x))
      R_2 <- -sum((dummy_df$p_pound - dummy_df$p_n)^2 * dx)

      dummy_df_d <- data.frame(
        id = 1:length(density_boot$x),
        x = density_boot$x,
        p_density = density_boot$p_density
      )
      dummy_df_d <- dummy_df_d[order(dummy_df_d$x), ]
      dx_d <- c(0, diff(dummy_df_d$x))
      PnDstar <- 2 * mean(density_boot$p_density - sum(dummy_df_d$p_density^2 * dx_d))
      P0Dstar <- 2 * mean(dummy_df$p_pound - sum(dummy_df$p_pound^2 * dx))
      return(data.frame(
        reg = bootstrapOnestepFit$Psi - self$pointTMLE$Psi,
        sec_ord = bootstrapOnestepFit$Psi - R_2 - self$pointTMLE$Psi,
        sec_ord_paper = PnDstar - P0Dstar + R_2
      ))
    },
    run_bootstrap = function(
                                 REPEAT_BOOTSTRAP = 2e2,
                                 ALPHA = 0.05,
                                 kind = NULL,
                                 inflate_lambda = 1,
                                 to_parallel = FALSE) {
      # all bootstrap
      library(foreach)
      `%mydo%` <- ifelse(to_parallel, `%dopar%`, `%do%`)
      if (to_parallel) {
        library(doSNOW)
        library(tcltk)
        nw <- parallel:::detectCores() # number of workers
        cl <- makeSOCKcluster(nw)
        registerDoSNOW(cl)
      }
      if (is.null(self$bootstrap_estimates)) {
        all_bootstrap_estimates <- foreach(
          it2 = 1:REPEAT_BOOTSTRAP,
          .combine = "rbind",
          .inorder = FALSE,
          .errorhandling = "remove",
          .export = c("self")
        ) %mydo% {
          self$bootstrap_once(
            self = self,
            data = self$x,
            epsilon_step = self$epsilon_step,
            population_x = self$x,
            population_tmle = self$pointTMLE,
            inflate_lambda = inflate_lambda
          )
        }
        self$bootstrap_estimates <- all_bootstrap_estimates
        dim(self$bootstrap_estimates)
      } else {
        message("invoke cached bootstrap results")
      }
      if (to_parallel) stopCluster(cl)
      Z_quantile <- quantile(
        self$bootstrap_estimates[, kind],
        probs = c(ALPHA / 2, 1 - ALPHA / 2)
      )
      normal_CI <- self$pointTMLE$CI
      boot1_CI <- c(
        self$pointTMLE$Psi - Z_quantile[2], self$pointTMLE$Psi - Z_quantile[1]
      )
      self$CI_all <- list(normal_CI, boot1_CI)
    }
  )
)
