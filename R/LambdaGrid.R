#' @export
LambdaGrid <- R6Class("LambdaGrid",
  # lambda plateau method
  public = list(
    data = NULL,
    REPEAT_BOOTSTRAP = NULL,
    inflate_lambda = NULL,

    dict_boot = list(), # named list of comprehensiveBootstrap
    df_lambda_width = NULL,

    # OPTIONAL: slots to save lambdas
    lambdaPlateau = NULL,
    lambdaCV = NULL,
    lambdaOracle = NULL,
    initialize = function(data, REPEAT_BOOTSTRAP = 2e2, inflate_lambda = 1) {
      self$data <- data
      self$REPEAT_BOOTSTRAP <- REPEAT_BOOTSTRAP
      self$inflate_lambda <- inflate_lambda
    },
    get_lambda = function() {
      as.numeric(names(self$dict_boot))
    },
    get_value = function() {
      unname(self$dict_boot)
    },
    plot_CI = function(Psi = NULL) {
      # plot CI v.s. log(lambda)
      lambdas <- self$get_lambda()
      CI_list <- self$get_value()
      CI_all <- lapply(CI_list, function(x) x$CI_all)

      df_ls <- list()
      count <- 1
      for (i in 1:length(CI_all)) {
        for (j in 1:length(CI_all[[i]])) {
          df_ls[[count]] <- data.frame(
            lambda = lambdas[i],
            CI_low = CI_all[[i]][[j]][1],
            CI_upp = CI_all[[i]][[j]][2],
            kindCI = names(CI_all[[i]][j])
          )
          count <- count + 1
        }
      }
      df2 <- do.call(rbind, df_ls)
      # not plot centered version
      df2 <- df2[grep("ctr", df2$kindCI, invert = TRUE), ]
      df2$logLambda <- log10(df2$lambda)
      df2$center <- (df2$CI_low + df2$CI_upp) / 2

      library(ggplot2)
      p <- ggplot(df2, aes(
        x = logLambda,
        y = center,
        group = interaction(logLambda, kindCI),
        color = kindCI
      )) +
        geom_point(position = position_dodge(0.5)) +
        geom_errorbar(
          aes(ymin = CI_low, ymax = CI_upp),
          width = .5,
          position = position_dodge(0.5)
        ) +
        ylab("")
      if (!is.null(Psi)) p <- p + geom_hline(yintercept = Psi, linetype = 2)
      return(p)
    },
    get_lambda_df = function() {
      lambdas <- self$get_lambda()
      CI_list <- self$get_value()
      width_all <- lapply(CI_list, function(x) x$width_all)

      df_ls <- list()
      count <- 1
      for (i in 1:length(width_all)) {
        for (j in 1:length(width_all[[i]])) {
          df_ls[[count]] <- data.frame(
            lambda = lambdas[i],
            width = width_all[[i]][[j]],
            kindCI = names(width_all[[i]][j])
          )
          count <- count + 1
        }
      }
      df2 <- do.call(rbind, df_ls)
      # not use centered version
      self$df_lambda_width <- df2[grep("ctr", df2$kindCI, invert = TRUE), ]
      return(self$df_lambda_width)
    },
    plot_width = function(type_CI = NULL) {
      # plot CIwidth v.s. log(lambda)
      if (is.null(type_CI)) type_CI <- unique(self$df_lambda_width$kindCI)
      library(ggplot2)
      p1 <- ggplot(
        self$df_lambda_width[self$df_lambda_width$kindCI %in% type_CI, ],
        aes(x = lambda, y = width, group = kindCI, color = kindCI)
      ) +
        geom_point() +
        geom_line() +
        ylab("width of interval") +
        scale_x_log10() +
        theme_bw()
      p2 <- ggplot(
        self$df_lambda_width[self$df_lambda_width$kindCI %in% type_CI, ],
        aes(x = l1, y = width, group = kindCI, color = kindCI)
      ) +
        geom_point() +
        geom_line() +
        ylab("width of interval") +
        theme_bw()
      # return(p1)
      return(list(p1 = p1, p2 = p2))
    },
    select_lambda_pleateau_wald = function(df_lambda_width) {
      # grab pleateau (when wald plateaus)
      if (!is.null(self$lambdaCV)) {
        # don't check for plateau for lambda >= lambda_CV
        df_lambda_width <- df_lambda_width[df_lambda_width$lambda <= self$lambdaCV, ]
      }
      df_ls <- list()
      count <- 1
      for (kind in "wald") {
        tempDf <- df_lambda_width[df_lambda_width$kindCI == kind, ]
        tempDf <- tempDf[order(tempDf$lambda), ]

        platOut <- grabPlateau$new(x = log10(tempDf$lambda), y = tempDf$width)
        coordOut <- platOut$find_plateau(find_plateau_start = TRUE)
        # coordOut <- platOut$find_plateau(find_plateau_start = FALSE)
        coordOut$kindCI <- kind

        df_ls[[count]] <- coordOut
        count <- count + 1
      }
      allPlateaus <- do.call(rbind, df_ls)
      return(10^allPlateaus$x)
    },
    get_Psi_df = function() {
      lambdas <- self$get_lambda()
      CI_list <- self$get_value()
      Psi_all <- sapply(CI_list, function(x) x$Psi)
      se_Psi <- sapply(CI_list, function(x) x$bootOut$pointTMLE$se_Psi)
      df_Psi <- data.frame(lambda = lambdas, Psi = Psi_all, se_Psi = se_Psi)
      return(df_Psi)
    },
    select_lambda_psi_grad = function(df_Psi) {
      grad_Psi <- diff(df_Psi$Psi)
      grad_se <- diff(df_Psi$se_Psi)

      monotonicfy_the_grad <- function(grad) {
        sign_first_IsPositive <- tail(grad, 1) > 0
        if (sign_first_IsPositive) {
          grad[grad < 0] <- 0
        } else {
          grad[grad > 0] <- 0
        }
        return(grad)
      }
      grad_Psi <- monotonicfy_the_grad(grad_Psi)
      grad_se <- monotonicfy_the_grad(grad_se)

      alpha <- 1.96
      objective1 <- abs(grad_Psi - alpha * grad_se)
      objective2 <- abs(grad_Psi + alpha * grad_se)
      # -1 and +1 is to make the indexing start from 0,
      # so that mod operator will minus the length when needed
      lambda_index <- (which.min(c(objective1, objective2)) - 1) %% length(objective1) + 1
      lambda_balanceMSE <- df_Psi$lambda[lambda_index]
      return(lambda_balanceMSE)
    }
  )
)

#' @export
avgDensity_LambdaGrid <- R6Class("avgDensity_LambdaGrid",
  inherit = LambdaGrid,
  public = list(
    bin_width = NULL,
    epsilon_step = NULL,
    initialize = function(bin_width, epsilon_step, ...) {
      super$initialize(...)
      self$bin_width <- bin_width
      self$epsilon_step <- epsilon_step
      return(self)
    },
    add_lambda = function(lambda_grid = NULL, to_parallel = FALSE, ...) {
      library(foreach)
      `%mydo%` <- ifelse(to_parallel, `%dopar%`, `%do%`)
      new_ls <- foreach(
        lambda = lambda_grid,
        .combine = c,
        .inorder = TRUE,
        .errorhandling = "pass",
        .export = c("self"),
        .verbose = T
      ) %mydo% {
        boot_here <- comprehensiveBootstrap$new(
          parameter = avgDensityBootstrap,
          x = self$data$x,
          bin_width = self$bin_width,
          lambda_grid = lambda,
          epsilon_step = self$epsilon_step
        )
        boot_here$compute_wald_width()
        return(boot_here)
      }
      # named list. the name is the lambda used for fitting
      names(new_ls) <- formatC(lambda_grid, format = "e", digits = 5)
      self$dict_boot <- c(self$dict_boot, new_ls)
    }
  )
)

#' @export
ATE_LambdaGrid <- R6Class("ATE_LambdaGrid",
  inherit = LambdaGrid,
  public = list(
    data = NULL,
    initialize = function(...) {
      super$initialize(...)
      return(self)
    },
    add_lambda = function(lambda_grid = NULL, to_parallel = FALSE, ...) {
      library(foreach)
      `%mydo%` <- ifelse(to_parallel, `%dopar%`, `%do%`)
      new_ls <- foreach(
        lambda1 = lambda_grid,
        .combine = c,
        .inorder = TRUE,
        .errorhandling = "pass",
        .export = c("self"),
        .verbose = T
      ) %mydo% {
        boot_here <- comprehensiveBootstrap$new(
          parameter = ateBootstrap,
          data = self$data,
          lambda1 = lambda1,
          ...
        )
        boot_here$compute_wald_width()
        return(boot_here)
      }
      # named list. the name is the lambda used for fitting
      names(new_ls) <- formatC(lambda_grid, format = "e", digits = 5)
      self$dict_boot <- c(self$dict_boot, new_ls)
    }
  )
)

#' @export
blipVar_contY_LambdaGrid <- R6Class("blipVar_contY_LambdaGrid",
  inherit = LambdaGrid,
  public = list(
    data = NULL,
    initialize = function(...) {
      super$initialize(...)
      return(self)
    },
    add_lambda = function(lambda_grid = NULL, to_parallel = FALSE, ...) {
      library(foreach)
      `%mydo%` <- ifelse(to_parallel, `%dopar%`, `%do%`)
      new_ls <- foreach(
        lambda1 = lambda_grid,
        .combine = c,
        .inorder = TRUE,
        .errorhandling = "pass",
        .export = c("self"),
        .verbose = T
      ) %mydo% {
        boot_here <- comprehensiveBootstrap$new(
          parameter = blipVarianceBootstrap_contY,
          data = self$data,
          lambda1 = lambda1,
          ...
        )
        boot_here$compute_wald_width()
        return(boot_here)
      }
      # named list. the name is the lambda used for fitting
      names(new_ls) <- formatC(lambda_grid, format = "e", digits = 5)
      self$dict_boot <- c(self$dict_boot, new_ls)
    }
  )
)

#' @export
grabPlateau <- R6Class("grabPlateau",
  # grab a plateau of a function y = f(x)
  public = list(
    x = NULL, # x needs to be sorted
    y = NULL,
    initialize = function(x, y) {
      self$x <- x
      self$y <- y
    },
    find_plateau = function(find_plateau_start = TRUE) {
      # argmin(sec diff)
      dx <- diff(self$x)
      dx_shift <- c(dx[2:length(dx)], NA)
      denom <- dx * dx_shift
      denom <- denom[!is.na(denom)]
      sec_diff <- diff(diff(self$y)) / (denom)
      if (find_plateau_start) {
        idx <- which.max(sec_diff)
      } else {
        idx <- which.min(sec_diff)
      }
      # fix when there is no plateau
      if (length(idx) == 0) {
        idx <- 1
      } else if (idx == 0) {
        idx <- 1
      }
      return(data.frame(y = self$y[idx], x = self$x[idx]))
    }
  )
)
