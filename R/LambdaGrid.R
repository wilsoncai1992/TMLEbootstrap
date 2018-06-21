#' @export
LambdaGrid <- R6Class("LambdaGrid",
  # lambda plateau method
  public = list(
    REPEAT_BOOTSTRAP = NULL,
    inflate_lambda = NULL,

    dict_boot = list(), # named list of comprehensiveBootstrap
    df_lambda_width = NULL,

    # OPTIONAL: slots to save lambdas
    lambdaPlateau = NULL,
    lambdaCV = NULL,
    lambdaOracle = NULL,
    initialize = function() {
    },
    get_lambda = function() {
      as.numeric(names(self$dict_boot))
    },
    get_value = function() {
      unname(self$dict_boot)
    },
    plot_CI = function(Psi = NULL){
      # plot CI v.s. log(lambda)
      lambdas <- self$get_lambda()
      CI_list <- self$get_value()
      CI_all <- lapply(CI_list, function(x) x$CI_all)

      df_ls <- list()
      count <- 1
      for (i in 1:length(CI_all)) {
        for (j in 1:length(CI_all[[i]])) {
          df_ls[[count]] <- data.frame(lambda = lambdas[i], CI_low = CI_all[[i]][[j]][1], CI_upp = CI_all[[i]][[j]][2], kindCI = names(CI_all[[i]][j]))
          count <- count + 1
        }
      }
      df2 <- do.call(rbind, df_ls)
      df2 <- df2[grep('ctr', df2$kindCI, invert = TRUE),] # not plot centered version
      df2$logLambda <- log10(df2$lambda)
      df2$center <- (df2$CI_low + df2$CI_upp)/2

      library(ggplot2)
      p <- ggplot(df2, aes(x = logLambda, y = center, group=interaction(logLambda, kindCI), color=kindCI)) +
        geom_point(position = position_dodge(0.5)) +
        geom_errorbar(aes(ymin = CI_low, ymax = CI_upp), width = .5, position = position_dodge(0.5)) +
        ylab('')
      if(!is.null(Psi)) p <- p + geom_hline(yintercept = Psi, linetype = 2)
      return(p)
    },
    get_lambda_df = function(){
      lambdas <- self$get_lambda()
      CI_list <- self$get_value()
      width_all <- lapply(CI_list, function(x) x$width_all)

      df_ls <- list()
      count <- 1
      for (i in 1:length(width_all)) {
        for (j in 1:length(width_all[[i]])) {
          df_ls[[count]] <- data.frame(lambda = lambdas[i], width = width_all[[i]][[j]], kindCI = names(width_all[[i]][j]))
          count <- count + 1
        }
      }
      df2 <- do.call(rbind, df_ls)
      self$df_lambda_width <- df2[grep('ctr', df2$kindCI, invert = TRUE),] # not use centered version
    },
    plot_width = function(){
      # plot CIwidth v.s. log(lambda)
      library(ggplot2)
      p <- ggplot(self$df_lambda_width, aes(x = log10(lambda), y = width, group=kindCI, color=kindCI)) +
        geom_point() +
        geom_line() +
        ylab('width of interval')
      return(p)
    },
    select_lambda_pleateau_wald = function(){
      # grab pleateau (when wald plateaus)
      df_ls <- list()
      count <- 1
      for (kind in 'wald') {
        tempDf <- self$df_lambda_width[self$df_lambda_width$kindCI == kind,]
        tempDf <- tempDf[order(tempDf$lambda),]

        platOut <- grabPlateau$new(x = log10(tempDf$lambda), y = tempDf$width)
        # coordOut <- platOut$plateau_1()
        coordOut <- platOut$plateau_2()
        coordOut$kindCI <- kind

        df_ls[[count]] <- coordOut
        count <- count + 1
      }
      allPlateaus <- do.call(rbind, df_ls)
      return(10^allPlateaus$x)
    }
    # select_lambda_pleateau_reg = function(){
    #   # grab pleateau (when reg plateaus)
    #   df_ls <- list()
    #   count <- 1
    #   for (kind in 'reg_scale') {
    #     tempDf <- self$df_lambda_width[self$df_lambda_width$kindCI == kind,]
    #     tempDf <- tempDf[order(tempDf$lambda),]

    #     platOut <- grabPlateau$new(x = log10(tempDf$lambda), y = tempDf$width)
        # coordOut <- platOut$plateau_1()
    #     coordOut <- platOut$plateau_2()
    #     coordOut$kindCI <- kind

    #     df_ls[[count]] <- coordOut
    #     count <- count + 1
    #   }
    #   allPlateaus <- do.call(rbind, df_ls)
    #   return(10^allPlateaus$x)
    # },
    # select_lambda_pleateau_secscalepen = function(){
    #   # grab pleateau (when secOrd + pen + scale plateaus)
    #   df_ls <- list()
    #   count <- 1
    #   for (kind in 'secOrd_scale_pen') {
    #     tempDf <- self$df_lambda_width[self$df_lambda_width$kindCI == kind,]
    #     tempDf <- tempDf[order(tempDf$lambda),]

    #     platOut <- grabPlateau$new(x = log10(tempDf$lambda), y = tempDf$width)
        # coordOut <- platOut$plateau_1()
    #     coordOut <- platOut$plateau_2()
    #     coordOut$kindCI <- kind

    #     df_ls[[count]] <- coordOut
    #     count <- count + 1
    #   }
    #   allPlateaus <- do.call(rbind, df_ls)
    #   return(10^allPlateaus$x)
    # },
    # select_lambda_pleateau_all = function(){
    #   # grab pleateau (when average width plateaus)
    #   df_ls <- list()
    #   count <- 1
    #   for (kind in unique(self$df_lambda_width$kindCI)) {
    #     tempDf <- self$df_lambda_width[self$df_lambda_width$kindCI == kind,]
    #     tempDf <- tempDf[order(tempDf$lambda),]

    #     platOut <- grabPlateau$new(x = log10(tempDf$lambda), y = tempDf$width)
        # coordOut <- platOut$plateau_1()
    #     coordOut <- platOut$plateau_2()
    #     # browser()
    #     coordOut$kindCI <- kind
    #     df_ls[[count]] <- coordOut
    #     count <- count + 1
    #   }
    #   allPlateaus <- do.call(rbind, df_ls)
    #   return(10^median(allPlateaus$x))
    # }
  )
)

#' @export
avgDensity_LambdaGrid <- R6Class("avgDensity_LambdaGrid",
  inherit = LambdaGrid,
  public = list(
    data = NULL,
    bin_width = NULL,
    epsilon_step = NULL,
    initialize = function(data,
                          bin_width,
                          epsilon_step,
                          REPEAT_BOOTSTRAP = 2e2,
                          inflate_lambda = 1) {
      self$data <- data
      self$bin_width <- bin_width
      self$epsilon_step <- epsilon_step
      self$REPEAT_BOOTSTRAP <- REPEAT_BOOTSTRAP
      self$inflate_lambda <- inflate_lambda
    },
    add_lambda = function(lambda_grid = NULL) {
      new_ls <- list()
      for (lambda in lambda_grid) {
        boot_here <- comprehensiveBootstrap$new(parameter = avgDensityBootstrap,
                                                x = self$data$x,
                                                bin_width = self$bin_width,
                                                lambda_grid = lambda,
                                                epsilon_step = self$epsilon_step)
        boot_here$bootstrap(REPEAT_BOOTSTRAP = self$REPEAT_BOOTSTRAP,
                            inflate_lambda = self$inflate_lambda)
        boot_here$all_CI()
        boot_here$compute_width()
        new_ls <- c(new_ls, boot_here)
        message(paste(lambda, 'is added'))
      }
      names(new_ls) <- lambda_grid # named list. the name is the lambda used for fitting
      self$dict_boot <- c(self$dict_boot, new_ls)
    }
  )
)

#' @export
ATE_LambdaGrid <- R6Class("ATE_LambdaGrid",
  inherit = LambdaGrid,
  public = list(
    data = NULL,
    initialize = function(data,
                          REPEAT_BOOTSTRAP = 2e2,
                          inflate_lambda = 1) {
      self$data <- data
      self$REPEAT_BOOTSTRAP <- REPEAT_BOOTSTRAP
      self$inflate_lambda <- inflate_lambda
    },
    add_lambda = function(lambda_grid = NULL) {
      new_ls <- list()
      for (lambda1 in lambda_grid) {
        boot_here <- comprehensiveBootstrap$new(parameter = ateBootstrap,
                                                data = self$data,
                                                lambda1 = lambda1)
        boot_here$bootstrap(REPEAT_BOOTSTRAP = self$REPEAT_BOOTSTRAP)
        boot_here$all_CI()
        boot_here$compute_width()
        new_ls <- c(new_ls, boot_here)
        message(paste(lambda1, 'is added'))
      }
      names(new_ls) <- lambda_grid # named list. the name is the lambda1 used for fitting
      self$dict_boot <- c(self$dict_boot, new_ls)
    }
  )
)

#' @export
blipVar_contY_LambdaGrid <- R6Class("blipVar_contY_LambdaGrid",
  inherit = LambdaGrid,
  public = list(
    data = NULL,
    initialize = function(data,
                          REPEAT_BOOTSTRAP = 2e2,
                          inflate_lambda = 1) {
      self$data <- data
      self$REPEAT_BOOTSTRAP <- REPEAT_BOOTSTRAP
      self$inflate_lambda <- inflate_lambda
    },
    add_lambda = function(lambda_grid = NULL) {
      new_ls <- list()
      for (lambda1 in lambda_grid) {
        boot_here <- comprehensiveBootstrap$new(parameter = blipVarianceBootstrap_contY,
                                                data = self$data,
                                                lambda1 = lambda1)
        boot_here$bootstrap(REPEAT_BOOTSTRAP = self$REPEAT_BOOTSTRAP)
        boot_here$all_CI()
        boot_here$compute_width()
        new_ls <- c(new_ls, boot_here)
        message(paste(lambda1, 'is added'))
      }
      names(new_ls) <- lambda_grid # named list. the name is the lambda1 used for fitting
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
    plateau_1 = function() {
      # immediate after largest cliff
      firstDiff <- diff(self$y)/diff(self$x)
      # plot(firstDiff)
      idx <- which.min(firstDiff) - 1
      if (idx == 0) idx <- 1 # fix when there is no plateau
      return(data.frame(y = self$y[idx], x = self$x[idx]))
    },
    plateau_2 = function() {
      # argmin(sec diff)
      dx <- diff(self$x)
      dx_shift <- c(dx[2:length(dx)], NA)
      denom <- dx * dx_shift
      denom <- denom[!is.na(denom)]
      secDiff <- diff(diff(self$y))/(denom)
      # plot(secDiff)
      idx <- which.min(secDiff)
      if (idx == 0) idx <- 1 # fix when there is no plateau
      return(data.frame(y = self$y[idx], x = self$x[idx]))
    }
  )
)