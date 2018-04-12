#' @export
avgDensity_LambdaGrid <- R6Class("avgDensity_LambdaGrid",
  public = list(
    data = NULL,
    bin_width = NULL,
    epsilon_step = NULL,
    REPEAT_BOOTSTRAP = NULL,
    inflate_lambda = NULL,

    # dict_boot = dict::dict(), # dictionary of comprehensiveBootstrap
    dict_boot = list(), # named list of comprehensiveBootstrap
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
        # self$dict_boot[[lambda]] <- boot_here
        message(paste(lambda, 'is added'))
      }
      names(new_ls) <- lambda_grid # named list. the name is the lambda used for fitting
      self$dict_boot <- c(self$dict_boot, new_ls)
    },
    get_lambda = function() {
      as.numeric(names(self$dict_boot))
    },
    get_value = function() {
      unname(self$dict_boot)
    },
    plot_CI = function(Psi = NULL){
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
    plot_width = function(){
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
      df2 <- df2[grep('ctr', df2$kindCI, invert = TRUE),] # not plot centered version

      library(ggplot2)
      p <- ggplot(df2, aes(x = log10(lambda), y = width, group=kindCI, color=kindCI)) +
        geom_point() +
        geom_line() +
        ylab('width of interval')
      return(p)
    }
  )
)