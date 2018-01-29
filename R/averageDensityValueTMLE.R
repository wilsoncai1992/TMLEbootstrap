library(R6)

#' @export
avgDensityTMLE <- R6Class("avgDensityTMLE",
  public = list(
    x = NULL,
    p_hat = NULL,
    Psi = NULL,
    EIC = NULL,
    epsilon_step = 1e-3,
    tol = NULL,
    CI = NULL,
    verbose = FALSE,
    max_iter = 1e2,
    longDataOut = NULL,
    HAL_tuned = NULL,
    initialize = function(x, epsilon_step = NULL, verbose = NULL) {
      self$x <- x
      self$tol <- 1/length(x)
      if (!is.null(epsilon_step)) self$epsilon_step <- epsilon_step
      if (!is.null(verbose)) self$verbose <- verbose
    },
    fit_density = function(bin_width = .2) {
      browser()
      library(SuperLearner)
      library(hal9001)
      self$longDataOut <- longiData$new(x = self$x, bin_width = bin_width)
      longDFOut <- self$longDataOut$generate_df_compress()

      verbose <- FALSE
      # tune HAL for density
      HAL_tuned <- hal9001::fit_hal(X = longDFOut[,'box'],
        Y = longDFOut$Y,
        weights = longDFOut$Freq,
        use_min = TRUE,
        yolo = FALSE,
        fit_type = 'glmnet',
        family = "binomial",
        n_folds = 3)
      yhat <- rje::expit(predict(HAL_tuned, new_data = self$longDataOut$x))
      density_intial <- empiricalDensity$new(p_density = yhat, x = self$x)
      self$p_hat <- density_intial$normalize()

      # SL_fit <- SuperLearner(Y = longDFOut$Y, X = longDFOut[,'box',F], newX = data.frame(box = self$longDataOut$x),
      # family = 'binomial',
      # SL.library = "SL.hal9001",
      # cvControl = list(V = 3),
      # verbose = verbose)
      # HAL_tuned <- SL_fit$fitLibrary$SL.hal9001_All$object
      self$HAL_tuned <- hal9001::squash_hal_fit(HAL_tuned)

      # density_intial <- empiricalDensity$new(p_density = SL_fit$SL.predict, x = self$x)
      # self$p_hat <- density_intial$normalize()

      # foo2 <- function(x) {(.5*dnorm(x, mean = 2) + .5*dnorm(x, mean = -2))}
      # density_intial$display(foo2)
    },
    calc_Psi = function(){
      # self$Psi <- mean(self$p_hat$p_density)

      dummy_df <- data.frame(id = 1:length(self$x), x = self$x, p_density = self$p_hat$p_density)
      dummy_df <- dummy_df[order(dummy_df$x),]
      dx <- c(0,diff(dummy_df$x))
      self$Psi <- sum(dummy_df$p_density^2 * dx)
    },
    calc_EIC = function() {
      self$EIC <- 2 * (self$p_hat$p_density - self$Psi)
    },
    updateOnce = function() {
      if (mean(self$EIC) < 0) self$epsilon_step <- -self$epsilon_step
      self$p_hat$p_density <- self$p_hat$p_density * exp(self$epsilon_step * self$EIC)
      self$p_hat$normalize()
    },
    onestepTarget = function(verbose = FALSE) {
      n_iter <- 0
      meanEIC_prev <- abs(mean(self$EIC))
      while(abs(mean(self$EIC)) >= self$tol){
      # while(abs(mean(self$EIC)) >= 1e-20){
      # while(TRUE){
        meanEIC_prev <- abs(mean(self$EIC))
        self$calc_Psi()
        self$calc_EIC()
        # self$p_hat$display()
        self$updateOnce()
        if (self$verbose | verbose) print(c(mean(self$EIC), self$Psi))
        n_iter <- n_iter + 1
        if (abs(mean(self$EIC)) > meanEIC_prev){
          self$epsilon_step <- -self$epsilon_step
          # message('not stable!')
        }
        if (n_iter >= self$max_iter){
          break()
          message('max iteration number reached!')
        }
      }
    },
    inference = function(){
      sd_EIC <- sd(self$EIC)
      upper <- self$Psi + 1.96/sqrt(length(self$EIC))*sd_EIC
      lower <- self$Psi - 1.96/sqrt(length(self$EIC))*sd_EIC
      self$CI <- c(lower, upper)
    }
  )
)