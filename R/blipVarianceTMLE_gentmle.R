library(R6)
#' @export
blipVarianceTMLE_gentmle <- R6Class("blipVarianceTMLE_gentmle",
  public = list(
    data = NULL,
    # initial fit
    Q_fit = NULL,
    Q_AW = NULL,
    Q_1W = NULL,
    Q_0W = NULL,
    g_fit = NULL,
    # g_AW = NULL,
    g_1W = NULL,

    Psi = NULL,
    gentmle_object = NULL,
    se_Psi = NULL,
    CI = NULL,
    verbose = FALSE,
    initialize = function(data, epsilon_step = NULL, verbose = NULL) {
      self$data <- data
      if(class(data$W) != 'data.frame') message('W not data.frame')
      if (!is.null(verbose)) self$verbose <- verbose
    },
    initial_fit = function(){
      library(hal9001)
      # Q fit
      self$Q_fit <- hal9001::fit_hal(X = data.frame(self$data$A, self$data$W),
                                     Y = self$data$Y,
                                     n_folds = 3,
                                     fit_type = 'glmnet',
                                     family = 'binomial',
                                     yolo = FALSE)
      self$Q_AW <- plogis(hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(self$data$A, self$data$W)))
      self$Q_1W <- plogis(hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(1, self$data$W)))
      self$Q_0W <- plogis(hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(0, self$data$W)))
      # g fit
      self$g_fit <- hal9001::fit_hal(X = self$data$W,
                                     Y = self$data$A,
                                     n_folds = 3,
                                     fit_type = 'glmnet',
                                     family = 'binomial',
                                     yolo = FALSE)
      # self$g_AW <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
      self$g_1W <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
    },
    target = function() {
      library(gentmle2)
      initdata <- data.frame(A = self$data$A,
                             Y = self$data$Y,
                             gk = self$g_1W,
                             Qk = self$Q_AW,
                             Q1k = self$Q_1W,
                             Q0k = self$Q_0W)
      self$gentmle_object <- gentmle(initdata = initdata,
                                     params = list(param_sigmaATE),
                                     approach = "full",
                                     # max_iter = 1e5,
                                     submodel = submodel_logit)
      self$Psi <- self$gentmle_object$tmleests

      self$se_Psi <- sd(self$gentmle_object$Dstar)/sqrt(length(self$data$A))
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
    }
))

#' @export
blipVarianceTMLE_gentmle_contY <- R6Class("blipVarianceTMLE_gentmle_contY",
  inherit = blipVarianceTMLE_gentmle,
  public = list(
    # continuous Y
    scale_Y = NULL,
    Y_rescale = NULL,

    scale_Q = NULL,
    Q_AW_rescale = NULL,
    Q_1W_rescale = NULL,
    Q_0W_rescale = NULL,
    scaleY = function(){
      # scale Y to (0,1)
      self$scale_Y <- scaleX$new(X = self$data$Y)
      self$Y_rescale <- self$scale_Y$scale01(newX = self$data$Y)
    },
    scaleBack_afterTMLE = function(){
      self$Psi <- self$Psi * self$scale_Y$rangeX^2 # variance is rescaled
      self$se_Psi <- self$se_Psi * self$scale_Y$rangeX^2 # EIC is scaled by the same amount
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
    },
    initial_fit = function(){
      # message('continuous Y')
      library(hal9001)
      # Q fit
      self$Q_fit <- hal9001::fit_hal(X = data.frame(self$data$A, self$data$W),
                                     Y = self$data$Y,
                                     n_folds = 3,
                                     fit_type = 'glmnet',
                                     family = 'gaussian',
                                     yolo = FALSE)
      self$Q_AW <- hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(self$data$A, self$data$W))
      self$Q_1W <- hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(1, self$data$W))
      self$Q_0W <- hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(0, self$data$W))
      # g fit
      self$g_fit <- hal9001::fit_hal(X = self$data$W,
                                     Y = self$data$A,
                                     n_folds = 3,
                                     fit_type = 'glmnet',
                                     family = 'binomial',
                                     yolo = FALSE)
      # self$g_AW <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
      self$g_1W <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
      # scale Q to (0,1)
      self$scale_Q <- scaleX$new(X = c(self$Q_1W, self$Q_0W))
      self$Q_AW_rescale <- self$scale_Q$scale01(newX = self$Q_AW)
      self$Q_1W_rescale <- self$scale_Q$scale01(newX = self$Q_1W)
      self$Q_0W_rescale <- self$scale_Q$scale01(newX = self$Q_0W)
    },
    target = function() {
      # message('continuous Y')
      library(gentmle2)
      initdata <- data.frame(A = self$data$A,
                             Y = self$Y_rescale,
                             gk = self$g_1W,
                             Qk = self$Q_AW_rescale,
                             Q1k = self$Q_1W_rescale,
                             Q0k = self$Q_0W_rescale)
      self$gentmle_object <- gentmle(initdata = initdata,
                                     params = list(param_sigmaATE),
                                     approach = "full",
                                     # max_iter = 1e5,
                                     submodel = submodel_logit)
      self$Psi <- self$gentmle_object$tmleests
      self$se_Psi <- sd(self$gentmle_object$Dstar)/sqrt(length(self$data$A))
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
    }
))
