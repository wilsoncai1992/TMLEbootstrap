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
    g_AW = NULL,
    g_1W = NULL,

    Psi = NULL,
    gentmle_object = NULL,
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
      self$g_AW <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
      self$g_1W <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
    },
    target = function() {
      library(gentmle2)
      initdata <- data.frame(A = self$data$A,
                             Y = self$data$Y,
                             gk = self$g_AW,
                             Qk = self$Q_AW,
                             Q1k = self$Q_1W,
                             Q0k = self$Q_0W)
      self$gentmle_object <- gentmle(initdata = initdata,
                                     params = list(param_sigmaATE),
                                     approach = "full",
                                     # max_iter = 1e5,
                                     submodel = submodel_logit)
      self$Psi <- self$gentmle_object$tmleests

      se_Psi <- sd(self$gentmle_object$Dstar)/sqrt(length(self$data$A))
      self$CI <- self$Psi + c(-1.96, 1.96) * se_Psi
    }
))

#' @export
blipVarianceTMLE_gentmle_contY <- R6Class("blipVarianceTMLE_gentmle_contY",
  inherit = blipVarianceTMLE_gentmle,
  public = list(
    # initial_fit = function(){
    #   library(hal9001)
    #   # Q fit
    #   self$Q_fit <- hal9001::fit_hal(X = data.frame(self$data$A, self$data$W),
    #                                  Y = self$data$Y,
    #                                  n_folds = 3,
    #                                  fit_type = 'glmnet',
    #                                  family = 'binomial',
    #                                  yolo = FALSE)
    #   self$Q_AW <- plogis(hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(self$data$A, self$data$W)))
    #   self$Q_1W <- plogis(hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(1, self$data$W)))
    #   self$Q_0W <- plogis(hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(0, self$data$W)))
    #   # g fit
    #   self$g_fit <- hal9001::fit_hal(X = self$data$W,
    #                                  Y = self$data$A,
    #                                  n_folds = 3,
    #                                  fit_type = 'glmnet',
    #                                  family = 'binomial',
    #                                  yolo = FALSE)
    #   self$g_AW <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
    #   self$g_1W <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
    # },
    # target = function() {
    #   library(gentmle2)
    #   initdata <- data.frame(A = self$data$A,
    #                          Y = self$data$Y,
    #                          gk = self$g_AW,
    #                          Qk = self$Q_AW,
    #                          Q1k = self$Q_1W,
    #                          Q0k = self$Q_0W)
    #   self$gentmle_object <- gentmle(initdata = initdata,
    #                                  params = list(param_sigmaATE),
    #                                  approach = "full",
    #                                  # max_iter = 1e5,
    #                                  submodel = submodel_logit)
    #   self$Psi <- self$gentmle_object$tmleests
    # 
    #   se_Psi <- sd(self$gentmle_object$Dstar)/sqrt(length(self$data$A))
    #   self$CI <- self$Psi + c(-1.96, 1.96) * se_Psi
    # }
))
