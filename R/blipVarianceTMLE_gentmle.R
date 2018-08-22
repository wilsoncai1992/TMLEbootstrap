library(R6)
library(hal9001)
library(gentmle2)



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
      # HAL initial fit for blip variance TMLE (iterative); targeting is done in gentmle2
      self$data <- data
      if(class(data$W) != 'data.frame') message('W not data.frame')
      if (!is.null(verbose)) self$verbose <- verbose
    },
    initial_fit = function(lambda1 = NULL, lambda2 = NULL){
      # hal9001 to fit binary Q(Y|A,W), and g(A|W); save the fit object
      # Q fit
      if(is.null(lambda1)){ # use CV
        self$Q_fit <- hal9001::fit_hal(X = data.frame(self$data$A, self$data$W),
                                       Y = self$data$Y,
                                       n_folds = 3,
                                       family = 'binomial',
                                       fit_type = 'glmnet',
                                       yolo = FALSE)
      }else{
        self$Q_fit <- hal9001::fit_hal_single_lambda(X = data.frame(self$data$A, self$data$W),
                                       Y = self$data$Y,
                                       family = 'binomial',
                                       fit_type = 'glmnet',
                                       lambda = lambda1,
                                       yolo = FALSE)
      }
      self$Q_AW <- plogis(hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(self$data$A, self$data$W)))
      self$Q_1W <- plogis(hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(1, self$data$W)))
      self$Q_0W <- plogis(hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(0, self$data$W)))
      # g fit
      if(is.null(lambda2)){ # use CV
        self$g_fit <- hal9001::fit_hal(X = self$data$W,
                                       Y = self$data$A,
                                       n_folds = 3,
                                       fit_type = 'glmnet',
                                       family = 'binomial',
                                       yolo = FALSE)
      }else{
        self$g_fit <- hal9001::fit_hal(X = self$data$W,
                                       Y = self$data$A,
                                       fit_type = 'glmnet',
                                       family = 'binomial',
                                       lambda = lambda2,
                                       yolo = FALSE)
      }
      # self$g_AW <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
      self$g_1W <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
    },
    target = function() {
      # use the logistic submodel to target blip variance; output Psi, EIC, CI
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
    },
    compute_EIC = function(Y, A, Q1k, Q0k, Qk, gk, psi){
      # helper function to compute EIC values. shared among blipVarianceTMLE class
      HA = 2 * (Q1k - Q0k - mean(Q1k-Q0k)) * (A/gk - (1 - A)/(1 - gk))
      IC = HA*(Y-Qk)+(Q1k - Q0k - mean(Q1k-Q0k))^2 - psi
      return(IC)
    },
    inference_without_target = function(){
      # apply parameter mapping without doing targeting step
      self$Psi <- var(self$Q_1W - self$Q_0W)
      EIC <- self$compute_EIC(Y = self$data$Y,
                              A = self$data$A,
                              Q1k = self$Q_1W,
                              Q0k = self$Q_0W,
                              Qk = self$Q_AW,
                              gk = self$g_1W,
                              psi = self$Psi)
      self$se_Psi <- sqrt(var(EIC)/length(EIC))
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
      # sub-class of `blipVarianceTMLE_gentmle`;
      # replace initial fit with continuous hal9001;
      # rescales Y into (0,1) before TMLE;
      # scales the Psi, EIC, CI back to the original scale after the TMLE
      # scale Y to (0,1)
      self$scale_Y <- scaleX$new(X = self$data$Y)
      self$Y_rescale <- self$scale_Y$scale01(newX = self$data$Y)
    },
    scaleBack_afterTMLE = function(){
      self$Psi <- self$Psi * self$scale_Y$rangeX^2 # variance is rescaled
      self$se_Psi <- self$se_Psi * self$scale_Y$rangeX^2 # EIC is scaled by the same amount
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
    },
    initial_fit = function(lambda1 = NULL, lambda2 = NULL){
      # message('continuous Y')
      # Q fit
      if(is.null(lambda1)){ # use CV
        self$Q_fit <- hal9001::fit_hal(X = data.frame(self$data$A, self$data$W),
                                       Y = self$data$Y,
                                       n_folds = 3,
                                       fit_type = 'glmnet',
                                       family = 'gaussian',
                                       yolo = FALSE)
      }else{
        self$Q_fit <- hal9001::fit_hal_single_lambda(X = data.frame(self$data$A, self$data$W),
                                       Y = self$data$Y,
                                       fit_type = 'glmnet',
                                       family = 'gaussian',
                                       lambda = lambda1,
                                       yolo = FALSE)
      }
      self$Q_AW <- hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(self$data$A, self$data$W))
      self$Q_1W <- hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(1, self$data$W))
      self$Q_0W <- hal9001:::predict.hal9001(self$Q_fit, new_data = data.frame(0, self$data$W))
      # g fit
      if(is.null(lambda2)){ # use CV
        self$g_fit <- hal9001::fit_hal(X = self$data$W,
                                       Y = self$data$A,
                                       n_folds = 3,
                                       fit_type = 'glmnet',
                                       family = 'binomial',
                                       yolo = FALSE)
      }else{
        self$g_fit <- hal9001::fit_hal_single_lambda(X = self$data$W,
                                       Y = self$data$A,
                                       fit_type = 'glmnet',
                                       family = 'binomial',
                                       lambda = lambda2,
                                       yolo = FALSE)
      }
      self$g_1W <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$W)))
      # scale Q to (0,1)
      self$scale_Q <- scaleX$new(X = c(self$Q_1W, self$Q_0W))
      self$Q_AW_rescale <- self$scale_Q$scale01(newX = self$Q_AW)
      self$Q_1W_rescale <- self$scale_Q$scale01(newX = self$Q_1W)
      self$Q_0W_rescale <- self$scale_Q$scale01(newX = self$Q_0W)
    },
    target = function() {
      # message('continuous Y')
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
      # self$se_Psi <- self$gentmle_object$ED2/sqrt(length(self$data$A)) # this is from jeremy
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
    },
    inference_without_target = function(){
      # apply parameter mapping without doing targeting step
      self$Psi <- var(self$Q_1W_rescale - self$Q_0W_rescale)
      EIC <- self$compute_EIC(Y = self$Y_rescale,
                              A = self$data$A,
                              Q1k = self$Q_1W_rescale,
                              Q0k = self$Q_0W_rescale,
                              Qk = self$Q_AW_rescale,
                              gk = self$g_1W,
                              psi = self$Psi)
      self$se_Psi <- sqrt(var(EIC)/length(EIC))
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
    }
))
