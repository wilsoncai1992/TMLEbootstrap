library(R6)
#' @export
ateTMLE <- R6Class("ateTMLE",
  public = list(
    data = NULL,
    # initial fit
    Q_fit = NULL,
    Q_1W = NULL,
    Q_0W = NULL,
    g_fit = NULL,
    g1_W = NULL,

    tmle_object = NULL,
    Psi = NULL,
    se_Psi = NULL,
    CI = NULL,
    # verbose = FALSE,
    lambda1 = NULL,
    lambda2 = NULL,
    initialize = function(data) {
      self$data <- data
      if(class(data$W) != 'data.frame') message('W not data.frame')
    },
    initial_fit = function(lambda1 = NULL, lambda2 = NULL){
      self$lambda1 <- lambda1
      self$lambda2 <- lambda2
      # hal9001 to fit binary Q(Y|A,W), and g(A|W); save the fit object
      library(hal9001)
      # Q fit
      if(is.null(lambda1)){ # use CV
        self$Q_fit <- hal9001::fit_hal(X = data.frame(self$data$A, self$data$W),
                                       Y = self$data$Y,
                                       family = 'gaussian',
                                       fit_type = 'glmnet',
                                       n_folds = 3,
                                       use_min = TRUE,
                                       yolo = FALSE)
      }else{ # use manual lambda1
        self$Q_fit <- hal9001::fit_hal_single_lambda(X = data.frame(self$data$A, self$data$W),,
                                                      Y = self$data$Y,
                                                     family = 'gaussian',
                                                      lambda = lambda1,
                                                      fit_type = 'glmnet',
                                                      use_min = TRUE, #useless
                                                      yolo = FALSE)
      }
      # Q_HAL_tuned <- squash_hal_fit(Qfit)
      # g fit
      if(is.null(lambda2)){ # use CV
        self$g_fit <- hal9001::fit_hal(X = data.frame(self$data$W),
                                       Y = self$data$A,
                                       family = 'binomial',
                                       fit_type = 'glmnet',
                                       n_folds = 3,
                                       use_min = TRUE,
                                       yolo = FALSE)
      }else{ # use manual lambda1
        self$g_fit <- hal9001::fit_hal_single_lambda(X = data.frame(self$data$W),
                                               Y = self$data$A,
                                               family = 'binomial',
                                                lambda = lambda2,
                                                fit_type = 'glmnet',
                                                use_min = TRUE, #useless
                                                yolo = FALSE)
      }
      # g_HAL_tuned <- hal9001::squash_hal_fit(gfit)
      # get Q_1W, Q_0W
      self$Q_1W <- stats::predict(object = self$Q_fit, new_data = data.frame(1, self$data$W))
      self$Q_0W <- stats::predict(object = self$Q_fit, new_data = data.frame(0, self$data$W))
      # get g1_W
      self$g1_W <- plogis(stats::predict(object = self$g_fit, new_data = data.frame(self$data$W)))
    },
    plot_Q1W = function(foo = NULL){
      plot(self$Q_1W ~ self$data$W[,1], col = 'blue')
      if(!is.null(foo)) curve(expr = foo, from = -10, to = 10, add = TRUE, lty = 2, n=1e3)
      points(self$data$Y[self$data$A == 1] ~ self$data$W[self$data$A == 1,1], col = 'grey')

      # plot(Q_1W - Q_0W ~ self$data$W)
    },
    target = function() {
      library(tmle)
      self$tmle_object <- tmle(Y = self$data$Y, A = self$data$A, W = as.matrix(self$data$W),
                               Q = cbind(self$Q_0W, self$Q_1W),
                               g1W = self$g1_W,
                               family = 'gaussian',
                               fluctuation = 'linear',
                               V = 3, verbose = FALSE)
      self$Psi <- self$tmle_object$estimates$ATE$psi
      self$se_Psi <- sqrt(self$tmle_object$estimates$ATE$var.psi)
      self$CI <- self$tmle_object$estimates$ATE$CI
    }
))
