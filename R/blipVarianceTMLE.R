library(R6)
library(hal9001)
library(gentmle2)

#' @export
#' @importFrom Matrix colMeans
blipVarianceTMLE <- R6Class("blipVarianceTMLE",
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
    EIC = NULL,
    verbose = FALSE,
    initialize = function(data, epsilon_step = NULL, verbose = NULL) {
      # HAL initial fit for blip variance TMLE (iterative);
      # targeting is done in gentmle2
      self$data <- data
      if (class(data$W) != "data.frame") message("W not data.frame")
      if (!is.null(verbose)) self$verbose <- verbose
    },
    initial_fit = function(lambda1 = NULL, lambda2 = NULL, M1 = NULL, M2 = NULL) {
      use_penalized_mode <- any(c(!is.null(lambda1), !is.null(lambda2)))
      use_constrained_mode <- any(c(!is.null(M1), !is.null(M2)))

      if (use_penalized_mode & use_constrained_mode) {
        stop("cannot do two modes!")
      }
      if (use_penalized_mode) {
        self$initial_fit_pen_likeli(lambda1 = lambda1, lambda2 = lambda2)
      }
      if (use_constrained_mode) {
        stop("not implemented")
      }
      if (!use_penalized_mode & !use_constrained_mode) {
        # when user have NULL for everything: default to CV
        self$initial_fit_pen_likeli(NULL, NULL)
      }

      # get Q_1W, Q_0W
      self$Q_AW <- hal9001:::predict.hal9001(
        self$Q_fit,
        new_data = data.frame(self$data$A, self$data$W)
      )
      self$Q_1W <- hal9001:::predict.hal9001(
        self$Q_fit,
        new_data = data.frame(1, self$data$W)
      )
      self$Q_0W <- hal9001:::predict.hal9001(
        self$Q_fit,
        new_data = data.frame(0, self$data$W)
      )
      # get g1_W
      self$g_1W <- hal9001:::predict.hal9001(
        self$g_fit,
        new_data = data.frame(self$data$W)
      )
    },
    initial_fit_pen_likeli = function(lambda1 = NULL, lambda2 = NULL) {
      # hal9001 to fit binary Q(Y|A,W), and g(A|W); save the fit object
      # Q fit
      if (is.null(lambda1)) {
        # use CV
        self$Q_fit <- hal9001::fit_hal(
          X = data.frame(self$data$A, self$data$W),
          Y = self$data$Y,
          n_folds = 3,
          family = "binomial",
          fit_type = "glmnet",
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE
        )
      } else {
        self$Q_fit <- hal9001::fit_hal(
          X = data.frame(self$data$A, self$data$W),
          Y = self$data$Y,
          family = "binomial",
          fit_type = "glmnet",
          lambda = lambda1,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE
        )
      }
      # g fit
      if (is.null(lambda2)) {
        # use CV
        self$g_fit <- hal9001::fit_hal(
          X = self$data$W,
          Y = self$data$A,
          n_folds = 3,
          fit_type = "glmnet",
          family = "binomial",
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE
        )
      } else {
        self$g_fit <- hal9001::fit_hal(
          X = self$data$W,
          Y = self$data$A,
          fit_type = "glmnet",
          family = "binomial",
          lambda = lambda2,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE
        )
      }
    },
    target = function() {
      # use the logistic submodel to target blip variance; output Psi, EIC, CI
      initdata <- data.frame(
        A = self$data$A,
        Y = self$data$Y,
        gk = self$g_1W,
        Qk = self$Q_AW,
        Q1k = self$Q_1W,
        Q0k = self$Q_0W
      )
      self$gentmle_object <- gentmle2::gentmle(
        initdata = initdata,
        params = list(gentmle2::param_sigmaATE),
        approach = "full",
        # max_iter = 1e5,
        submodel = gentmle2::submodel_logit
      )
      self$Psi <- self$gentmle_object$tmleests

      self$se_Psi <- sd(self$gentmle_object$Dstar) / sqrt(length(self$data$A))
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
    },
    compute_EIC = function(Y, A, Q1k, Q0k, Qk, gk, psi) {
      # helper function to compute EIC values. shared among blipVarianceTMLE class
      HA <- 2 * (Q1k - Q0k - mean(Q1k - Q0k)) * (A / gk - (1 - A) / (1 - gk))
      IC <- HA * (Y - Qk) + (Q1k - Q0k - mean(Q1k - Q0k))^2 - psi
      return(IC)
    },
    inference_without_target = function() {
      # apply parameter mapping without doing targeting step
      self$Psi <- var(self$Q_1W - self$Q_0W)
      EIC <- self$compute_EIC(
        Y = self$data$Y,
        A = self$data$A,
        Q1k = self$Q_1W,
        Q0k = self$Q_0W,
        Qk = self$Q_AW,
        gk = self$g_1W,
        psi = self$Psi
      )
      self$se_Psi <- sqrt(var(EIC) / length(EIC))
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
      self$EIC <- EIC
    },
    compute_min_phi_ratio = function() {
      # return the ratio of 1 in the basis.
      # argmin over all columns where there are non-zero beta value (intercept excluded)
      Qbasis_list <- self$Q_fit$basis_list
      Qcopy_map <- self$Q_fit$copy_map
      X <- data.frame(self$data$A, self$data$W)
      if (length(Qbasis_list) > 0) {
        x_basis <- hal9001::make_design_matrix(as.matrix(X), Qbasis_list)
        unique_columns <- as.numeric(names(Qcopy_map))
        # design matrix. each column correspond to Q_fit$coefs.
        # don't have intercept column
        x_basis <- x_basis[, unique_columns]
        phi_ratio <- Matrix::colMeans(x_basis)

        beta_non_intercept <- self$Q_fit$coefs[-1]
        beta_non_zero <- beta_non_intercept != 0
        nonzero_beta_phi_ratio <- phi_ratio[beta_non_zero]
      } else {
        # there is no coef left
        nonzero_beta_phi_ratio <- numeric()
      }

      # return NULL if: # all beta are zero OR # Qbasis has zero length
      if (length(nonzero_beta_phi_ratio) != 0) {
        return(min(nonzero_beta_phi_ratio))
      } else {
        return(NULL)
      }
    }
  )
)

#' @export
blipVarianceTMLEContinuousY <- R6Class("blipVarianceTMLEContinuousY",
  inherit = blipVarianceTMLE,
  public = list(
    # continuous Y
    scale_Y = NULL,
    Y_rescale = NULL,

    scale_Q = NULL,
    Q_AW_rescale = NULL,
    Q_1W_rescale = NULL,
    Q_0W_rescale = NULL,
    scale_outcome = function() {
      # sub-class of `blipVarianceTMLE`;
      # replace initial fit with continuous hal9001;
      # rescales Y into (0,1) before TMLE;
      # scales the Psi, EIC, CI back to the original scale after the TMLE
      # scale Y to (0,1)
      self$scale_Y <- scaleX$new(X = self$data$Y)
      self$Y_rescale <- self$scale_Y$scale01(newX = self$data$Y)
    },
    scale_back_after_tmle = function() {
      # variance is rescaled
      self$Psi <- self$Psi * self$scale_Y$rangeX^2
      # EIC is scaled by the same amount
      self$se_Psi <- self$se_Psi * self$scale_Y$rangeX^2
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
    },
    initial_fit = function(lambda1 = NULL, lambda2 = NULL, M1 = NULL, M2 = NULL) {
      use_penalized_mode <- any(c(!is.null(lambda1), !is.null(lambda2)))
      use_constrained_mode <- any(c(!is.null(M1), !is.null(M2)))

      if (use_penalized_mode & use_constrained_mode) {
        stop("cannot do two modes!")
      }
      if (use_penalized_mode) {
        self$initial_fit_pen_likeli(lambda1 = lambda1, lambda2 = lambda2)
      }
      if (use_constrained_mode) {
        stop("not implemented")
      }
      if (!use_penalized_mode & !use_constrained_mode) {
        # when user have NULL for everything: default to CV
        self$initial_fit_pen_likeli(NULL, NULL)
      }

      # get Q_1W, Q_0W
      self$Q_AW <- hal9001:::predict.hal9001(
        self$Q_fit,
        new_data = data.frame(self$data$A, self$data$W)
      )
      self$Q_1W <- hal9001:::predict.hal9001(
        self$Q_fit,
        new_data = data.frame(1, self$data$W)
      )
      self$Q_0W <- hal9001:::predict.hal9001(
        self$Q_fit,
        new_data = data.frame(0, self$data$W)
      )
      # get g1_W
      self$g_1W <- hal9001:::predict.hal9001(
        self$g_fit,
        new_data = data.frame(self$data$W)
      )
      # scale Q to (0,1)
      self$Q_AW_rescale <- self$scale_Y$scale01(newX = self$Q_AW)
      self$Q_1W_rescale <- self$scale_Y$scale01(newX = self$Q_1W)
      self$Q_0W_rescale <- self$scale_Y$scale01(newX = self$Q_0W)
      # self$scale_Q <- scaleX$new(X = c(self$Q_1W, self$Q_0W))
      # self$Q_AW_rescale <- self$scale_Q$scale01(newX = self$Q_AW)
      # self$Q_1W_rescale <- self$scale_Q$scale01(newX = self$Q_1W)
      # self$Q_0W_rescale <- self$scale_Q$scale01(newX = self$Q_0W)
    },
    initial_fit_pen_likeli = function(lambda1 = NULL, lambda2 = NULL) {
      # Q fit
      if (is.null(lambda1)) {
        # use CV
        self$Q_fit <- hal9001::fit_hal(
          X = data.frame(self$data$A, self$data$W),
          Y = self$data$Y,
          n_folds = 3,
          fit_type = "glmnet",
          family = "gaussian",
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE
        )
      } else {
        self$Q_fit <- hal9001::fit_hal(
          X = data.frame(self$data$A, self$data$W),
          Y = self$data$Y,
          fit_type = "glmnet",
          family = "gaussian",
          lambda = lambda1,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE
        )
      }
      # g fit
      if (is.null(lambda2)) {
        # use CV
        self$g_fit <- hal9001::fit_hal(
          X = self$data$W,
          Y = self$data$A,
          n_folds = 3,
          fit_type = "glmnet",
          family = "binomial",
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE
        )
      } else {
        self$g_fit <- hal9001::fit_hal(
          X = self$data$W,
          Y = self$data$A,
          fit_type = "glmnet",
          family = "binomial",
          lambda = lambda2,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE
        )
      }
    },
    target = function() {
      initdata <- data.frame(
        A = self$data$A,
        Y = self$Y_rescale,
        gk = self$g_1W,
        Qk = self$Q_AW_rescale,
        Q1k = self$Q_1W_rescale,
        Q0k = self$Q_0W_rescale
      )
      self$gentmle_object <- gentmle2::gentmle(
        initdata = initdata,
        params = list(gentmle2::param_sigmaATE),
        approach = "full",
        submodel = gentmle2::submodel_logit
      )
      self$Psi <- self$gentmle_object$tmleests
      self$se_Psi <- sd(self$gentmle_object$Dstar) / sqrt(length(self$data$A))
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
    },
    inference_without_target = function() {
      # apply parameter mapping without doing targeting step
      self$Psi <- var(self$Q_1W_rescale - self$Q_0W_rescale)
      EIC <- self$compute_EIC(
        Y = self$Y_rescale,
        A = self$data$A,
        Q1k = self$Q_1W_rescale,
        Q0k = self$Q_0W_rescale,
        Qk = self$Q_AW_rescale,
        gk = self$g_1W,
        psi = self$Psi
      )
      self$se_Psi <- sqrt(var(EIC) / length(EIC))
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
      self$EIC <- EIC
    }
  )
)
