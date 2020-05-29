#' Compute TMLE on the variance of CATE (binary Y)
#'
#' inherit the initial fit routine from ate
#' @export
#' @importFrom Matrix colMeans
blipVarianceTMLE <- R6Class("blipVarianceTMLE",
  inherit = ateTMLE,
  public = list(
    gentmle_object = NULL,
    initial_fit = function(family_y = "binomial", ...) {
      super$initial_fit(family_y = family_y, ...)
    },
    #' @description
    #' use the logistic submodel to target blip variance; output Psi, EIC, CI
    #'
    #' @return NULL
    target = function() {
      initdata <- data.frame(
        A = self$data$A,
        Y = self$data$Y,
        gk = self$g1W,
        Qk = self$QAW,
        Q1k = self$Q1W,
        Q0k = self$Q0W
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
    #' @description
    #' Compute EIC of blip variance
    #'
    #' @param Y continuous target variable
    #' @param A binary treatment
    #' @param Q1k E(Y|A=1,W)
    #' @param Q0k E(Y|A=0,W)
    #' @param Qk E(Y|A=a,W)
    #' @param gk E(A=a|W)
    #' @param psi var(E(Y|A=1,W) - E(Y|A=0,W))
    #'
    #' @return EIC vector
    compute_EIC = function(Y, A, Q1k, Q0k, Qk, gk, psi) {
      # helper function to compute EIC values. shared among blipVarianceTMLE class
      HA <- 2 * (Q1k - Q0k - mean(Q1k - Q0k)) * (A / gk - (1 - A) / (1 - gk))
      IC <- HA * (Y - Qk) + (Q1k - Q0k - mean(Q1k - Q0k))^2 - psi
      return(IC)
    },
    inference_without_target = function() {
      # apply parameter mapping without doing targeting step
      self$Psi <- var(self$Q1W - self$Q0W)
      EIC <- self$compute_EIC(
        Y = self$data$Y,
        A = self$data$A,
        Q1k = self$Q1W,
        Q0k = self$Q0W,
        Qk = self$QAW,
        gk = self$g1W,
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

#' Compute TMLE on the variance of CATE (continuous Y)
#'
#' inherit the initial fit routine from ateTMLE
#' replace initial fit with continuous hal9001;
#' rescales Y into (0,1) before TMLE;
#' scales the Psi, EIC, CI back to the original scale after the TMLE
#' scale Y to (0,1)
#' @export
blipVarianceTMLEContinuousY <- R6Class("blipVarianceTMLEContinuousY",
  inherit = blipVarianceTMLE,
  public = list(
    # continuous Y
    scale_Y = NULL,
    Y_rescale = NULL,

    scale_Q = NULL,
    QAW_rescale = NULL,
    Q1W_rescale = NULL,
    Q0W_rescale = NULL,
    scale_outcome = function() {
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
    initial_fit = function(...) {
      super$initial_fit(family_y = "gaussian", ...)
      # scale Q to (0,1)
      self$QAW_rescale <- self$scale_Y$scale01(newX = self$QAW)
      self$Q1W_rescale <- self$scale_Y$scale01(newX = self$Q1W)
      self$Q0W_rescale <- self$scale_Y$scale01(newX = self$Q0W)
    },
    target = function() {
      initdata <- data.frame(
        A = self$data$A,
        Y = self$Y_rescale,
        gk = self$g1W,
        Qk = self$QAW_rescale,
        Q1k = self$Q1W_rescale,
        Q0k = self$Q0W_rescale
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
      self$Psi <- var(self$Q1W_rescale - self$Q0W_rescale)
      EIC <- self$compute_EIC(
        Y = self$Y_rescale,
        A = self$data$A,
        Q1k = self$Q1W_rescale,
        Q0k = self$Q0W_rescale,
        Qk = self$QAW_rescale,
        gk = self$g1W,
        psi = self$Psi
      )
      self$se_Psi <- sqrt(var(EIC) / length(EIC))
      self$CI <- self$Psi + c(-1.96, 1.96) * self$se_Psi
      self$EIC <- EIC
    }
  )
)
