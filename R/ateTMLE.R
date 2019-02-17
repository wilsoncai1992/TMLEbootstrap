library(R6)
library(tmle)
library(hal9001)
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
    EIC = NULL,
    # verbose = FALSE,
    lambda1 = NULL,
    lambda2 = NULL,
    initialize = function(data) {
      self$data <- data
      if (class(data$W) != "data.frame") message("W not data.frame")
    },
    initial_fit = function(
      lambda1 = NULL, lambda2 = NULL, M1 = NULL, M2 = NULL, ...
    ) {
      use_penalized_mode <- any(c(!is.null(lambda1), !is.null(lambda2)))
      use_constrained_mode <- any(c(!is.null(M1), !is.null(M2)))
      if (use_penalized_mode & use_constrained_mode) {
        stop("cannot do two modes!")
      }
      if (use_penalized_mode) {
        self$initial_fit_pen_likeli(lambda1 = lambda1, lambda2 = lambda2, ...)
      }
      if (use_constrained_mode) {
        self$initial_fit_constrained_form(M1 = M1, M2 = M2, ...)
      }
      if (!use_penalized_mode & !use_constrained_mode) {
        # when user have NULL for everything: default to CV
        self$initial_fit_pen_likeli(NULL, NULL, ...)
      }

      # get Q_1W, Q_0W
      self$Q_1W <- stats::predict(object = self$Q_fit, new_data = data.frame(1, self$data$W))
      self$Q_0W <- stats::predict(object = self$Q_fit, new_data = data.frame(0, self$data$W))
      # get g1_W
      self$g1_W <- stats::predict(object = self$g_fit, new_data = data.frame(self$data$W))
    },
    initial_fit_pen_likeli = function(
      lambda1 = NULL, lambda2 = NULL, lambda_min_ratio = NULL, n_folds = 3, ...
    ) {
      # lambda1 for Q fit
      # lambda2 for g fit
      self$lambda1 <- lambda1
      self$lambda2 <- lambda2
      # hal9001 to fit binary Q(Y|A,W), and g(A|W); save the fit object
      # Q fit
      if (is.null(lambda1)) {
        # use CV
        self$Q_fit <- hal9001::fit_hal(
          X = data.frame(self$data$A, self$data$W),
          Y = self$data$Y,
          family = "gaussian",
          fit_type = "glmnet",
          n_folds = 3,
          use_min = TRUE,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE,
          ...
        )
        # WILSON hack the lambda_min_ratio
        if (!is.null(lambda_min_ratio)) {
          lambda_max = max(self$Q_fit$hal_lasso$lambda)
          lambda_min = lambda_max * lambda_min_ratio
          log_lambda_range <- log(c(lambda_min, lambda_max))
          lambda_grid_new <- exp(seq(log_lambda_range[2], log_lambda_range[1], length.out = 1e2))
          self$Q_fit <- hal9001::fit_hal(
            X = data.frame(self$data$A, self$data$W),
            Y = self$data$Y,
            family = "gaussian",
            fit_type = "glmnet",
            n_folds = n_folds,
            use_min = TRUE,
            lambda = lambda_grid_new,
            return_lasso = TRUE,
            return_x_basis = FALSE,
            yolo = FALSE,
            ...
          )
        }
      } else if (lambda1 >= 0) {
        # use manual lambda1
        self$Q_fit <- hal9001::fit_hal(
          X = data.frame(self$data$A, self$data$W),
          Y = self$data$Y,
          family = "gaussian",
          lambda = lambda1,
          fit_type = "glmnet",
          use_min = TRUE, # useless
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE,
          ...
        )
      }
      # g fit
      if (is.null(lambda2)) {
        # use CV
        self$g_fit <- hal9001::fit_hal(
          X = data.frame(self$data$W),
          Y = self$data$A,
          family = "binomial",
          fit_type = "glmnet",
          n_folds = n_folds,
          use_min = TRUE,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE,
          ...
        )
      } else {
        # use manual lambda1
        self$g_fit <- hal9001::fit_hal(
          X = data.frame(self$data$W),
          Y = self$data$A,
          family = "binomial",
          lambda = lambda2,
          fit_type = "glmnet",
          use_min = TRUE, # useless
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE,
          ...
        )
      }
    },
    initial_fit_constrained_form = function(
      M1 = NULL, M2 = NULL, n_folds = 3, ...
    ) {
      # M1 for Q fit
      # M2 for g fit
      # self$M1 <- M1
      # self$M2 <- M2

      # hal9001 to fit binary Q(Y|A,W), and g(A|W); save the fit object
      # Q fit
      if (is.null(M1)) {
        self$Q_fit <- hal9001::fit_hal(
          X = data.frame(self$data$A, self$data$W),
          Y = self$data$Y,
          family = "gaussian",
          fit_type = "glmnet",
          n_folds = n_folds,
          use_min = TRUE,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE,
          ...
        )
      } else if (M1 >= 0) { # use manual M1
        self$Q_fit <- hal9001::fit_hal_constraint_form(
          X = data.frame(self$data$A, self$data$W),
          Y = self$data$Y,
          family = "gaussian",
          fit_type = "glmnet",
          yolo = FALSE,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          M = M1,
          ...
        )
      }
      # g fit
      if (is.null(M2)) { # use CV
        self$g_fit <- hal9001::fit_hal(
          X = data.frame(self$data$W),
          Y = self$data$A,
          family = "binomial",
          fit_type = "glmnet",
          n_folds = n_folds,
          use_min = TRUE,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          yolo = FALSE,
          ...
        )
      } else { # use manual M1
        self$g_fit <- hal9001::fit_hal_constraint_form(
          X = data.frame(self$data$W),
          Y = self$data$A,
          family = "binomial",
          fit_type = "glmnet",
          yolo = FALSE,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          M = M2,
          ...
        )
      }
    },
    plot_Q1W = function(foo = NULL) {
      # plot the Q(1,W) function (optional: against a foo function)
      plot(self$Q_1W ~ self$data$W[, 1], col = "blue")
      if (!is.null(foo)) {
        curve(expr = foo, from = -10, to = 10, add = TRUE, lty = 2, n = 1e3)
      }
      points(
        self$data$Y[self$data$A == 1] ~ self$data$W[self$data$A == 1, 1],
        col = "grey"
      )
    },
    target = function() {
      # perform iterative TMLE
      self$tmle_object <- tmle::tmle(
        Y = self$data$Y,
        A = self$data$A,
        W = as.matrix(self$data$W),
        Q = cbind(self$Q_0W, self$Q_1W),
        g1W = self$g1_W,
        family = "gaussian",
        fluctuation = "logistic",
        V = 3,
        verbose = FALSE
      )
      self$Psi <- self$tmle_object$estimates$ATE$psi
      self$se_Psi <- sqrt(self$tmle_object$estimates$ATE$var.psi)
      self$CI <- self$tmle_object$estimates$ATE$CI
      self$EIC <- self$tmle_object$estimates$IC$IC.ATE
    },
    compute_EIC = function(A, gk, Y, Qk, Q1k, Q0k, psi) {
      HA <- A / gk - (1 - A) / (1 - gk)
      EIC <- HA * (Y - Qk) + Q1k - Q0k - psi
      return(EIC)
    },
    inference_without_target = function(data, Q_fit = NULL, g_fit = NULL, to_return = TRUE) {
      # input: data from `simulate_data`: named list with W, A, Y
      if (class(Q_fit) == 'hal9001') {
        # if the Q is a hal9001 fit, do the prediction routine
        Q_1W <- stats::predict(object = Q_fit, new_data = data.frame(1, data$W))
        Q_0W <- stats::predict(object = Q_fit, new_data = data.frame(0, data$W))
      } else if (class(Q_fit) == 'function'){
        # if the Q is the true fit (in function format). evaluate the true Q function
        # Q_1W <- Q_fit(w = data$W$W, a = 1)
        # Q_0W <- Q_fit(w = data$W$W, a = 0)
        Q_1W <- Q_fit(w = data$W, a = 1) #WILSON HACK
        Q_0W <- Q_fit(w = data$W, a = 0) #WILSON HACK
      } else if (is.null(Q_fit)) {
        # if NULL, assume Q_1W has been manually input
        Q_1W <- self$Q_1W
        Q_0W <- self$Q_0W
      }

      if (class(g_fit) == 'hal9001') {
        # if the g is a hal9001 fit, do the prediction routine
        g1_W <- stats::predict(object = g_fit, new_data = data.frame(data$W))
      } else if (class(g_fit) == 'function') {
        # if the g is the true fit (in function format). evaluate the true Q function
        # g1_W <- g_fit(w = data$W$W)
        g1_W <- g_fit(w = data$W) #WILSON HACK
      } else if (is.null(g_fit)) {
        # if NULL, assume g1_W has been manually input
        g1_W <- self$g1_W
      }
      # apply parameter mapping without doing any targeting
      # A has to be 0/1 coding
      Psi <- mean(Q_1W - Q_0W)
      EIC <- self$compute_EIC(
        A = data$A,
        gk = g1_W,
        Y = data$Y,
        Qk = Q_1W * data$A + Q_0W * (1 - data$A),
        Q1k = Q_1W,
        Q0k = Q_0W,
        psi = Psi
      )
      se_Psi <- sqrt(var(EIC) / length(EIC))
      CI <- Psi + c(-1.96, 1.96) * se_Psi
      EIC <- EIC
      if (to_return){
        return(list(Psi = Psi, EIC = EIC))
      } else {
        self$Psi <- Psi
        self$se_Psi <- se_Psi
        self$CI <- CI
        self$EIC <- EIC
        if (!is.null(Q_fit)) self$Q_fit <- Q_fit
      }
    },
    compute_min_phi_ratio = function() {
      # return the ratio of 1 in the basis.
      # argmin over all columns where there are non-zero beta value (intercept excluded)
      Qbasis_list <- self$Q_fit$basis_list
      Qcopy_map <- self$Q_fit$copy_map
      X <- data.frame(self$data$A, self$data$W)
      if (length(Qbasis_list) > 0) {
        x_basis <- hal9001:::make_design_matrix(as.matrix(X), Qbasis_list)
        unique_columns <- as.numeric(names(Qcopy_map))
        # design matrix. each column correspond to Q_fit$coefs. don't have intercept column
        x_basis <- x_basis[, unique_columns]
        phi_ratio <- Matrix::colMeans(x_basis)

        length(self$Q_fit$coefs)
        beta_nonIntercept <- self$Q_fit$coefs[-1]
        beta_nonzero <- beta_nonIntercept != 0
        nonzeroBeta_phiRatio <- phi_ratio[beta_nonzero]
      } else {
        # there is no coef left
        nonzeroBeta_phiRatio <- numeric()
      }

      # return NULL if:
      # all beta are zero
      # Qbasis has zero length
      if (length(nonzeroBeta_phiRatio) != 0) {
        return(min(nonzeroBeta_phiRatio))
      } else {
        return(NULL)
      }
    }
  )
)
