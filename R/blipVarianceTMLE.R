library(R6)
#' @export
blipVarianceTMLE <- R6Class("blipVarianceTMLE",
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
    # blip
    blip = NULL,
    Psi = NULL,
    # targeting
    H_AW = NULL,
    H_1W = NULL,
    H_0W = NULL,
    EIC = NULL,
    epsilon_n = NULL,

    H_AW_ATE = NULL,
    H_1W_ATE = NULL,
    H_0W_ATE = NULL,

    tol = NULL,
    CI = NULL,
    verbose = FALSE,
    max_iter = 1e2,
    initialize = function(data, epsilon_step = NULL, verbose = NULL) {
      self$data <- data
      if(class(data$W) != 'data.frame') message('W not data.frame')
      self$tol <- 1/nrow(data)
      if (!is.null(epsilon_step)) self$epsilon_step <- epsilon_step
      if (!is.null(verbose)) self$verbose <- verbose
    },
    initial_fit = function(){
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
      self$g_AW <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(self$data$A, self$data$W)))
      self$g_1W <- plogis(hal9001:::predict.hal9001(self$g_fit, new_data = data.frame(1, self$data$W)))
    },
    compute_blip = function(){
      self$blip <- self$Q_1W - self$Q_0W
    },
    compute_Psi = function(){
      self$Psi <- var(self$blip)
    },
    compute_EIC_H = function(){
      # input: blip; Q; g; psi; data
      self$H_AW <- 2 * (self$blip - mean(self$blip)) * (2*self$data$A - 1) / self$g_AW
      self$H_1W <- 2 * (self$blip - mean(self$blip)) * (2 - 1) / self$g_1W
      self$H_0W <- 2 * (self$blip - mean(self$blip)) * (0 - 1) / (1-self$g_1W) #g_0W
      self$EIC <- self$H_AW * (self$data$Y - self$Q_AW) + (self$blip - mean(self$blip))^2 - self$Psi
      # clever covariate for ATE
      self$H_AW_ATE <- (2*self$data$A - 1) / self$g_AW
      self$H_1W_ATE <- (2 - 1) / self$g_1W
      self$H_0W_ATE <- (0 - 1) / (1-self$g_1W) #g_0W
    },
    target_once = function(){
      # fit epsilon; linear fluctuation model
      fluctuation_fit <- glm.fit(x = self$H_AW, y = self$data$Y, offset = self$Q_AW, family = gaussian(), intercept = FALSE)
      self$epsilon_n <- fluctuation_fit$coefficients
      # update Q_n using TMLE
      # self$Q_1W <- self$Q_1W + self$H_1W * self$epsilon_n
      # self$Q_0W <- self$Q_0W + self$H_0W * self$epsilon_n
      # WRONG
      self$Q_1W <- self$Q_1W - self$H_1W * self$epsilon_n
      self$Q_0W <- self$Q_0W - self$H_0W * self$epsilon_n
    },
    compute_stopping = function(){
      return(mean(self$EIC))
    }
))


