#' Fit a 1-d density using HAL regression
#'
#' perform binning first then run a HAL classifier
#' @import hal9001
#' @export
densityHAL <- R6Class("densityHAL",
  # for fixed lambda; fit a hal fit on longitudinal dataframe; outout density estimate
  # input: longiData on whole data, for split schema
  # output: tuned HAL
  public = list(
    longiData = NULL,
    x = NULL,
    hal_fit = NULL, # hal9001 object; fit using optimal lambda; for output
    initialize = function(x, longiData) {
      self$longiData <- longiData
      self$x <- x
    },
    fit = function(lambda = 2e-5, ...) {
      # use `longiData` to transform `x`; fit a single lambda HAL on the data
      df_compressed <- self$longiData$generate_df_compress(x = self$x)
      self$hal_fit <- hal9001::fit_hal(
        X = df_compressed[, "box"],
        Y = df_compressed$Y,
        weights = df_compressed$Freq,
        family = "binomial",
        lambda = lambda,
        fit_type = "glmnet",
        use_min = TRUE, # useless
        return_lasso = TRUE,
        return_x_basis = FALSE,
        cv_select = FALSE,
        yolo = FALSE,
        ...
      )
    },
    predict = function(new_x = NULL) {
      # predict density p_hat on `new_x`
      if (length(new_x) > 1e4) return(self$predict_long(new_x = new_x))
      return(predict(self$hal_fit, new_data = new_x))
    },
    predict_long = function(new_x = NULL) {
      # make prediction faster on larger x
      message("use long prediction routine...")
      x_list <- split(new_x, ceiling(seq_along(new_x) / 20))
      out <- lapply(x_list, function(x) self$predict(new_x = x))
      return(do.call("c", out))
    },
    eval_misclass_loss = function(new_x = NULL) {
      # evaluate misclassification error loss on `new_x`
      df_valid <- self$longiData$generate_df_compress(x = new_x)
      yhat <- self$predict(new_x = df_valid$box)
      yhat_class <- (yhat > .5) + 0L
      return(sum(abs(yhat_class - df_valid$Y) * df_valid$Freq) / sum(df_valid$Freq))
    },
    eval_crossentropy_loss = function(new_x = NULL) {
      # evaluate cross-entropy loss on `new_x`
      df_valid <- self$longiData$generate_df_compress(x = new_x)
      yhat <- self$predict(new_x = df_valid$box)
      not_weighted <- cross_entropy(y = df_valid$Y, yhat = yhat)
      return(sum(not_weighted * df_valid$Freq) / sum(df_valid$Freq))
    }
  )
)
#' Fit a 1-d density using HAL regression; automatic tuning of L1 penalty
#'
#' perform binning first then run a HAL classifier
#' @export
cvDensityHAL <- R6Class("cvDensityHAL",
  # cross validate a grid of `densityHAL` with a grid of lambda
  public = list(
    longiData = NULL,
    x = NULL,
    folds = NULL,
    # results
    results = NULL,
    lambda.min = NULL,
    initialize = function(x, longiData) {
      self$longiData <- longiData
      self$x <- x
    },
    assign_fold = function(n_fold = 3) {
      # create `origami` fold; subsample `x`
      self$folds <- origami::make_folds(n = length(self$x), V = n_fold)
    },
    cv = function(lambda = 2e-5, verbose = FALSE, ...) {
      # fix one lambda; evaluate the validation loss on folds; take the mean loss
      cv_once <- function(fold, data, longiData, lambda, verbose = FALSE, ...) {
        if (verbose) message(paste("fitting lambda =", lambda))
        # define training and validation sets based on input object of class "folds"
        train_data <- origami::training(data)
        valid_data <- origami::validation(data)

        HALfit <- densityHAL$new(x = train_data, longiData = longiData)
        HALfit$fit(lambda = lambda, ...)
        return(list(loss = HALfit$eval_crossentropy_loss(new_x = valid_data)))
      }
      cv_results <- origami::cross_validate(
        cv_fun = cv_once,
        folds = self$folds,

        data = self$x,
        longiData = self$longiData,
        lambda = lambda,
        verbose = verbose,
        ...
      )
      return(mean(cv_results$loss))
    },
    cv_lambda_grid = function(lambda_grid = NULL, lambda_min_ratio = NULL, ...) {
      # repeat `cv` with a grid of lambda; store the error for each lambda;
      # pick the lambda minimizer of validation loss
      # c(1e-6,2e-5)
      # OPTIONAL: glmnet to get lambda_grid
      get_lambda_grid_from_cv <- is.null(lambda_grid)
      if ((!get_lambda_grid_from_cv) & (!is.null(lambda_min_ratio))) {
        stop("do not set `lambda_min_ratio` if you do not CV!")
      }
      if (get_lambda_grid_from_cv) {
        df_compressed <- self$longiData$generate_df_compress(x = self$x)
        # WILSON: has bug. There can be all Y=0 in the training sample
        hal_for_lambda <- hal9001::fit_hal(
          X = df_compressed[, "box"],
          Y = df_compressed$Y,
          weights = df_compressed$Freq,
          family = "binomial",
          fit_type = "glmnet",
          use_min = TRUE,
          return_lasso = TRUE,
          return_x_basis = FALSE,
          cv_select = TRUE,
          yolo = FALSE
        )
        lambda_grid <- hal_for_lambda$hal_lasso$lambda
        if (!is.null(lambda_min_ratio)) {
          # manually increase the range of lambda grid
          lambda_grid <- create_lambda_grid_by_ratio(
            lambda_grid, lambda_min_ratio
          )
        }
      }

      list_df <- list()
      b <- 1
      for (i in lambda_grid) {
        loss <- self$cv(lambda = i, ...)
        list_df[[b]] <- c(i, loss)
        b <- b + 1
      }
      results <- as.data.frame(do.call(rbind, list_df))
      colnames(results) <- c("lambda", "loss")
      self$results <- results
      self$lambda.min <- results$lambda[which.min(results$loss)]
    },
    compute_model_full_data = function(lambda) {
      # re-fit the best lambda model on the entire dataset;
      # output the `densityHAL` object
      HALfit_out <- densityHAL$new(x = self$x, longiData = self$longiData)
      HALfit_out$fit(lambda = lambda)
      return(HALfit_out)
    }
  )
)
#' cross-entropy loss
#'
#' @keywords internal
cross_entropy <- function(y, yhat) {
  -log(yhat) * as.numeric(y == 1) - log(1 - yhat) * as.numeric(y == 0)
}

#' @keywords internal
expit <- function(x) exp(x) / (1 + exp(x))


#' @keywords internal
create_lambda_grid_by_ratio <- function(lambda_grid, lambda_min_ratio) {
  lambda_max <- max(lambda_grid)
  lambda_min <- lambda_max * lambda_min_ratio
  log_lambda_range <- log(c(lambda_min, lambda_max))
  lambda_grid_new <- exp(seq(
    log_lambda_range[2],
    log_lambda_range[1],
    length.out = 1e2
  ))
  return(lambda_grid_new)
}
