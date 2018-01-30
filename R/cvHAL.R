require(R6)
#' @export
densityHAL <- R6Class("densityHAL",
  public = list(
    longiData = NULL,
    x = NULL,
    hal_fit = NULL,
    initialize = function(x, longiData) {
      # input: longiData on whole data, for split schema
      # output: tuned HAL
      self$longiData <- longiData
      self$x <- x
    },
    fit = function(lambda = 2e-5){
      df_compressed <- self$longiData$generate_df_compress(x = self$x)
      self$hal_fit <- hal9001::fit_hal_single_lambda(X = df_compressed[,'box'],
        Y = df_compressed$Y,
        weights = df_compressed$Freq,
        family = "binomial",
        lambda = lambda,
        fit_type = 'glmnet',
        use_min = TRUE, #useless
        yolo = FALSE)
      # self$hal_fit$lambda_star
    },
    predict = function(new_x = NULL){
      return(rje::expit(predict(self$hal_fit, new_data = new_x)))
    },
    eval_misclass_loss = function(new_x = NULL){
      df_valid <- self$longiData$generate_df_compress(x = new_x)
      yhat <- self$predict(new_x = df_valid$box)
      yhat_class <- (yhat > .5) + 0L
      return(sum(abs(yhat_class - df_valid$Y) * df_valid$Freq) / sum(df_valid$Freq))
    }
  )
)

#' @export
cv_densityHAL <- R6Class("cv_densityHAL",
  public = list(
    longiData = NULL,
    x = NULL,
    folds = NULL,
    initialize = function(x, longiData) {
      self$longiData <- longiData
      self$x <- x
    },
    assign_fold = function(n_fold = 3){
      sample_per_fold <- ceiling(length(self$x) / n_fold)
      fold_assign <- sample(head(rep(1:n_fold, each = sample_per_fold), n = length(self$x)))
      x_folds <- list()
      for (i in 1:n_fold) {
        x_folds[[i]] <- self$x[fold_assign == i]
      }
      self$folds <- x_folds
    },
    cv = function(n_fold = 3){
      folds <- origami::make_folds(n = length(self$x), V = n_fold)
      cv_results <- origami::cross_validate(
        cv_fun = cv_once, folds = folds,
        data = x
      )

    }
  )
)

cv_once <- function(fold, data, longiData, lambda){
  # define training and validation sets based on input object of class "folds"
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  HALfit <- densityHAL$new(x = train_data, longiData = longiData)
  HALfit$fit(lambda = lambda)
  HALfit$eval_misclass_loss(new_x = valid_data)
  return()
}