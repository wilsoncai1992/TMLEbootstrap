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
    },
    eval_crossentropy_loss = function(new_x = NULL){
      df_valid <- self$longiData$generate_df_compress(x = new_x)
      yhat <- self$predict(new_x = df_valid$box)
      not_weighted <- cross_entropy(y = df_valid$Y, yhat = yhat)
      return(sum(not_weighted * df_valid$Freq) / sum(df_valid$Freq))
    }
    # density_intial <- empiricalDensity$new(p_density = yhat, x = x)
    # p_hat <- density_intial$normalize()
  )
)

#' @export
cv_densityHAL <- R6Class("cv_densityHAL",
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
    assign_fold = function(n_fold = 3){
      self$folds <- origami::make_folds(n = length(self$x), V = n_fold)
    },
    cv = function(lambda = 2e-5){
      cv_once <- function(fold, data, longiData, lambda){
        # define training and validation sets based on input object of class "folds"
        train_data <- origami::training(data)
        valid_data <- origami::validation(data)

        HALfit <- densityHAL$new(x = train_data, longiData = longiData)
        HALfit$fit(lambda = lambda)
        # return(HALfit$eval_misclass_loss(new_x = valid_data))
        return(list(loss = HALfit$eval_crossentropy_loss(new_x = valid_data)))
      }
      cv_results <- origami::cross_validate(
        cv_fun = cv_once, folds = self$folds,
        data = self$x, longiData = self$longiData, lambda = lambda
      )
      return(mean(cv_results$loss))
    },
    cv_lambda_grid = function(lambda_grid = c(1e-6,2e-5)){
      list_df <- list()
      b <- 1
      for (i in lambda_grid) {
        loss <- self$cv(lambda = i)
        list_df[[b]] <- c(i, loss)
        b <- b+1
      }
      results <- as.data.frame(do.call(rbind, list_df))
      colnames(results) <- c('lambda', 'loss')
      self$results <- results
      self$lambda.min <- results$lambda[which.min(results$loss)]
    },
    compute_best_model = function(){
      HALfit_out <- densityHAL$new(x = self$x, longiData = self$longiData)
      HALfit_out$fit(lambda = self$lambda.min)
      return(HALfit_out)
    }
  )
)

#' @export
cross_entropy <- function(y, yhat) -log(yhat) * as.numeric(y == 1) -log(1 - yhat) * as.numeric(y == 0)
