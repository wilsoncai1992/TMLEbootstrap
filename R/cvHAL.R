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
      # new_longiData <- self$longiData$generate_df(x = new_x)
      return(rje::expit(predict(self$hal_fit, new_data = new_x)))
    },
    eval_misclass_loss = function(new_x = NULL, new_y = NULL){
      yhat <- self$predict(new_x = new_x)
      yhat_class <- (yhat > .5) + 0L
      return(mean(abs(yhat - new_y)))
    }
  )
)

# #' @export
# cv_densityHAL <- R6Class("cv_densityHAL",
#   public = list(
#     longiData = NULL,
#     x = NULL,
#     initialize = function(x, longiData) {
#       self$longiData <- longiData
#       self$x <- x
#     }
#   )
# )
