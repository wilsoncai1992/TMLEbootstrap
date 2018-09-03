library(R6)
#' @export
scaleX <- R6Class("scaleX",
  # perform standardization based on (x - min)/(max - min)
  public = list(
    X = NULL,
    minX = NULL,
    maxX = NULL,
    rangeX = NULL,
    initialize = function(X) {
      self$X <- X
      self$minX <- min(X)
      self$maxX <- max(X)
      self$rangeX <- self$maxX - self$minX
    },
    scale01 = function(newX = NULL) {
      return((newX - self$minX) / self$rangeX)
    },
    scaleOriginal = function(newX = NULL) {
      return(newX * self$rangeX + self$minX)
    }
  )
)
