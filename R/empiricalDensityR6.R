# library(R6)
#' @export
empiricalDensity <- R6Class("empiricalDensity",
  # helper for avgDens TMLE; hold univariate density and normalization
  public = list(
    p_density = NULL,
    x = NULL,
    initialize = function(p_density, x) {
      # class for convenient normalizing density estimate to 1
      self$p_density <- p_density
      self$x <- x
    },
    normalize = function() {
      dummy_df <- data.frame(id = 1:length(self$x), x = self$x, p_density = self$p_density)
      dummy_df <- dummy_df[order(dummy_df$x),]
      dx <- c(0,diff(dummy_df$x))
      dummy_df$p_density <- dummy_df$p_density/sum(dummy_df$p_density * dx)
      dummy_df <- dummy_df[order(dummy_df$id),]
      self$p_density <- dummy_df$p_density
      return(self)
    },
    display = function(p_truth = NULL, ...) {
      # plot the density; can give another function to overlay two densities
      plot(self$p_density ~ self$x, ...)
      # overlay a true p function on the plot
      if (!is.null(p_truth)) curve(p_truth, from = -10, to = 10, n = 1e3, add = TRUE)
    }
    )
)