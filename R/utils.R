#' Store a 1-dimensional density function
#'
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
      dummy_df <- data.frame(
        id = 1:length(self$x), x = self$x, p_density = self$p_density
      )
      dummy_df <- dummy_df[order(dummy_df$x), ]
      dx <- c(0, diff(dummy_df$x))
      dummy_df$p_density <- dummy_df$p_density / sum(dummy_df$p_density * dx)
      dummy_df <- dummy_df[order(dummy_df$id), ]
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

#' Convert univariate series to longitudinal format
#'
#' helper for avgDens TMLE.
#' @export
longiData <- R6Class("longiData",
  public = list(
    x = NULL,
    bin_width = NULL,
    # longiDF = NULL,
    grids = NULL,
    initialize = function(x, bin_width = .5) {
      self$x <- x
      self$bin_width <- bin_width
      # x <- floor(x)
      bounds <- range(x)
      self$grids <- seq(from = bounds[1], to = bounds[2], by = bin_width)
    },
    generate_df = function(x = NULL) {
      if (is.null(x)) x <- self$x
      all_df <- list()
      b <- 1
      for (i in self$grids) {
        Y <- ((i - .5 * self$bin_width <= x) & (x < i + .5 * self$bin_width)) + 0L
        all_df[[b]] <- data.frame(id = 1:length(x), Y = Y, box = i)
        b <- b + 1
      }
      all_df <- do.call("rbind", all_df)
      all_df$Y <- as.numeric(as.character(all_df$Y)) # turn factor to numeric
      # self$longiDF <- all_df[order( all_df[,1], all_df[,3] ),2:3]
      return(all_df[order(all_df[, 1], all_df[, 3]), 2:3])
    },
    generate_df_compress = function(x = NULL) {
      if (is.null(x)) x <- self$x
      all_df <- list()
      b <- 1
      for (i in self$grids) {
        Y <- ((i - .5 * self$bin_width <= x) & (x < i + .5 * self$bin_width)) + 0L
        df_counts <- as.data.frame(table(Y))
        all_df[[b]] <- data.frame(box = i, df_counts)
        b <- b + 1
      }
      all_df <- do.call("rbind", all_df)
      all_df$Y <- as.numeric(as.character(all_df$Y)) # turn factor to numeric
      return(all_df)
    }
  )
)

#' Perform min/max standardization based on (x - min)/(max - min)
#'
#' helper for blipVarianceTMLEContinuousY
#' @keywords internal
scaleX <- R6Class("scaleX",
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
