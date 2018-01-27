require(R6)
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
        # Y <- ((i <= x) & (x < i+bin_width)) + 0L
        Y <- ((i-.5*self$bin_width <= x) & (x < i+.5*self$bin_width)) + 0L
        all_df[[b]] <- data.frame(id = 1:length(x), Y = Y, box = i)
        b <- b+1
      }
      all_df <- do.call("rbind", all_df)
      # self$longiDF <- all_df[order( all_df[,1], all_df[,3] ),2:3]
      return(all_df[order( all_df[,1], all_df[,3] ),2:3])
    },
    generate_df_compress = function(x = NULL) {
      if (is.null(x)) x <- self$x
      all_df <- list()
      b <- 1
      for (i in self$grids) {
        Y <- ((i-.5*self$bin_width <= x) & (x < i+.5*self$bin_width)) + 0L
        df_counts <- as.data.frame(table(Y))
        all_df[[b]] <- data.frame(box = i, df_counts)
        b <- b+1
      }
      all_df <- do.call("rbind", all_df)
      return(all_df)
    }
  )
)