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

require(R6)
#' @export
longiData <- R6Class("longiData",
  # helper for avgDens TMLE. convert univariate series to longitudinal format
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

#' @export
longiData_resample <- R6Class("longiData_resample",
  # perform bootstrap on longitudinal format; respecting weights
  public = list(
    df_compressed = NULL,
    n_row = NULL,
    n_sample = NULL,
    initialize = function(df_compressed) {
      self$df_compressed <- df_compressed
      self$n_row <- nrow(self$df_compressed)
      self$n_sample <- sum(self$df_compressed$Freq)
    },
    create_folds = function(n_fold = 3) {
      artificial_ind <- rep(1:self$n_row, self$df_compressed$Freq)
      sample_per_fold <- ceiling(self$n_sample / n_fold)
      fold_assign <- sample(
        head(rep(1:n_fold, each = sample_per_fold), n = self$n_sample)
      )

      df_compressed_folds <- list()
      for (i in 1:n_fold) {
        table_here <- as.data.frame(table(artificial_ind[fold_assign == i]))
        table_here$Var1 <- as.numeric(as.character(table_here$Var1))
        temp_df <- self$df_compressed[table_here$Var1, ]
        temp_df$Freq <- table_here$Freq
        df_compressed_folds[[i]] <- longiData_resample$new(df_compressed = temp_df)
      }
      return(df_compressed_folds)
    },
    bootstrap_with_replacement = function(n = self$n_sample) {
      samp_idx <- sample(
        seq_len(self$n_row), n,
        prob = self$df_compressed$Freq, replace = TRUE
      )
      samp_idx_tbl <- as.data.frame(table(samp_idx))
      samp_idx_tbl$samp_idx <- as.numeric(as.character(samp_idx_tbl$samp_idx))
      new_df_compressed <- self$df_compressed[samp_idx_tbl$samp_idx, ]
      new_df_compressed$Freq <- samp_idx_tbl$Freq
      return(new_df_compressed)
    }
  )
)

#' @export
sum_longiData_resample <- function(list) {
  # combine a list of longiData
  # input: list of longiData_resample
  # output: a new longiData_resample
  library(dplyr)
  list_df <- lapply(list, function(x) x$df_compressed)
  df_long <- do.call(rbind, list_df)
  out <- data.frame(
    df_long %>%
      group_by(box, Y) %>%
      summarise(Freq = sum(Freq)) %>%
      ungroup()
  )
  return(longiData_resample$new(df_compressed = out))
}
