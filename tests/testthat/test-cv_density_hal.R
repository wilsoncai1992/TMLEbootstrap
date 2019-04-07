context("cv_density_hal should estimate the true density reasonably")

simulate_data <- function(n_sim, n_mode) {
  # simulate 1d gaussian mixture density
  modes <- seq(from = -4, to = 4, length.out = n_mode)
  sigma <- 10 / n_mode / 6
  # sample object
  x <- rnorm(n = n_sim, mean = sample(x = modes, size = n_sim, replace = TRUE), sd = sigma)
  # function object
  dmixture <- function(x, modes, sigma) {
    all_d <- c()
    for (mode_here in modes) all_d <- c(all_d, dnorm(x, mean = mode_here, sd = sigma))
    return(mean(all_d))
  }
  foo2 <- Vectorize(function(x) dmixture(x = x, modes = modes, sigma = sigma))
  foo <- Vectorize(function(x) dmixture(x = x, modes = modes, sigma = sigma)^2)
  return(list(x = x, f_d = foo2, f_dsquare = foo))
}

fit_density_hal_one_setting <- function(n_sim, n_mode, bin_width) {
  data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
  # specify grid of lambda for HAL fit
  lambda_grid <- 10^seq(-1, -5, by = -1)
  longDataOut <- longiData$new(x = data_out$x, bin_width = bin_width)
  # tune HAL for density
  cvHAL_fit <- cv_densityHAL$new(x = data_out$x, longiData = longDataOut)
  cvHAL_fit$assign_fold(n_fold = 4)
  cvHAL_fit$cv_lambda_grid(lambda_grid = lambda_grid)
  hal_out <- cvHAL_fit$compute_model_full_data(cvHAL_fit$lambda.min)
  # normalize the fit to plot
  density_intial <- empiricalDensity$new(
    p_density = hal_out$predict(new_x = data_out$x),
    x = data_out$x
  )
  density_intial$normalize()
  # plot the fitted density on top of true density
  density_intial$display(p_truth = data_out$f_d)
}

n_sim <- 1e3
n_modes <- seq(1, 10, by = 1)
bin_width <- 1e-1
for (n_mode in n_modes) {
  fit_density_hal_one_setting(n_sim, n_mode, bin_width)
}
