################################################################################
# setup
################################################################################
# set.seed(628957)
context("avgDensityBootstrap results should not be NA")
simulate_data <- function(n_sim, n_mode) {
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
################################################################################
# simulation
################################################################################
n_sim <- 5e2
n_mode <- 2
bin_width <- 5e-1

data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
bootstrapFit <- avgDensityBootstrap$new(x = data_out$x, bin_width = bin_width)
bootstrapFit$bootstrap(n_bootstrap = 2e1)

bootstrapFitExact <- bootstrapFit$clone(deep = TRUE)
bootstrapFitExact$exact_bootstrap(n_bootstrap = 2e1)
bootstrapFitExact_paper <- bootstrapFit$clone(deep = TRUE)
bootstrapFitExact_paper$exact_bootstrap_paper(n_bootstrap = 2e1)

regularCI <- bootstrapFit$all_boot_CI()
exact_CI <- bootstrapFitExact$all_boot_CI()
exact_CI_paper <- bootstrapFitExact_paper$all_boot_CI()
out <- list(
  wald = regularCI$wald,
  reg = regularCI$boot,
  reg_pen = regularCI$penalized,
  reg_scale = regularCI$scale,
  reg_scale_pen = regularCI$scale_penalized,
  taylor = exact_CI$boot,
  taylor_pen = exact_CI$penalized,
  taylor_scale = exact_CI$scale,
  taylor_scale_pen = exact_CI$scale_penalized,
  taylor2 = exact_CI_paper$boot,
  taylor2_pen = exact_CI_paper$penalized,
  taylor2_scale = exact_CI_paper$scale,
  taylor2_scale_pen = exact_CI_paper$scale_penalized
)

# without targeting
bootOut_HALMLE <- avgDensityBootstrap$new(x = data_out$x, bin_width = bin_width, targeting = FALSE)
bootOut_HALMLE$bootstrap(2e1)
halmleCI <- bootOut_HALMLE$all_boot_CI()


CVOut <- comprehensiveBootstrap$new(
  parameter = avgDensityBootstrap,
  x = data_out$x,
  bin_width = bin_width,
  lambda_grid = NULL,
  epsilon_step = 1e-1
)
CVOut$bootstrap(n_bootstrap = 2e1)
CVOut$all_CI()

test_that("avgDensityBootstrap results should not be NA", {
  expect_true(all(!sapply(out, is.na)))
})

test_that("HAL-MLE bootstrap results should not be NA", {
  expect_true(all(!sapply(halmleCI, is.na)))
})

test_that("HALselect some beta", {
  expect_true(!is.null(bootOut_HALMLE$pointTMLE$compute_min_phi_ratio()))
})

test_that("comprehensiveBootstrap intervals working", {
  expect_true(all(!sapply(CVOut$CI_all, is.na)))
})

################################################################################
tune_param_fit <- avgDensityTuneHyperparam$new(
  data = data_out, bin_width = bin_width, epsilon_step = 1e-2, n_bootstrap = 2e2
)
tune_param_fit$add_lambda(lambda_grid = 10 ^ seq(0, -3, length.out = 2))
test_that("plateau tuning parameter is working", {
  expect(!is.null(
    tune_param_fit$select_lambda_pleateau_wald(tune_param_fit$get_lambda_df())
  ))
})
# tune_param_fit$plot_CI()
tune_param_fit$plot_width(type_CI = "wald")
