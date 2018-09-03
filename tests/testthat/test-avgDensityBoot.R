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
n_sim <- 1e2
n_mode <- 1
bin_width <- 8e-1

data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
bootstrapFit <- avgDensityBootstrap$new(x = data_out$x, bin_width = bin_width)
bootstrapFitExact <- bootstrapFit$clone(deep = TRUE) # deep copy point tmle, less repeat

bootstrapFit$bootstrap(REPEAT_BOOTSTRAP = 2e2)
bootstrapFitExact$exact_bootstrap(REPEAT_BOOTSTRAP = 2e2)
regularCI <- bootstrapFit$all_boot_CI()
taylorCI <- bootstrapFitExact$all_boot_CI()

out <- list(
  wald = regularCI$wald,
  reg = regularCI$boot,
  reg_pen = regularCI$penalized,
  reg_scale = regularCI$scale,
  reg_scale_pen = regularCI$scale_penalized,
  taylor = taylorCI$boot,
  taylor_pen = taylorCI$penalized,
  taylor_scale = taylorCI$scale,
  taylor_scale_pen = taylorCI$scale_penalized
)

# without targeting
bootOut_HALMLE <- avgDensityBootstrap$new(x = data_out$x, bin_width = bin_width, targeting = FALSE)
bootOut_HALMLE$bootstrap(2e1)
halmleCI <- bootOut_HALMLE$all_boot_CI()
################################################################################
test_that("avgDensityBootstrap results should not be NA", {
  expect_true(all(!sapply(out, is.na)))
})

test_that("HAL-MLE bootstrap results should not be NA", {
  expect_true(all(!sapply(halmleCI, is.na)))
})

test_that("HALselect some beta", {
  expect_true(!is.null(bootOut_HALMLE$pointTMLE$compute_min_phi_ratio()))
})
