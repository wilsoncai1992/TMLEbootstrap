################################################################################
# setup
################################################################################
# set.seed(628957)
context("blip Variance bootstrap results should not be NA")
simulate_data <- function(n_sim, n_mode = 1, a1 = 2) {
  dmixture <- function(x, modes, sigma) {
    all_d <- c()
    for (mode_here in modes) all_d <- c(all_d, dnorm(x, mean = mode_here, sd = sigma))
    return(mean(all_d))
  }
  modes <- seq(from = -2, to = 2, length.out = n_mode)
  sigma <- 10 / n_mode / 8
  fW <- Vectorize(function(x) a1 * dmixture(x = x, modes = modes, sigma = sigma))
  f2W <- Vectorize(function(x) (a1 * dmixture(x = x, modes = modes, sigma = sigma))^2)

  A <- rbinom(n = n_sim, size = 1, prob = .5)
  W <- runif(n = n_sim, min = -4, max = 4)
  Y <- rnorm(n = n_sim, mean = A * fW(W), sd = .1)
  df <- list(Y = Y, A = A, W = data.frame(W))

  Y1 <- fW(W)
  Y0 <- 0
  mean(Y1 - Y0) # this is true ATE
  psi_true <- var(Y1 - Y0) # this is true blip variance
  psi_true

  E_blip <- mean(fW(seq(-4, 4, by = 1e-3)))
  E2_blip <- mean(f2W(seq(-4, 4, by = 1e-3)))
  true_var_blip <- E2_blip - E_blip^2
  return(list(df = df, psi_true = true_var_blip, fW = fW))
}

################################################################################
# simulation
################################################################################
n_sim <- 2e2
n_mode <- 1
a1 <- 2

data <- simulate_data(n_sim = n_sim, n_mode = n_mode, a1 = a1)
# data$psi_true
df <- data$df
bootstrapFit <- blipVarianceBootstrapContinuousY$new(data = df)
bootstrapFit$bootstrap(n_bootstrap = 2e1)

bootstrapFitExact <- bootstrapFit$clone(deep = TRUE)
bootstrapFitExact$exact_bootstrap(n_bootstrap = 2e1)
bootstrapFitExact_paper <- bootstrapFit$clone(deep = TRUE)
bootstrapFitExact_paper$exact_bootstrap_paper(n_bootstrap = 2e1)

regular_ci <- bootstrapFit$all_boot_CI()
exact_ci <- bootstrapFitExact$all_boot_CI()
exact_ci_paper <- bootstrapFitExact_paper$all_boot_CI()

out <- list(
  wald = regular_ci$wald,
  reg = regular_ci$boot,
  reg_pen = regular_ci$penalized,
  reg_scale = regular_ci$scale,
  reg_scale_pen = regular_ci$scale_penalized,
  taylor = exact_ci$boot,
  taylor_pen = exact_ci$penalized,
  taylor_scale = exact_ci$scale,
  taylor_scale_pen = exact_ci$scale_penalized,
  taylor2 = exact_ci_paper$boot,
  taylor2_pen = exact_ci_paper$penalized,
  taylor2_scale = exact_ci_paper$scale,
  taylor2_scale_pen = exact_ci_paper$scale_penalized
)

# without targeting
bootOut_HALMLE <- blipVarianceBootstrapContinuousY$new(data = df, targeting = FALSE)
bootOut_HALMLE$bootstrap(2e1)
halmleCI <- bootOut_HALMLE$all_boot_CI()

CVOut <- comprehensiveBootstrap$new(
  parameter = blipVarianceBootstrapContinuousY,
  data = df
)
CVOut$bootstrap(n_bootstrap = 2e1)
CVOut$all_CI()


################################################################################
test_that("blipVarianceBoot results should not be NA", {
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
