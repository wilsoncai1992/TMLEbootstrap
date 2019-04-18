################################################################################
# setup
################################################################################
# set.seed(628957)
library(truncnorm)
context("ateBootstrap results should not be NA")
simulate_data <- function(n_sim, a1, a2, b1) {
  thresholding <- function(x, min, max) pmin(pmax(x, min), max)

  W <- truncnorm::rtruncnorm(n = n_sim, a = -10, b = 10, mean = 0, sd = 4)
  A <- rbinom(
    n_sim,
    size = 1,
    prob = thresholding(.3 + 0.1 * W * sin(a2 * W), 0.3, 0.7) +
      rnorm(n_sim, mean = 0, sd = 0.05)
  )
  Y <- b1 * sin(W * a1) + A + rnorm(n_sim, 0, 1)

  X_matrix_0 <- data.frame(A, W)
  all_df <- data.frame(Y, A, W)

  # append one value of Z
  Z <- rep(1, n_sim)
  X_matrix <- cbind(Z, X_matrix_0)
  all_df <- cbind(Z, all_df)
  return(list(Y = Y, A = A, W = data.frame(W), X_matrix = X_matrix, all_df = all_df))
}
################################################################################
# simulation
################################################################################
n_sim <- 2e2
a1 <- 1
b1 <- 3
a2 <- .1
INFLATE_LAMBDA <- 1

data_sim <- simulate_data(n_sim = n_sim, a1 = a1, a2 = a2, b1 = b1)
boot_output <- ateBootstrap$new(data = data_sim)
boot_output$bootstrap(2e1)

boot_output_exact <- boot_output$clone(deep = TRUE)
boot_output_exact$exact_bootstrap(2e1)
boot_output_exact_paper <- boot_output$clone(deep = TRUE)
boot_output_exact_paper$exact_bootstrap_paper(2e1)

# combine exact boot with existing results
regular_ci <- boot_output$all_boot_CI()
exact_ci <- boot_output_exact$all_boot_CI()
exact_ci_paper <- boot_output_exact_paper$all_boot_CI()
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
boot_output_HALMLE <- ateBootstrap$new(data = data_sim, targeting = FALSE)
boot_output_HALMLE$bootstrap(2e1)
halmleCI <- boot_output_HALMLE$all_boot_CI()

# comprehensive bootstrap
CVOut <- comprehensiveBootstrap$new(parameter = ateBootstrap, data = data_sim)
CVOut$bootstrap(2e1)
CVOut$all_CI()

################################################################################
test_that("ateBootstrap results should not be NA", {
  expect_true(all(!sapply(out, is.na)))
})

test_that("HAL-MLE bootstrap results should not be NA", {
  expect_true(all(!sapply(halmleCI, is.na)))
})

test_that("HALselect some beta", {
  expect_true(!is.null(boot_output_HALMLE$pointTMLE$compute_min_phi_ratio()))
})

test_that("comprehensiveBootstrap intervals working", {
  expect_true(all(!sapply(CVOut$CI_all, is.na)))
})
