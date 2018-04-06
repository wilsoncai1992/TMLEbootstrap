################################################################################
# setup
################################################################################
# set.seed(628957)
# library(adaptest)
library(truncnorm)
context("ateBootstrap results should not be NA")
simulate_data <- function(n_sim, a1, a2, b1) {
  thresholding <- function(x, min, max) pmin(pmax(x, min), max)

  W <- truncnorm::rtruncnorm(n = n_sim, a = -10, b = 10, mean = 0, sd = 4)
  A <- rbinom(n_sim, size = 1, prob = thresholding(.3 + 0.1*W*sin(a2*W), 0.3, 0.7) + rnorm(n_sim, mean = 0, sd = 0.05))

  # Y <- 0.05*W^2 + b1*sin(W*a1) + A + rnorm(n_sim, 0, 1)
  Y <- b1*sin(W*a1) + A + rnorm(n_sim, 0, 1)

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
n_sim = 2e2
a1 = 1
b1 = 3
a2 = .1
INFLATE_LAMBDA = 1

data_sim <- simulate_data(n_sim = n_sim, a1 = a1, a2 = a2, b1 = b1)
bootOut <- ateBootstrap$new(data = data_sim)
bootOutExact <- bootOut$clone(deep = TRUE) # deep copy point tmle, less repeat

bootOut$bootstrap(2e2)
bootOutExact$exact_bootstrap(2e2)
# combine exact boot with existing results
regularCI <- bootOut$all_boot_CI()
taylorCI <- bootOutExact$all_boot_CI()
out <- list(wald = regularCI$wald,
            reg = regularCI$boot,
            reg_pen = regularCI$penalized,
            reg_scale = regularCI$scale,
            reg_scale_pen = regularCI$scale_penalized,
            taylor = taylorCI$boot,
            taylor_pen = taylorCI$penalized,
            taylor_scale = taylorCI$scale,
            taylor_scale_pen = taylorCI$scale_penalized)
################################################################################
test_that("ateBootstrap results should not be NA", {
  expect_true(all(!sapply(out, is.na)))
})

