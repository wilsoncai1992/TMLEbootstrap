################################################################################
# setup
################################################################################
# set.seed(628957)
# library(adaptest)
context("blipVarianceBoot results should not be NA")
simulate_data <- function(n_sim, n_mode = 1, a1 = 2) {
  dmixture <- function(x, modes, sigma) {
    all_d <- c()
    for(mode_here in modes) all_d <- c(all_d, dnorm(x, mean = mode_here, sd = sigma))
    return(mean(all_d))
  }
  modes <- seq(from = -2, to = 2, length.out = n_mode)
  sigma <- 10/n_mode/8
  fW <- Vectorize(function(x) a1*dmixture(x = x, modes = modes, sigma = sigma)) # vary magnitude

  A <- rbinom(n = n_sim, size = 1, prob = .5)
  W <- runif(n = n_sim, min = -4, max = 4)
  Y <- rnorm(n = n_sim, mean = A*fW(W), sd = .1)
  df <- list(Y = Y, A = A, W = data.frame(W))

  Y1 <- fW(W)
  Y0 <- 0
  mean(Y1 - Y0) # this is true ATE
  psi_true <- var(Y1 - Y0) # this is true blip variance
  psi_true
  return(list(df = df, psi_true = psi_true, fW = fW))
}
################################################################################
# simulation
################################################################################
n_sim = 2e2
n_mode = 1
a1 = 2

data <- simulate_data(n_sim = n_sim, n_mode = n_mode, a1 = a1)
df <- data$df
bootstrapFit <- blipVarianceBootstrap_contY$new(data = df)
bootstrapFitExact <- bootstrapFit$clone(deep = TRUE) # deep copy point tmle, less repeat

bootstrapFit$bootstrap(REPEAT_BOOTSTRAP = 2e2)
bootstrapFitExact$exact_bootstrap(REPEAT_BOOTSTRAP = 2e2)
regularCI <- bootstrapFit$all_boot_CI()
taylorCI <- bootstrapFitExact$all_boot_CI()

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
test_that("blipVarianceBoot results should not be NA", {
  expect_true(all(!sapply(out, is.na)))
})

