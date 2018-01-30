library("simcausal")
D <- DAG.empty()
D <-
  D + node("W1", distr = "rbern", prob = 0.5) +
  node("sA.mu", distr = "rconst", const = 4*W1 - 2) +
  node("sA", distr = "rnorm", mean = sA.mu, sd = 1)
D <- set.DAG(D, n.test = 10)
# datO <- sim(D, n = 1e4, rndseed = 12345)
datO <- sim(D, n = 1e3, rndseed = 12345)

x <- datO$sA

bin_width <- .3
longDataOut <- longiData$new(x = x, bin_width = bin_width)
longDFOut <- longDataOut$generate_df_compress()


cvHAL_fit <- cv_densityHAL$new(x = x, longiData = longDataOut)
cvHAL_fit$assign_fold(n_fold = 3)
# cvHAL_fit$cv(lambda = 2e-5)
cvHAL_fit$cv_lambda_grid(lambda_grid = c(1e-6,2e-5))
# longDFResample <- longiData_resample$new(df_compressed = longDFOut)
# df_boot_with_replace <- longDFResample$bootstrap_with_replacement()
# folds_out <- longDFResample$create_folds(n_fold = 10)

# summed <- sum_longiData_resample(folds_out)
# all.equal(summed$df_compressed, longDFResample$df_compressed)
# all.equal(summed, longDFResample)

# train_1 <- sum_longiData_resample(folds_out[1:2])
# train_1$n_sample

# for(fold in 1:length(folds_out)){
#   train_here <- sum_longiData_resample(folds_out[which_are_train])
#   # train_here$n_sample
#   test_here <- sum_longiData_resample(folds_out[fold])
#   # test_here$n_sample
#   
#   cvHAL$new()
# }

# hal_fit <- hal9001::fit_hal_single_lambda(X = longDFOut[,'box'],
#   Y = longDFOut$Y,
#   weights = longDFOut$Freq,
#   family="binomial",
#   lambda = 2e-5,
#   fit_type = 'glmnet',
#   use_min = TRUE, #useless
#   yolo = FALSE)
# yhat <- rje::expit(predict(hal_fit, new_data = longDataOut$x))
# density_intial <- empiricalDensity$new(p_density = yhat, x = x)
# p_hat <- density_intial$normalize()
# foo2 <- function(x) {(.5*dnorm(x, mean = 2) + .5*dnorm(x, mean = -2))}
# p_hat$display(foo2)




















HAL_tuned <- hal9001::fit_hal(X = longDFOut[,'box'],
  Y = longDFOut$Y,
  weights = longDFOut$Freq,
  use_min = TRUE,
  yolo = FALSE,
  fit_type = 'glmnet',
  family = "binomial",
  n_folds = 3)
yhat <- rje::expit(predict(HAL_tuned, new_data = self$longDataOut$x))
density_intial <- empiricalDensity$new(p_density = yhat, x = self$x)
self$p_hat <- density_intial$normalize()


tmleOut <- avgDensityTMLE$new(x = datO$sA, verbose = FALSE)
tmleOut$fit_density(bin_width = .3)
foo2 <- function(x) {(.5*dnorm(x, mean = 2) + .5*dnorm(x, mean = -2))}
tmleOut$p_hat$display(foo2)
tmleOut$calc_Psi()
tmleOut$calc_EIC()
tmleOut$Psi

bootOut <- avgDensityBootstrap$new(x = x)
bootOut$bootstrap(REPEAT_BOOTSTRAP = 1e2)
bootOut$CI_all

foo <- function(x) {(.5*dnorm(x, mean = 2) + .5*dnorm(x, mean = -2))^2}
true_Psi <-sum(foo(seq(-10, 10, 1e-3))*1e-3)
true_Psi
