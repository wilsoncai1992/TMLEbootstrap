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
longDataOut <- longiData$new(x = datO$sA, bin_width = .1)
longDFOut <- longDataOut$longiDF


verbose <- FALSE
library(SuperLearner)
library(hal9001)
# tune HAL for density
SL_fit <- SuperLearner(Y = longDFOut$Y, X = longDFOut[,'box',F], newX = data.frame(box = longDataOut$x),
                       family = 'binomial',
                       SL.library = "SL.hal9001",
                       cvControl = list(V = 3),
                       verbose = verbose)
Q_HAL_tuned <- SL_fit$fitLibrary$SL.hal9001_All$object
Q_HAL_tuned <- squash_hal_fit(Q_HAL_tuned)

density_intial <- empiricalDensity$new(p_density = SL_fit$SL.predict, x = x)
density_intial$normalize()
foo2 <- function(x) {(.5*dnorm(x, mean = 2) + .5*dnorm(x, mean = -2))}
density_intial$display(foo2)
