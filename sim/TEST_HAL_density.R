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

tmleOut <- avgDensityTMLE$new(x = datO$sA, verbose = FALSE)
tmleOut$fit_density()
foo2 <- function(x) {(.5*dnorm(x, mean = 2) + .5*dnorm(x, mean = -2))}
tmleOut$p_hat$display(foo2)
tmleOut$calc_Psi()
tmleOut$calc_EIC()
tmleOut$Psi

bootOut <- avgDensityBootstrap$new(x = x)
bootOut$bootstrap(REPEAT_BOOTSTRAP = 1e3)
bootOut$CI_all

foo <- function(x) {(.5*dnorm(x, mean = 2) + .5*dnorm(x, mean = -2))^2}
true_Psi <-sum(foo(seq(-10, 10, 1e-3))*1e-3)
true_Psi
