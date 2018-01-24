#' @export
avgDensityBootstrap <- R6Class("avgDensityBootstrap",
  public = list(
    x = NULL,
    pointTMLE = NULL,
    Psi = NULL,
    # EIC = NULL,
    epsilon_step = 1e-2,
    # tol = 1e-3,
    CI_all = NULL,
    bootstrap_estimates = NULL,
    initialize = function(x, epsilon_step = NULL) {
      self$x <- x
      if(!is.null(epsilon_step)) self$epsilon_step <- epsilon_step
      onestepFit <- avgDensityTMLE$new(x = self$x, epsilon_step = self$epsilon_step, verbose  = TRUE)
      onestepFit$fit_density()
      onestepFit$calc_Psi()
      onestepFit$calc_EIC()
      onestepFit$onestepTarget()
      onestepFit$inference()

      self$pointTMLE <- onestepFit
      self$Psi <- onestepFit$Psi
    },
    bootstrap = function(REPEAT_BOOTSTRAP = 2e2){
      browser()
      SAMPLE_PER_BOOTSTRAP <- length(self$x)

      betfun <- function(data, epsilon_step = self$epsilon_step){
        # indices is the random indexes for the bootstrap sample
        indices <- sample(1:length(data), size = SAMPLE_PER_BOOTSTRAP, replace = TRUE) # user specify sample size
        d = data[indices]

        bootstrapOnestepFit <- avgDensityTMLE$new(x = d, epsilon_step = epsilon_step)
        # fit new density
        HAL_wrapper <<- generate_SL.fixed_HAL(hal9001_object = self$pointTMLE$HAL_tuned)
        HAL_wrapper
        longDFOut_new <- self$pointTMLE$longDataOut$generate_df(x = d)
        SL_fit <- SuperLearner(Y = longDFOut_new$Y, X = longDFOut_new[,'box',F], newX = data.frame(box = d),
                               family = 'binomial',
                               SL.library = "HAL_wrapper",
                               cvControl = list(V = 3),
                               verbose = verbose)
        density_intial <- empiricalDensity$new(p_density = SL_fit$SL.predict, x = d)
        bootstrapOnestepFit$p_hat <- density_intial$normalize()
        # target new fit
        bootstrapOnestepFit$calc_Psi()
        bootstrapOnestepFit$calc_EIC()
        bootstrapOnestepFit$onestepTarget()
        # bootstrapOnestepFit$inference()

        return(c(bootstrapOnestepFit$Psi))
      }

      library(foreach)

      # library(doSNOW)
      # library(tcltk)
      # nw <- parallel:::detectCores()  # number of workers
      # cl <- makeSOCKcluster(nw)
      # registerDoSNOW(cl)

      # library(Rmpi)
      # library(doMPI)
      # cl = startMPIcluster()
      # registerDoMPI(cl)
      # clusterSize(cl) # just to check


      all_bootstrap_estimates <- foreach(it2 = 1:(REPEAT_BOOTSTRAP), .combine = c,
                                         .inorder = FALSE,
                                         .packages = c('R6'),
                                         .errorhandling = 'remove',
                                         .export = c('self'),
                                         .verbose = F) %do% {
                                         # .verbose = T) %dopar% {
        # source("./averageDensityValueTMLE.R")
        # source("./empiricalDensityR6.R")
        if(it2 %% 10 == 0) print(it2)
        betfun(self$x, self$epsilon_step)
      }

      # save(all_bootstrap_estimates, file = 'all_bootstrap_estimates.rda')
      # closeCluster(cl)
      stopCluster(cl)

      ALPHA <- 0.05
      # remove errors
      if( !all(sapply(all_bootstrap_estimates, class) == 'numeric') ) message(paste('Error happens.', sum(sapply(all_bootstrap_estimates, class) == 'numeric'), 'bootstraps are correct'))
      all_bootstrap_estimates <- as.numeric(all_bootstrap_estimates[sapply(all_bootstrap_estimates, class) == 'numeric'])
      self$bootstrap_estimates <- all_bootstrap_estimates

      boot1_CI <- quantile(all_bootstrap_estimates, probs = c(ALPHA/2, 1 - ALPHA/2))
      normal_CI <- self$pointTMLE$CI
      self$CI_all <- list(normal_CI, boot1_CI)
    },
    wider_boot_CI = function(){
      bootCI_minus_Psi <- self$CI_all[[2]] - self$Psi
      new_dist <- max(abs(bootCI_minus_Psi))
      new_upper <- self$Psi + new_dist
      new_lower <- self$Psi - new_dist
      new_CI <- c(new_lower, new_upper)
      return(list(self$CI_all[[1]], self$CI_all[[2]], new_CI))
    }
  )
)
