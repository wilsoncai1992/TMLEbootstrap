.onAttach <- function(...) {
  packageStartupMessage(paste(
    "TMLEbootstrap: ",
    "\n bootstrap confidence intervals for Targeted Maximum Likelihood Estimators"
  ))
  packageStartupMessage(
    "Version: ",
    utils::packageDescription("TMLEbootstrap")$Version
  )
}
