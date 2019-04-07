.onAttach <- function(...) {
  packageStartupMessage(paste(
    "TMLEbootstrap: bootstrap confidence intervals for Targeted Maximum Likelihood Estimators",
    "\n"
  ))
  packageStartupMessage(
    "Version: ",
    utils::packageDescription("TMLEbootstrap")$Version
  )
}
