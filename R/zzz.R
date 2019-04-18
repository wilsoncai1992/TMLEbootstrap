.onAttach <- function(...) {
  packageStartupMessage(paste(
    "TMLEbootstrap: Bootstrap Confidence Intervals For Targeted Maximum Likelihood Estimators"
  ))
  packageStartupMessage(
    "Version: ",
    utils::packageDescription("TMLEbootstrap")$Version
  )
}
