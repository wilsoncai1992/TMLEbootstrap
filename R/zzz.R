.onAttach <- function(...) {
  packageStartupMessage(paste("fixedHAL: TEMP",
                              "\n TEMP"))
  packageStartupMessage("Version: ",
                        utils::packageDescription("fixedHAL")$Version)
}

