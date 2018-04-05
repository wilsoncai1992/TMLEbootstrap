library(testthat)
library(fixedHAL)

Sys.setenv(R_TESTS = "")
test_check("fixedHAL")
