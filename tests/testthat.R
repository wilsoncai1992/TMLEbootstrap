library(testthat)
library(TMLEbootstrap)

Sys.setenv(R_TESTS = "")
test_check("TMLEbootstrap")
