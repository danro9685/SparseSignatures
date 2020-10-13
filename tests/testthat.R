Sys.setenv("R_TESTS" = "")

library("testthat")
library("SparseSignatures")

test_check("SparseSignatures")
