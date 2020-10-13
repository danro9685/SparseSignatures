context("SparseSignatures")

data("nmf_LassoK_example")
test_that("SparseSignatures produces correct output", {
    expect_equal(names(nmf_LassoK_example),c("alpha","beta","starting_beta","best_loglik","loglik_progression"))
})
