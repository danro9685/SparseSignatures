data(starting_betas_example)
beta = starting_betas_example[["5_signatures","Value"]]

context("nmfLasso")

test_that("nmfLasso produces correct output", {
    expect_equal(names(nmfLasso(x=patients,K=5,starting_beta=beta,background_signature=background,lambda_rate_alpha=0.00,lambda_rate_beta=0.05,iterations=5,num_processes=NA)),c("alpha","beta","starting_alpha","starting_beta","loglik_progression"))
})
