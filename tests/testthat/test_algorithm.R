data(starting_betas_example)
beta = starting_betas_example[["5_signatures","Value"]]

context("nmf.LassoK")

test_that("nmf.LassoK produces correct output", {
    expect_equal(names(nmf.LassoK(x=patients,K=5,beta=beta,background=background,lambda_rate=0.10,iterations=5,num_processes=NA)),c("alpha","beta","starting_beta","best_loglik","loglik_progression"))
})
