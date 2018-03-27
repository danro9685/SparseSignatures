data(cv_example)
data(lambda_range_example)
data(nmf_LassoK_example)

context("as.mean.squared.error")

test_that("as.mean.squared.error produces correct output", {
    expect_equal(length(as.mean.squared.error(cv_example)),2)
})

context("as.alpha.in.range")

test_that("as.alpha.in.range produces correct output", {
    expect_equal(length(as.alpha.in.range(lambda_range_example,lambda_value=0.10)),6160)
})

context("as.beta.in.range")

test_that("as.beta.in.range produces correct output", {
    expect_equal(length(as.beta.in.range(lambda_range_example,lambda_value=0.10)),1056)
})

context("as.starting.beta.in.range")

test_that("as.starting.beta.in.range produces correct output", {
    expect_equal(length(as.starting.beta.in.range(lambda_range_example,lambda_value=0.10)),1056)
})

context("as.loglik.progression.in.range")

test_that("as.loglik.progression.in.range produces correct output", {
    expect_equal(length(as.loglik.progression.in.range(lambda_range_example,lambda_value=0.10)),20)
})

context("as.alpha")

test_that("as.alpha produces correct output", {
    expect_equal(length(as.alpha(nmf_LassoK_example)),3360)
})

context("as.beta")

test_that("as.beta produces correct output", {
    expect_equal(length(as.beta(nmf_LassoK_example)),576)
})

context("as.starting.beta")

test_that("as.starting.beta produces correct output", {
    expect_equal(length(as.starting.beta(nmf_LassoK_example)),576)
})

context("as.loglik.progression")

test_that("as.loglik.progression produces correct output", {
    expect_equal(length(as.loglik.progression(nmf_LassoK_example)),5)
})
