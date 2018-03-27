data(ssm560_reduced)
library("BSgenome.Hsapiens.1000genomes.hs37d5")
bsg = BSgenome.Hsapiens.1000genomes.hs37d5
data(mutation_categories)

context("import.data")

test_that("import.data produces correct data", {
    expect_output(import.data(input=ssm560_reduced,bsg=bsg,mutation_categories=mutation_categories),"import.data")
})
