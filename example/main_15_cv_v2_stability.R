# set the working directory
work.dir = "/Users/daniele/Documents/Stanford-github/SparseSignatures"
setwd(work.dir)

# load the required libraries and sources
library("NMF")
library("nnlasso")
library("nnls")
library("parallel")
source("R/mutation_factor.R")

# load the data
load("data/patients.RData")
load("data/genome.RData")
load("data/allbg.RData")

# set the number of signatures and lambda to be considered
K = 3:9
lambda_values = c(0.10,0.15)
cross_validation_entries = 0.15
cross_validation_iterations = 5
num_processes = 15
num_iterations = 100

# ANALYSIS USING THE GERMLINE SIGNATURE AS BACKGROUND

# load the previous results with the same configurations
load(file="data/signatures_nmfLasso_germline_15.RData")
initial_betas = signatures_nmfLasso_germline[["starting_beta"]][paste0(K,"_signatures"),,drop=FALSE]

# fit the signatures with the given background noise model
stability_cv_15 = list()
for(i in 1:num_iterations) {
    signatures_nmfLasso_germline = nmfLasso(x=patients,K=K,starting_beta=initial_betas,background_signature=allbg$freq,nmf_runs=10,lambda_values=lambda_values,cross_validation_entries=cross_validation_entries,cross_validation_iterations=cross_validation_iterations,iterations=20,max_iterations_lasso=10000,num_processes=num_processes,seed=(i*1000*2),verbose=TRUE)
    stability_cv_15[[i]] = signatures_nmfLasso_germline[["best_configuration"]]
}
save(stability_cv_15,file="data/stability_cv_15.RData")
