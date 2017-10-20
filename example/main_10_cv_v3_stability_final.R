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
load("data/allbg.RData")
background_genome = allbg$freq

# set the number of signatures and lambda to be considered
K = 3:9
lambda_values = c(0.10,0.15)
nmf_runs = 100
cross_validation_entries = 0.10
cross_validation_iterations = 10
num_processes = 48
num_iterations = 100
my_seed_starting_beta = 28954

# ANALYSIS USING THE GERMLINE SIGNATURE FROM THE PAPER AS BACKGROUND

# fit the initial betas for each configuration
initial_betas_germline_cv_10 = startingBetasEstimation(x=patients,K=K,nmf_runs=nmf_runs,seed=my_seed_starting_beta,verbose=TRUE)
save(initial_betas_germline_cv_10,file="data/initial_betas_germline_cv_10.RData")

# fit the signatures with the given background noise model
stability_germline_cv_10 = list()
for(i in 1:num_iterations) {
    signatures_nmfLasso_germline = nmfLasso(x=patients,K=K,starting_beta=initial_betas_germline_cv_10,background_signature=background_genome,nmf_runs=nmf_runs,lambda_values=lambda_values,cross_validation_entries=cross_validation_entries,cross_validation_iterations=cross_validation_iterations,iterations=20,max_iterations_lasso=10000,num_processes=num_processes,seed=(i*1000*1),verbose=TRUE)
    stability_germline_cv_10[[i]] = signatures_nmfLasso_germline[["best_configuration"]]
}
save(stability_germline_cv_10,file="data/stability_germline_cv_10.RData")
