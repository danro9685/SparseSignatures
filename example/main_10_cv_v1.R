# set the working directory
work.dir = "/Users/daniele/Documents/Stanford-github/SparseSignatures"
setwd(work.dir)

# load the required libraries and sources
library("NMF")
library("nnlasso")
library("nnls")
library("parallel")
source("R/mutation_factor.R")

# libraries required to plot the extracted signatures
library("data.table")
library("ggplot2")
library("gridExtra")

# function to plot extracted signatures
"plotSignatures" <- function( beta, patients_ids = colnames(beta), backgroundFreq = TRUE ) {
    
    # set rownames and colnames
    if(backgroundFreq) {
        rownames(beta) = c("Background",paste0("S",1:(nrow(beta)-1)))
    }
    else {
        rownames(beta) = paste0("S",1:nrow(beta))
    }
    colnames(beta) = patients_ids
    
    # make the plot
    x = as.data.table(melt(beta))
    x[,context := paste0(substr(Var2,1,1),".",substr(Var2,7,7))]
    x[,alt := paste0(substr(Var2,3,3),">",substr(Var2,5,5))]
    
    glist = list()
    
    for(i in 1:nrow(beta)) {
        sig = rownames(beta)[i]
        glist[[i]] = ggplot(x[Var1 == sig]) + 
        geom_bar(aes(x = context, y = value, fill = alt), stat = "identity", position = "identity") + 
        facet_grid(.~alt) + 
        theme(legend.position="none", 
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            panel.background = element_rect(fill = "white", colour = NA),
            panel.grid.major = element_line(colour = "grey80"),
            panel.grid.minor = element_blank()) + 
        ggtitle(sig)
    }
    grid.arrange(grobs=glist, ncol=ceiling(nrow(beta)/4))
    
}

# load the data
load("data/patients.RData")
load("data/genome.RData")
load("data/allbg.RData")

# set the number of signatures and lambda to be considered
K = 2:15
lambda_values = c(0.01,0.05,0.10,0.15)
cross_validation_entries = 0.10
cross_validation_iterations = 5
num_processes = 10
my_seed_starting_beta = 93211
my_seed_nmfLasso = 74833

# ANALYSIS USING THE GENOME FREQUENCY AS BACKGROUND

# fit the initial betas for each configuration
initial_betas_genome = startingBetasEstimation(x=patients,K=K,nmf_runs=10,seed=my_seed_starting_beta,verbose=TRUE)
save(initial_betas_genome,file="data/initial_betas_genome_10.RData")

# fit the signatures with the given background noise model
signatures_nmfLasso_genome = nmfLasso(x=patients,K=K,starting_beta=initial_betas_genome,background_signature=genome$freq,nmf_runs=10,lambda_values=lambda_values,cross_validation_entries=cross_validation_entries,cross_validation_iterations=cross_validation_iterations,iterations=20,max_iterations_lasso=10000,num_processes=num_processes,seed=my_seed_nmfLasso,verbose=TRUE)
save(signatures_nmfLasso_genome,file="data/signatures_nmfLasso_genome_10.RData")

# plot the resulting signatures
initial_betas_genome = rbind(signatures_nmfLasso_genome$best_configuration$background_signature,signatures_nmfLasso_genome$best_configuration$starting_beta)
initial_betas_genome = initial_betas_genome / rowSums(initial_betas_genome)
plotSignatures(initial_betas_genome,patients_ids=colnames(patients),backgroundFreq=TRUE)
plotSignatures(signatures_nmfLasso_genome$best_configuration$beta,patients_ids=colnames(patients),backgroundFreq=TRUE)

# plot the log-likelihood values
plot(signatures_nmfLasso_genome$best_configuration$loglik_progression)

# ANALYSIS USING THE GERMLINE SIGNATURE AS BACKGROUND

# fit the initial betas for each configuration
initial_betas_germline = initial_betas_genome
save(initial_betas_germline,file="data/initial_betas_germline_10.RData")

# fit the signatures with the given background noise model
signatures_nmfLasso_germline = nmfLasso(x=patients,K=K,starting_beta=initial_betas_germline,background_signature=allbg$freq,nmf_runs=10,lambda_values=lambda_values,cross_validation_entries=cross_validation_entries,cross_validation_iterations=cross_validation_iterations,iterations=20,max_iterations_lasso=10000,num_processes=num_processes,seed=my_seed_nmfLasso,verbose=TRUE)
save(signatures_nmfLasso_germline,file="data/signatures_nmfLasso_germline_10.RData")

# plot the resulting signatures
initial_betas_germline = rbind(signatures_nmfLasso_germline$best_configuration$background_signature,signatures_nmfLasso_germline$best_configuration$starting_beta)
initial_betas_germline = initial_betas_germline / rowSums(initial_betas_germline)
plotSignatures(initial_betas_germline,patients_ids=colnames(patients),backgroundFreq=TRUE)
plotSignatures(signatures_nmfLasso_germline$best_configuration$beta,patients_ids=colnames(patients),backgroundFreq=TRUE)

# plot the log-likelihood values
plot(signatures_nmfLasso_germline$best_configuration$loglik_progression)
