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

# set the number of signatures and lambda to be considered
num_processes = 5
K = 5
my_seed_best_1 = 16599
my_seed_best_2 = 54200
lambda_values = c(0.10,0.20)

# ANALYSIS USING THE GERMLINE SIGNATURE ESTIMATED BY US AS BACKGROUND

# STEP 4: perform the signature discovery for the best configurations
load(file="paper_final_experiments/initial_betas_our_germline_nmf_standard_cv_10.RData")
final_solution_our_germline_nmf_standard_cv_10_lambda_10 = nmfLassoK(x=patients,K=K,beta=initial_betas_germline_cv_10[3,1][[1]],lambda_rate=lambda_values[1],iterations=100,max_iterations_lasso=100000,num_processes=num_processes,seed=my_seed_best_1,verbose=TRUE)
final_solution_our_germline_nmf_standard_cv_10_lambda_15 = nmfLassoK(x=patients,K=K,beta=initial_betas_germline_cv_10[3,1][[1]],lambda_rate=lambda_values[2],iterations=100,max_iterations_lasso=100000,num_processes=num_processes,seed=my_seed_best_2,verbose=TRUE)

# save the results
save(final_solution_our_germline_nmf_standard_cv_10_lambda_10,file="paper_final_experiments/final_solution_our_germline_nmf_standard_cv_10_lambda_10.RData")
save(final_solution_our_germline_nmf_standard_cv_10_lambda_15,file="paper_final_experiments/final_solution_our_germline_nmf_standard_cv_10_lambda_15.RData")
