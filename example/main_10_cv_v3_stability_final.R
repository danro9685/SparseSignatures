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

# perform the analysis with random starting points

# K = 4
set.seed(56477)
curr_starting_beta_4 = initial_betas_germline_cv_10["4_signatures",][[1]]
curr_starting_beta_4[1,] = runif(96)
curr_starting_beta_4[2,] = runif(96)
curr_starting_beta_4[3,] = runif(96)
curr_starting_beta_4[4,] = runif(96)
curr_res_4 = nmfLassoDecomposition(x=patients,beta=rbind(background_genome,curr_starting_beta_4),lambda_rate=0.05,iterations=100,verbose=TRUE)
curr_alpha_nmf_lasso_4 = curr_res_4$alpha
curr_beta_nmf_lasso_4 = curr_res_4$beta
curr_predicted_counts_4 = round(curr_alpha_nmf_lasso_4%*%curr_beta_nmf_lasso_4)
curr_mse_4 = sum((patients-curr_predicted_counts_4)^2)/nrow(patients) # 92673.87
plot(curr_res_4$loglik_progression)
plotSignatures(curr_res_4$beta)

# K = 5
set.seed(75833)
curr_starting_beta_5 = initial_betas_germline_cv_10["5_signatures",][[1]]
curr_starting_beta_5[1,] = runif(96)
curr_starting_beta_5[2,] = runif(96)
curr_starting_beta_5[3,] = runif(96)
curr_starting_beta_5[4,] = runif(96)
curr_starting_beta_5[5,] = runif(96)
curr_res_5 = nmfLassoDecomposition(x=patients,beta=rbind(background_genome,curr_starting_beta_5),lambda_rate=0.05,iterations=100,verbose=TRUE)
curr_alpha_nmf_lasso_5 = curr_res_5$alpha
curr_beta_nmf_lasso_5 = curr_res_5$beta
curr_predicted_counts_5 = round(curr_alpha_nmf_lasso_5%*%curr_beta_nmf_lasso_5)
curr_mse_5 = sum((patients-curr_predicted_counts_5)^2)/nrow(patients) # 92673.87
plot(curr_res_5$loglik_progression)
plotSignatures(curr_res_5$beta)

# K = 6
set.seed(36411)
curr_starting_beta_6 = initial_betas_germline_cv_10["6_signatures",][[1]]
curr_starting_beta_6[1,] = runif(96)
curr_starting_beta_6[2,] = runif(96)
curr_starting_beta_6[3,] = runif(96)
curr_starting_beta_6[4,] = runif(96)
curr_starting_beta_6[5,] = runif(96)
curr_starting_beta_6[6,] = runif(96)
curr_res_6 = nmfLassoDecomposition(x=patients,beta=rbind(background_genome,curr_starting_beta_6),lambda_rate=0.05,iterations=100,verbose=TRUE)
curr_alpha_nmf_lasso_6 = curr_res_6$alpha
curr_beta_nmf_lasso_6 = curr_res_6$beta
curr_predicted_counts_6 = round(curr_alpha_nmf_lasso_6%*%curr_beta_nmf_lasso_6)
curr_mse_6 = sum((patients-curr_predicted_counts_6)^2)/nrow(patients) # 44471.28
plot(curr_res_6$loglik_progression)
plotSignatures(curr_res_6$beta)
