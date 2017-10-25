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
load("data/finalbg2.RData")
background_germline = finalbg$freq/sum(finalbg$freq)

# set the number of signatures and lambda to be considered
num_processes = 48
K = 3:9
nmf_runs = 100
my_seed_starting_beta = 76509
cross_validation_entries = 0.10
cross_validation_iterations = 10
lambda_values = c(0.10,0.15)
num_iterations_cv = 100

# ANALYSIS USING THE GERMLINE SIGNATURE ESTIMATED BY THE PAPER AS BACKGROUND

# STEP 1: fit the initial betas for each configuration
initial_betas_germline_cv_10 = startingBetasEstimation(x=patients,K=K,background_signature= background_germline,nmf_method="nmf_standard",nmf_runs=nmf_runs,num_processes=num_processes,seed=my_seed_starting_beta,verbose=TRUE)
save(initial_betas_germline_cv_10,file="paper_final_experiments/initial_betas_paper_germline_nmf_standard_cv_10.RData")
