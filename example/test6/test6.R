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

# library for comparison metrics
library("MLmetrics")

# function to plot extracted signatures
"plotSignatures" <- function( beta, patients_ids = colnames(beta), backgroundFreq = TRUE, signatures_names = NULL ) {
    
    # set rownames and colnames
    if(backgroundFreq) {
        if(is.null(signatures_names)) {
            rownames(beta) = c("Background",paste0("S",1:(nrow(beta)-1)))
        }
        else {
            rownames(beta) = c("Background", signatures_names)
        }
    }
    else {
        if(is.null(signatures_names)) {
        rownames(beta) = paste0("S",1:nrow(beta))
        }
        else {
            rownames(beta) = signatures_names
        }
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
load("data/brca_status.RData")
load("data/allbg.RData")
genome_background = allbg$freq

# get the cohort of interest
BRCA_negative_patients = brca_status$Sample[!brca_status$Gene%in%c("BRCA1","BRCA2")]
patients = patients[BRCA_negative_patients,]

load("data/initial_betas_germline_cv_10_control_only.RData")

# consider all the values for initial_betas_germline_cv_10
set.seed(58375)
nmf_standard_results = list()
nmf_lasso_results = list()
cat(0,"\n")
for(i in 1:nrow(initial_betas_germline_cv_10)) {
    # get the starting beta for the current configuration
    curr_starting_beta = initial_betas_germline_cv_10[[i,1]]
    # NMF standard
    curr_alpha_nmf_standard = array(0,c(469,(i+3)))
    curr_beta_nmf_standard = rbind(genome_background,curr_starting_beta/rowSums(curr_starting_beta))
    for(j in 1:469) {
        curr_alpha_nmf_standard[j,] = nnls(t(curr_beta_nmf_standard),as.vector(patients[j,]))$x
    }
    curr_predicted_counts = round(curr_alpha_nmf_standard%*%curr_beta_nmf_standard)
    curr_mse = sum((patients-curr_predicted_counts)^2)/nrow(patients)
    curr_mse_loss = MSE(curr_predicted_counts,as.matrix(patients))
    curr_results_nmf_standard = list(alpha=curr_alpha_nmf_standard,beta=curr_beta_nmf_standard,mse=curr_mse,mse_loss=curr_mse_loss)
    # save the results for the current configuration
    nmf_standard_results[[i]] = curr_results_nmf_standard
    cat(i/nrow(initial_betas_germline_cv_10),"\n")
}

# NMF lasso
i = 1
set.seed(65788)
curr_starting_beta = initial_betas_germline_cv_10[[i,1]]
curr_res = nmfLassoDecomposition(x=patients,beta=rbind(genome_background,curr_starting_beta),lambda_rate=0.15,verbose=FALSE)
curr_alpha_nmf_lasso = curr_res$alpha
curr_beta_nmf_lasso = curr_res$beta
curr_predicted_counts = round(curr_alpha_nmf_lasso%*%curr_beta_nmf_lasso)
curr_mse = sum((patients-curr_predicted_counts)^2)/nrow(patients)
curr_mse_loss = MSE(curr_predicted_counts,as.matrix(patients))
curr_results_nmf_lasso = list(alpha=curr_alpha_nmf_lasso,beta=curr_beta_nmf_lasso,mse=curr_mse,mse_loss=curr_mse_loss)
nmf_lasso_results[[i]] = curr_results_nmf_lasso

i = 2
set.seed(96700)
curr_starting_beta = initial_betas_germline_cv_10[[i,1]]
curr_res = nmfLassoDecomposition(x=patients,beta=rbind(genome_background,curr_starting_beta),lambda_rate=0.15,verbose=FALSE)
curr_alpha_nmf_lasso = curr_res$alpha
curr_beta_nmf_lasso = curr_res$beta
curr_predicted_counts = round(curr_alpha_nmf_lasso%*%curr_beta_nmf_lasso)
curr_mse = sum((patients-curr_predicted_counts)^2)/nrow(patients)
curr_mse_loss = MSE(curr_predicted_counts,as.matrix(patients))
curr_results_nmf_lasso = list(alpha=curr_alpha_nmf_lasso,beta=curr_beta_nmf_lasso,mse=curr_mse,mse_loss=curr_mse_loss)
nmf_lasso_results[[i]] = curr_results_nmf_lasso

i = 3
set.seed(54677)
curr_starting_beta = initial_betas_germline_cv_10[[i,1]]
curr_res = nmfLassoDecomposition(x=patients,beta=rbind(genome_background,curr_starting_beta),lambda_rate=0.15,verbose=FALSE)
curr_alpha_nmf_lasso = curr_res$alpha
curr_beta_nmf_lasso = curr_res$beta
curr_predicted_counts = round(curr_alpha_nmf_lasso%*%curr_beta_nmf_lasso)
curr_mse = sum((patients-curr_predicted_counts)^2)/nrow(patients)
curr_mse_loss = MSE(curr_predicted_counts,as.matrix(patients))
curr_results_nmf_lasso = list(alpha=curr_alpha_nmf_lasso,beta=curr_beta_nmf_lasso,mse=curr_mse,mse_loss=curr_mse_loss)
nmf_lasso_results[[i]] = curr_results_nmf_lasso

i = 4
set.seed(12344)
curr_starting_beta = initial_betas_germline_cv_10[[i,1]]
curr_res = nmfLassoDecomposition(x=patients,beta=rbind(genome_background,curr_starting_beta),lambda_rate=0.15,verbose=FALSE)
curr_alpha_nmf_lasso = curr_res$alpha
curr_beta_nmf_lasso = curr_res$beta
curr_predicted_counts = round(curr_alpha_nmf_lasso%*%curr_beta_nmf_lasso)
curr_mse = sum((patients-curr_predicted_counts)^2)/nrow(patients)
curr_mse_loss = MSE(curr_predicted_counts,as.matrix(patients))
curr_results_nmf_lasso = list(alpha=curr_alpha_nmf_lasso,beta=curr_beta_nmf_lasso,mse=curr_mse,mse_loss=curr_mse_loss)
nmf_lasso_results[[i]] = curr_results_nmf_lasso

i = 5
set.seed(43566)
curr_starting_beta = initial_betas_germline_cv_10[[i,1]]
curr_res = nmfLassoDecomposition(x=patients,beta=rbind(genome_background,curr_starting_beta),lambda_rate=0.15,verbose=FALSE)
curr_alpha_nmf_lasso = curr_res$alpha
curr_beta_nmf_lasso = curr_res$beta
curr_predicted_counts = round(curr_alpha_nmf_lasso%*%curr_beta_nmf_lasso)
curr_mse = sum((patients-curr_predicted_counts)^2)/nrow(patients)
curr_mse_loss = MSE(curr_predicted_counts,as.matrix(patients))
curr_results_nmf_lasso = list(alpha=curr_alpha_nmf_lasso,beta=curr_beta_nmf_lasso,mse=curr_mse,mse_loss=curr_mse_loss)
nmf_lasso_results[[i]] = curr_results_nmf_lasso

i = 6
set.seed(96577)
curr_starting_beta = initial_betas_germline_cv_10[[i,1]]
curr_res = nmfLassoDecomposition(x=patients,beta=rbind(genome_background,curr_starting_beta),lambda_rate=0.15,verbose=FALSE)
curr_alpha_nmf_lasso = curr_res$alpha
curr_beta_nmf_lasso = curr_res$beta
curr_predicted_counts = round(curr_alpha_nmf_lasso%*%curr_beta_nmf_lasso)
curr_mse = sum((patients-curr_predicted_counts)^2)/nrow(patients)
curr_mse_loss = MSE(curr_predicted_counts,as.matrix(patients))
curr_results_nmf_lasso = list(alpha=curr_alpha_nmf_lasso,beta=curr_beta_nmf_lasso,mse=curr_mse,mse_loss=curr_mse_loss)
nmf_lasso_results[[i]] = curr_results_nmf_lasso

i = 7
set.seed(32144)
curr_starting_beta = initial_betas_germline_cv_10[[i,1]]
curr_res = nmfLassoDecomposition(x=patients,beta=rbind(genome_background,curr_starting_beta),lambda_rate=0.15,verbose=FALSE)
curr_alpha_nmf_lasso = curr_res$alpha
curr_beta_nmf_lasso = curr_res$beta
curr_predicted_counts = round(curr_alpha_nmf_lasso%*%curr_beta_nmf_lasso)
curr_mse = sum((patients-curr_predicted_counts)^2)/nrow(patients)
curr_mse_loss = MSE(curr_predicted_counts,as.matrix(patients))
curr_results_nmf_lasso = list(alpha=curr_alpha_nmf_lasso,beta=curr_beta_nmf_lasso,mse=curr_mse,mse_loss=curr_mse_loss)
nmf_lasso_results[[i]] = curr_results_nmf_lasso

nmf_results = list(standard=nmf_standard_results,lasso=nmf_lasso_results)

# save the results
save(nmf_results,file="example/test6/nmf_results.RData")

# compute the mse progression
mse_standard = NULL
mse_lasso = NULL
for(i in 1:7) {
    mse_standard = c(mse_standard,nmf_results$standard[[i]]$mse)
    mse_lasso = c(mse_lasso,nmf_results$lasso[[i]]$mse)
}

# plots standard
plotSignatures(nmf_results$standard[[1]]$beta)
plotSignatures(nmf_results$standard[[2]]$beta)
plotSignatures(nmf_results$standard[[3]]$beta)
plotSignatures(nmf_results$standard[[4]]$beta)
plotSignatures(nmf_results$standard[[5]]$beta)
plotSignatures(nmf_results$standard[[6]]$beta)
plotSignatures(nmf_results$standard[[7]]$beta)

# plots lasso
plotSignatures(nmf_results$lasso[[1]]$beta)
plotSignatures(nmf_results$lasso[[2]]$beta)
plotSignatures(nmf_results$lasso[[3]]$beta)
plotSignatures(nmf_results$lasso[[4]]$beta)
plotSignatures(nmf_results$lasso[[5]]$beta)
plotSignatures(nmf_results$lasso[[6]]$beta)
plotSignatures(nmf_results$lasso[[7]]$beta)

# make the final inference for the best configurations
i = 3
set.seed(43214)
curr_starting_beta = initial_betas_germline_cv_10[[i,1]]
curr_res = nmfLassoDecomposition(x=patients,beta=rbind(genome_background,curr_starting_beta),lambda_rate=0.25,iterations=100,verbose=TRUE)
curr_alpha_nmf_lasso = curr_res$alpha
curr_beta_nmf_lasso = curr_res$beta
curr_predicted_counts = round(curr_alpha_nmf_lasso%*%curr_beta_nmf_lasso)
curr_mse = sum((patients-curr_predicted_counts)^2)/nrow(patients) # 39570.87
curr_mse_loss = MSE(curr_predicted_counts,as.matrix(patients)) # 412.1965
plot(curr_res$loglik_progression)
plotSignatures(curr_res$beta)
