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
load("data/finalbg2.RData")
germline_signature = finalbg$freq/sum(finalbg$freq)

# option 1: our model for 5 signatures is correct
set.seed(43561)
load("data/signatures_nmfLasso_germline_15.RData")
alpha = signatures_nmfLasso_germline[["best_configuration"]][["alpha"]]
beta = signatures_nmfLasso_germline[["best_configuration"]][["beta"]]
predicted_counts = round(alpha%*%beta)
mse = sum((patients-predicted_counts)^2)/nrow(patients)
mse_loss = MSE(predicted_counts,as.matrix(patients))
print(mse) # 114677.9
print(mse_loss) # 1194.561

# option 2: our model for 7 signatures is correct
set.seed(51233)
load("data/signatures_nmfLasso_germline_10.RData")
alpha = signatures_nmfLasso_germline[["best_configuration"]][["alpha"]]
beta = signatures_nmfLasso_germline[["best_configuration"]][["beta"]]
predicted_counts = round(alpha%*%beta)
mse = sum((patients-predicted_counts)^2)/nrow(patients)
mse_loss = MSE(predicted_counts,as.matrix(patients))
print(mse) # 46141.05
print(mse_loss) # 480.6359

# check how much the loss drop by removing Stratton-like Signtures 3 and 8
predicted_counts = round(alpha[,-c(5,7)]%*%beta[-c(5,7),])
mse_reduced = sum((patients-predicted_counts)^2)/nrow(patients)
mse_loss_reduced = MSE(predicted_counts,as.matrix(patients))
print(mse_reduced) # 167634.5
print(mse_loss_reduced) # 1746.193
nmf_beta_2_v2 = basis(nmf(t(predicted_counts),rank=7,nrun=10))
nmf_beta_2_v2 = t(nmf_beta_2_v2)/colSums(nmf_beta_2_v2)
plotSignatures(nmf_beta_2_v2,backgroundFreq=FALSE)

# option 3: initial betas for 5 signatures is correct
set.seed(76800)
load("data/initial_betas_germline_15.RData")
beta = rbind(germline_signature,initial_betas_germline[[4,1]])
beta = beta/rowSums(beta)
alpha = array(0,c(560,6))
for(j in 1:560) {
    alpha[j,] = nnls(t(beta),as.vector(patients[j,]))$x
}
predicted_counts = round(alpha%*%beta)
mse = sum((patients-predicted_counts)^2)/nrow(patients)
mse_loss = MSE(predicted_counts,as.matrix(patients))
print(mse) # 225858.3
print(mse_loss) # 2352.69
nmf_beta_1 = basis(nmf(t(predicted_counts),rank=5,nrun=10))
nmf_beta_1 = t(nmf_beta_1)/colSums(nmf_beta_1)
plotSignatures(nmf_beta_1,backgroundFreq=FALSE)
nmf_beta_1_v2 = basis(nmf(t(predicted_counts),rank=7,nrun=10))
nmf_beta_1_v2 = t(nmf_beta_1_v2)/colSums(nmf_beta_1_v2)
plotSignatures(nmf_beta_1_v2,backgroundFreq=FALSE)

# option 4: initial betas for 7 signatures is correct
set.seed(12166)
load("data/initial_betas_germline_10.RData")
beta = rbind(germline_signature,initial_betas_germline[[6,1]])
beta = beta/rowSums(beta)
alpha = array(0,c(560,8))
for(j in 1:560) {
    alpha[j,] = nnls(t(beta),as.vector(patients[j,]))$x
}
predicted_counts = round(alpha%*%beta)
mse = sum((patients-predicted_counts)^2)/nrow(patients)
mse_loss = MSE(predicted_counts,as.matrix(patients))
print(mse) # 57880.83
print(mse_loss) # 602.9254
nmf_beta_2 = basis(nmf(t(predicted_counts),rank=7,nrun=10))
nmf_beta_2 = t(nmf_beta_2)/colSums(nmf_beta_2)
plotSignatures(nmf_beta_2,backgroundFreq=FALSE)

# Same analysis but using only Stratton's Signatures 3 and 8

# option 2: our model for 7 signatures is correct
set.seed(21344)
load("data/signatures_nmfLasso_germline_10.RData")
alpha = signatures_nmfLasso_germline[["best_configuration"]][["alpha"]]
beta = signatures_nmfLasso_germline[["best_configuration"]][["beta"]]
predicted_counts = round(alpha[,c(1,5,7)]%*%beta[c(1,5,7),])
mse = sum((patients-predicted_counts)^2)/nrow(patients)
mse_loss = MSE(predicted_counts,as.matrix(patients))
print(mse) # 7402204
print(mse_loss) # 77106.29

# option 4: initial betas for 7 signatures is correct
set.seed(96522)
load("data/initial_betas_germline_10.RData")
beta = rbind(germline_signature,initial_betas_germline[[6,1]])
beta = beta/rowSums(beta)
alpha = array(0,c(560,3))
for(j in 1:560) {
    alpha[j,] = nnls(t(beta[c(1,5,7)]),as.vector(patients[j,]))$x
}
predicted_counts = round(alpha%*%beta[c(1,5,7),])
mse = sum((patients-predicted_counts)^2)/nrow(patients)
mse_loss = MSE(predicted_counts,as.matrix(patients))
print(mse) # 7646175
print(mse_loss) # 79647.66

# Now let's directly use Stratton's signatures
cosmic_signatures_stratton = read.table(file="cosmic_signatures_stratton.txt",sep="\t",header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
rownames(cosmic_signatures_stratton) = cosmic_signatures_stratton[,3]
cosmic_signatures_stratton[,"Substitution Type"] = NULL
cosmic_signatures_stratton[,"Trinucleotide"] = NULL
cosmic_signatures_stratton[,"Somatic Mutation Type"] = NULL
cosmic_signatures_stratton = t(cosmic_signatures_stratton)

# get only the Signatures found in breast tumors
initial_betas = cosmic_signatures_stratton[c("Signature 1","Signature 2","Signature 3","Signature 5","Signature 6","Signature 8","Signature 13","Signature 17","Signature 18","Signature 20","Signature 26","Signature 30"),]
rownames(initial_betas) = c("Signature 1","Signature 2","Signature 3","Signature 5","Signature 6","Signature 8","Signature 13","Signature 17","Signature 18","Signature 20","Signature 26","Signature 30")
set.seed(43522)
res_Stratton_starts = nmfLassoDecomposition(x=patients,beta=rbind(germline_signature,initial_betas),lambda_rate=0.01)
plot(res_Stratton_starts$loglik_progression)
plotSignatures(res_Stratton_starts$beta,backgroundFreq=TRUE,signatures_names=c("Signature 1","Signature 2","Signature 3","Signature 5","Signature 6","Signature 8","Signature 13","Signature 17","Signature 18","Signature 20","Signature 26","Signature 30"))

# evaluate the mse for Stratton original signatures
beta_Stratton_original = rbind(germline_signature,initial_betas)
alpha_Stratton_original = array(0,c(560,13))
for(j in 1:560) {
    alpha_Stratton_original[j,] = nnls(t(beta_Stratton_original),as.vector(patients[j,]))$x
}
predicted_counts = round(alpha_Stratton_original%*%beta_Stratton_original)
mse = sum((patients-predicted_counts)^2)/nrow(patients)
mse_loss = MSE(predicted_counts,as.matrix(patients))
print(mse) # 104600.6
print(mse_loss) # 1089.59

# evaluate the mse for Stratton sparsified signatures
beta_Stratton_lasso = res_Stratton_starts$beta
alpha_Stratton_lasso = res_Stratton_starts$alpha
predicted_counts = round(alpha_Stratton_lasso%*%beta_Stratton_lasso)
mse = sum((patients-predicted_counts)^2)/nrow(patients)
mse_loss = MSE(predicted_counts,as.matrix(patients))
print(mse) # 89265.46
print(mse_loss) # 929.8486
