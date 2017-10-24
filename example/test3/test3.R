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
load("data/allbg.RData")
germline_background = allbg$freq
load("data/brca_status.RData")
patients = patients[brca_status$Sample[which(brca_status$Gene%in%c("BRCA1","BRCA2"))],]

# perform the inference only on carriers
set.seed(456711)
starting_beta_3 = basis(nmf(t(patients),rank=3,nrun=100))
set.seed(12455)
starting_beta_4 = basis(nmf(t(patients),rank=4,nrun=100))
set.seed(32499)
starting_beta_5 = basis(nmf(t(patients),rank=5,nrun=100))
set.seed(12100)
starting_beta_6 = basis(nmf(t(patients),rank=6,nrun=100))
set.seed(56211)
starting_beta_7 = basis(nmf(t(patients),rank=7,nrun=100))

# plot the starting beta
plotSignatures(t(starting_beta_3)/rowSums(t(starting_beta_3)),backgroundFreq=FALSE)
plotSignatures(t(starting_beta_4)/rowSums(t(starting_beta_4)),backgroundFreq=FALSE)
plotSignatures(t(starting_beta_5)/rowSums(t(starting_beta_5)),backgroundFreq=FALSE)
plotSignatures(t(starting_beta_6)/rowSums(t(starting_beta_6)),backgroundFreq=FALSE)
plotSignatures(t(starting_beta_7)/rowSums(t(starting_beta_7)),backgroundFreq=FALSE)

# run our method on these inputs
set.seed(32311)
beta_3 = nmfLassoDecomposition(x=patients,beta=rbind(germline_background,t(starting_beta_3)),lambda_rate=0.15,verbose=TRUE)
set.seed(65600)
beta_4 = nmfLassoDecomposition(x=patients,beta=rbind(germline_background,t(starting_beta_4)),lambda_rate=0.15,verbose=TRUE)
set.seed(32455)
beta_5 = nmfLassoDecomposition(x=patients,beta=rbind(germline_background,t(starting_beta_5)),lambda_rate=0.15,verbose=TRUE)
set.seed(12136)
beta_6 = nmfLassoDecomposition(x=patients,beta=rbind(germline_background,t(starting_beta_6)),lambda_rate=0.15,verbose=TRUE)
set.seed(75788)
beta_7 = nmfLassoDecomposition(x=patients,beta=rbind(germline_background,t(starting_beta_7)),lambda_rate=0.15,verbose=TRUE)

# plot the final beta
plotSignatures(beta_3$beta,backgroundFreq=TRUE)
plotSignatures(beta_4$beta,backgroundFreq=TRUE)
plotSignatures(beta_5$beta,backgroundFreq=TRUE)
plotSignatures(beta_6$beta,backgroundFreq=TRUE)
plotSignatures(beta_7$beta,backgroundFreq=TRUE)

# compute mse for final_beta
mse_3 = sum((patients-round(beta_3$alpha%*%beta_3$beta))^2)/nrow(patients)
mse_4 = sum((patients-round(beta_4$alpha%*%beta_4$beta))^2)/nrow(patients)
mse_5 = sum((patients-round(beta_5$alpha%*%beta_5$beta))^2)/nrow(patients)
mse_6 = sum((patients-round(beta_6$alpha%*%beta_6$beta))^2)/nrow(patients)
mse_7 = sum((patients-round(beta_7$alpha%*%beta_7$beta))^2)/nrow(patients)
plot(c(mse_3,mse_4,mse_5,mse_6,mse_7))
