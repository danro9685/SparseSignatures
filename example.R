# set the working directory
work.dir = "/Users/daniele/Documents/Stanford-github/SparseSignatures"
setwd(work.dir)

# load the required libraries and sources
library("SomaticSignatures")
library("nnlasso")
library("nnls")
source("R/mutation_factor.R")
library("data.table")
library("ggplot2")
library("gridExtra")

# function to plot extracted signatures
"plotSignatures" <- function( beta, patients_ids, genomeFreq = TRUE ) {
    
    # set rownames and colnames
    if(genomeFreq) {
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
        facet_wrap(~alt, nrow=1, scales = "free_x") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        ggtitle(sig) + 
        theme(legend.position="none")
    } 
    grid.arrange(grobs=glist,ncol=ceiling(nrow(beta)/6))
}

# read the data
clinical_data = read.csv("data/clinical_data.csv")
load("data/patient.RData")
load("data/genome.RData")

# set the number of signatures
K = 7

# fit the signatures with the genome frequencies as noise model
signatures_with_genome = nmfLasso(x=patient,K=K,background_signature=genome$freq,iterations=20,lambda_rate=0.01,seed=43899,verbose=TRUE)

# fit the signatures without the genome frequencies as noise model
signatures_without_genome = nmfLasso(x=patient,K=K,background_signature=NULL,iterations=20,lambda_rate=0.01,seed=43899,verbose=TRUE)

# plot the signatures
plotSignatures(signatures_with_genome$beta,patients_ids=colnames(patient),genomeFreq=TRUE)
plotSignatures(signatures_without_genome$beta,patients_ids=colnames(patient),genomeFreq=TRUE)

# plot the log-likelihood values
plot(signatures_with_genome$loglik)
plot(signatures_without_genome$loglik)
