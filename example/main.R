# set the working directory
work.dir = "/Users/daniele/Documents/Stanford-github/SparseSignatures"
setwd(work.dir)

# load the required libraries and sources
library("SomaticSignatures")
library("nnlasso")
library("nnls")
source("R/mutation_factor.R")

# libraries required to plot the extracted signatures
library("data.table")
library("ggplot2")
library("gridExtra")

# function to plot extracted signatures
"plotSignatures" <- function( beta, patients_ids = colnames(beta), genomeFreq = TRUE ) {
  
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
load("data/clinical_data.RData")
load("data/patients.RData")
load("data/genome.RData")

# set the number of signatures and lambda to be considered
K = 2:15
lambda_values = c(0.01,0.05,0.1,0.2)
cross_validation_entries = 0.15

# fit the signatures with the genome frequencies as noise model
signatures_with_genome = nmfLasso(x=patients,K=K,background_signature=genome$freq,lambda_values=lambda_values,cross_validation_entries= cross_validation_entries,iterations=20,seed=59040,verbose=TRUE)
save(signatures_with_genome,file="data/signatures_with_genome.RData")

# fit the signatures without the genome frequencies as noise model
signatures_without_genome = nmfLasso(x=patients,K=K,background_signature=NULL,lambda_values=lambda_values,cross_validation_entries= cross_validation_entries,iterations=20,seed=59040,verbose=TRUE)
save(signatures_without_genome,file="data/signatures_without_genome.RData")

# plot the signatures
plotSignatures(signatures_with_genome$grid_search[[10,2]]$beta,patientss_ids=colnames(patients),genomeFreq=TRUE)
plotSignatures(signatures_without_genome$grid_search[[14,5]]$beta,patientss_ids=colnames(patients),genomeFreq=TRUE)

# plot the log-likelihood values
plot(signatures_with_genome$grid_search[[10,2]]$loglik)
plot(signatures_without_genome$grid_search[[14,5]]$loglik)
