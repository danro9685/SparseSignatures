#Set working directory
setwd("~/Documents/GitHub/SparseSignatures/")

#Load data
<<<<<<< HEAD
load("../data/patients.RData")
load("../data/genome.RData")
=======
load("data/patients.RData")
load("data/genome.RData")
>>>>>>> f3c20e06bdc5e1ffaa0431adf6f5ad5e43969ec6

#Load libraries
library(NMF)

#Set seed
set.seed(111)

#Find signatures by NMF for varying values of K
for(K in 1:30){
<<<<<<< HEAD
  cat(paste0(K, "\n"))
  nmfsigs = basis(nmf(t(patients), rank = K))
  nmfsigs = nmfsigs[, order(rowMax(t(nmfsigs)), decreasing = TRUE)]
=======
  nmfsigs = basis(nmf(t(patients), rank = K))
>>>>>>> f3c20e06bdc5e1ffaa0431adf6f5ad5e43969ec6
  nmfsigs = t(nmfsigs)
  save(nmfsigs, file = paste0("input_nmf/nmfsigs_", K, ".RData"))
  cat(paste0("saved file for K = ", K, "\n"))
}
