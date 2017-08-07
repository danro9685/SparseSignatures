#Set working directory
setwd("~/Documents/GitHub/SparseSignatures/")

#Load data
load("data/patients.RData")
load("data/genome.RData")

#Load libraries
library(NMF)

#Set seed
set.seed(111)

#Find signatures by NMF for varying values of K
for(K in 1:30){
  cat(paste0(K, "\n"))
  nmfsigs = basis(nmf(t(patients), rank = K))
  nmfsigs = t(nmfsigs)
  save(nmfsigs, file = paste0("input_nmf/nmfsigs_", K, ".RData"))
  cat(paste0("saved file for K = ", K, "\n"))
}
