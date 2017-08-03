#Load data
load("../data/patient.RData")
load("../data/genome.RData")

#Load libraries
library(NMF)
library(nnls)

#Set seed
set.seed(111)

#Find signatures by NMF for varying values of K
for(K in 1:30){
  nmfsigs = basis(nmf(t(patient), rank = K))
  nmfsigs = nmfsigs[, order(rowMax(t(nmfsigs)), decreasing = TRUE)]
  nmfsigs = t(nmfsigs)
  save(nmfsigs, file = paste0("../input_nmf/nmfsigs_", K, ".RData"))
}
