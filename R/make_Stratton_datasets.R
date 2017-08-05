setwd("~/Documents/GitHub/SparseSignatures/")

library(data.table)


cna = fread("/Volumes/Seagate Bac/BRCA_genomics/additional_data/copy_number_somatic_mutation.tsv")
ssm = fread("/Volumes/Seagate Bac/BRCA_genomics/additional_data/simple_somatic_mutation_controlled_reduced_2.txt")

ssm = unique(ssm)

ssm = ssm[, c(1:3, 5:6)]
colnames(ssm) = c("patient", "chrom", "pos", "ref", "alt")
save(ssm, file = "data/ssm.RData")

cna = cna[, c(6, 12, 13, 14)]
cna = unique(cna)
colnames(cna) = c("patient", "chrom", "start", "end")
save(cna, file = "data/cna.RData")

require("deconstructSigs")

#Function to get snvs in specific intervals
intersectSnvCna = function(snvs, cna){
  snvs_in = c()
  for(i in 1:nrow(cna)){
    snvs_in = which(snvs[,1] == cna[i,1] & snvs[,2] == cna[i,2] & snvs[,3] >= cna[i,3] & snvs[,3] <= cna[i,4])
  }
  return(snvs_in)
}

#Get snvs outside copy number altered regions
filterDiploidSnvs = function(snvs, cna){
  which_in = intersectSnvCna(snvs, cna)
  return(snvs[-which_in,])
}