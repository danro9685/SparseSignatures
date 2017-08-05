#Set working directory
setwd("~/Documents/GitHub/SparseSignatures/")

#Load libraru
require("deconstructSigs")

#Load data
load("data/cna.RData")
load("data/ssm.RData")

#Function to get snvs in specific intervals
intersectSnvCna = function(ssm, cna){
  snvs_in = c()
  for(i in 1:nrow(cna)){
    snvs_in = c(snvs_in, which(ssm[,1] == cna[i,1] & 
                    ssm[,2] == cna[i,2] & 
                    ssm[,3] >= cna[i,3] & 
                    ssm[,3] <= cna[i,4]))
  }
  return(snvs_in)
}

#Get snvs outside copy number altered regions
filterDiploidSnvs = function(snvs, cna){
  which_in = intersectSnvCna(snvs, cna)
  return(snvs[-which_in,])
}

#Filter Stratton calls
ssm_filtered = filterDiploidSnvs(ssm, cna)

#Make filtered motif table
patients_filtered = mut.to.sigs.input(snvs_out, "sample", "chrom", "pos", "ref", "alt", bsg = BSgenome.Hsapiens.1000genomes.hs37d5)
patients_filtered = patients_filtered[sort(colnames(patients_filtered))]

#Save data
save(patients_filtered, file = "patients_filtered.RData")
