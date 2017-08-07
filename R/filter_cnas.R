#Set working directory
setwd("~/Documents/GitHub/SparseSignatures/")

#Load libraru
require("deconstructSigs")
require("BSgenome.Hsapiens.1000genomes.hs37d5")
require("data.table")

#Load data
load("data/cna.RData")
load("data/ssm.RData")

setkey(cna, patient, chrom)
setkey(ssm, patient, chrom)

#Function to get snvs in specific intervals
intersectSnvCna = function(ssm, cna){
  snvs_in = c()
  ssm[, index:= 1:nrow(ssm)]
  for(i in 1:nrow(cna)){
    snvs_in = c(snvs_in, ssm[.(cna[i, patient], cna[i, chrom])][pos >= cna[i,start] & pos <= cna[i,end], index])
    if(i %% 5000 == 0){
     cat(paste0("completed ", i, " of ", nrow(cna), " rows.\n"))
    }
  }
  ssm[, index := NULL]
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
patients_filtered = mut.to.sigs.input(ssm_filtered, "patient", "chrom", "pos", "ref", "alt", bsg = BSgenome.Hsapiens.1000genomes.hs37d5)
patients_filtered = patients_filtered[sort(colnames(patients_filtered))]

#Save data
save(patients_filtered, file = "data/patients_filtered.RData")
