library(ggplot2)
library(data.table)
library(nnls)

setwd("~/Documents/GitHub/SparseSignatures/")

load("data/signatures_nmfLasso_germline_15.RData")
#load("data/signatures_nmfLasso_new_background_v2.RData")
load("data/patients.RData")
load("data/clinical.RData")
load("data/brca_status.RData")

#Choose best alpha from grid search
best_alpha = as.data.table(signatures_nmfLasso_germline$best_configuration$alpha)
colnames(best_alpha)[2:ncol(best_alpha)] = paste0("S", 1:(ncol(best_alpha)-1))
best_alpha$patient = rownames(patients)

#Merge germline data
best_alpha = merge(best_alpha, brca_status[, c("Sample", "Gene")], by.x = "patient", by.y = "Sample")
best_alpha[, germline := ifelse(Gene %in% c("BRCA1", "BRCA2"), "BRCA+", "BRCA-WT"), by = 1:nrow(best_alpha)]
best_alpha = melt(best_alpha[, c(1:7, 9)], id.vars = c(1, 8), variable.name = "signature")

#Normalize alpha to sum to 1 for each patient
best_alpha[, norm:=value/sum(value), by = patient]

#Adjust alpha to equalize means for BRCA and non-BRCA
avg_mutations = best_alpha[, sum(value), by=.(patient, germline)][, mean(V1), by = germline]
adjustment_factor = avg_mutations[germline=="BRCA+", V1]/avg_mutations[germline=="BRCA-WT", V1]
best_alpha[, adjusted:=ifelse(germline=="BRCA+", value/adjustment_factor, value)]

#Plot without normalization
plt1 = ggplot(best_alpha) + geom_boxplot(aes(x = signature, y = value, fill = germline)) +ylim(c(0, 10000)) + ggtitle("Alpha values (raw)")

#Plot normalized alphas for carriers and controls
plt2= ggplot(best_alpha) +  geom_boxplot(aes(x = signature, y = norm, fill = germline)) + ggtitle("Alpha values (normalized for each patient)")

#Plot rescaled alpha
plt3 = ggplot(best_alpha) +  geom_boxplot(aes(x = signature, y = adjusted, fill = germline))+ylim(c(0, 7500)) + ggtitle("Alpha values (rescaled for BRCA+ and BRCA-WT)")

#####
#Construct fake beta
beta = as.data.table(signatures_nmfLasso_germline$best_configuration$beta)
peaks = c("T[C>A]A", "T[C>A]C", "T[C>A]G" , "T[C>A]T","T[C>G]A", "T[C>G]C", "T[C>G]G" , "T[C>G]T")
beta[5,which(colnames(beta) %in% peaks), with=F]

beta_1 = beta[5,]
beta_2 = beta[5,]

beta_1[,c("T[C>A]A", "T[C>A]C", "T[C>A]G" , "T[C>A]T","T[C>G]A", "T[C>G]C", "T[C>G]G" , "T[C>G]T")]=0
beta_2[,setdiff(colnames(beta_2), c("T[C>A]A", "T[C>A]C", "T[C>A]G" , "T[C>A]T","T[C>G]A", "T[C>G]C", "T[C>G]G" , "T[C>G]T"))]=0

beta_new = rbind(beta[1:4,], beta_1, beta_2, beta[6,])

#Calculate new alphas
alpha_new = matrix(0, nrow=560, ncol=7)
for(j in 1:560) {
  alpha_new[j,] = nnls(t(beta_new),as.vector(patients[j,]))$x
}



#Get p-values
for(sig in colnames(alpha_norm)[2:(ncol(alpha_norm)-2)]){
  print(sig)
  print(wilcox.test(alpha_norm[germline=="BRCA",get(sig)], alpha_norm[germline!="BRCA",get(sig)]))
}

#Include triple-negatives for comparison
clinical_merge = clinical[,.(sample_name, final.ER, final.PR, final.HER2)]
clinical_merge[, patient := paste0(sample_name, "a")]
clinical_merge[, TN := ifelse(final.ER=="negative" & final.PR=="negative" & final.HER2=="negative", "TN", "other"), by = 1:nrow(clinical_merge)]
alpha_norm = merge(alpha_norm, clinical_merge[, c("patient", "TN")], by = "patient")
alpha_norm[, status := "other"]
alpha_norm[TN=="TN", status := "TN"]
alpha_norm[germline=="BRCA", status := "BRCA"]
alpha_norm[, status := factor(status, levels = c("BRCA", "TN", "other"))]
ggplot(melt(alpha_norm[, c(1:7, 11)], id.vars = c(1, 8), variable.name = "signature")) + geom_boxplot(aes(x = signature, y = value, fill = status))
