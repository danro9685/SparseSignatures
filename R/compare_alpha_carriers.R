library(ggplot2)
library(data.table)

setwd("~/Documents/GitHub/SparseSignatures/")

load("data/signatures_nmfLasso_new_background.RData")
#load("data/signatures_nmfLasso_new_background_v2.RData")
#load("data/signatures_nmfLasso.RData")
load("data/patients.RData")
load("data/clinical.RData")
load("data/brca_status.RData")

#Choose best alpha from grid search
best_alpha = signatures_nmfLasso$best_configuration$alpha

#Normalize alpha to sum to 1 for each patient
alpha_norm = as.data.table(best_alpha/rowSums(best_alpha))
colnames(alpha_norm)[2:ncol(alpha_norm)] = paste0("S", 1:(ncol(alpha_norm)-1))
alpha_norm$patient = rownames(patients)

#Merge germline data
alpha_norm = merge(alpha_norm, brca_status[, c("Sample", "Gene")], by.x = "patient", by.y = "Sample")
alpha_norm[, germline := ifelse(Gene %in% c("BRCA1", "BRCA2"), "BRCA", "Control"), by = 1:nrow(alpha_norm)]

#Plot normalized alphas for carriers and controls
ggplot(melt(alpha_norm[, c(1:7, 9)], id.vars = c(1, 8), variable.name = "signature")) + 
  geom_boxplot(aes(x = signature, y = value, fill = germline)) 

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
