library(ggplot2)

setwd("~/Documents/GitHub/SparseSignatures/")

load("data/signatures_with_genome.RData")
load("data/patients.RData")
load("data/clinical_data.RData")

#Choose best alpha from grid search
best_alpha = signatures_with_genome$best_configuration$alpha

#Normalize alpha to sum to 1 for each patient
alpha_norm = as.data.frame(best_alpha/rowSums(best_alpha))
colnames(alpha_norm)[2:ncol(alpha_norm)] = paste0("S", 1:(ncol(alpha_norm)-1))
alpha_norm$patient = rownames(patients)

#Merge germline data
alpha_norm = merge(alpha_norm, clinical_data[, c("Sample", "Gene")], by.x = "patient", by.y = "Sample")
alpha_norm$germline = ifelse(alpha_norm$Gene %in% c("BRCA1", "BRCA2"), "BRCA", "Control")

#Plot normalized alphas for carriers and controls
ggplot(melt(alpha_norm[, c(1:7, 9)], id.vars = c(1, 8), variable.name = "signature")) + geom_boxplot(aes(x = signature, y = value, fill = germline))
