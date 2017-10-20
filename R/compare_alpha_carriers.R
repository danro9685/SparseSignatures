library(ggplot2)
library(data.table)
library(nnls)

setwd("~/Documents/GitHub/SparseSignatures/")

load("data/signatures_nmfLasso_germline_15.RData")
#load("data/signatures_nmfLasso_new_background_v2.RData")
load("data/patients.RData")
load("data/clinical.RData")
load("data/brca_status.RData")
load("data/genome.RData")

#Function to plot signatures
"plotSignatures" <- function( beta, patients_ids = colnames(beta), backgroundFreq = TRUE ) {
  
  # set rownames and colnames
  if(backgroundFreq) {
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

#Choose best alpha from grid search
best_alpha = as.data.table(signatures_nmfLasso_germline$best_configuration$alpha)
colnames(best_alpha)[2:ncol(best_alpha)] = paste0("S", 1:(ncol(best_alpha)-1))
best_alpha$patient = rownames(patients)

#Merge germline data
best_alpha = merge(best_alpha, brca_status[, c("Sample", "Gene")], by.x = "patient", by.y = "Sample")
best_alpha[, germline := ifelse(Gene %in% c("BRCA1", "BRCA2"), "BRCA+", "BRCA-WT"), by = 1:nrow(best_alpha)]
best_alpha = melt(best_alpha, id.vars = c(1, 8, 9), variable.name = "signature")

#Normalize alpha to sum to 1 for each patient
best_alpha[, norm:=value/sum(value), by = patient]

#Adjust alpha to equalize means for BRCA and non-BRCA
avg_mutations = best_alpha[, sum(value), by=.(patient, germline)][, median(V1), by = germline]
adjustment_factor = avg_mutations[germline=="BRCA+", V1]/avg_mutations[germline=="BRCA-WT", V1]
best_alpha[, adjusted:=ifelse(germline=="BRCA+", value/adjustment_factor, value)]

#Plot without normalization
plt1 = ggplot(best_alpha) + geom_boxplot(aes(x = signature, y = value, fill = Gene)) +ylim(c(0, 10000)) + ggtitle("Alpha values (raw)")
plt1
#Plot normalized alphas for carriers and controls
plt2= ggplot(best_alpha) +  geom_boxplot(aes(x = signature, y = norm, fill = germline)) + ggtitle("Alpha values (normalized for each patient)")
plt2
#Plot rescaled alpha
plt3 = ggplot(best_alpha) +  geom_boxplot(aes(x = signature, y = adjusted, fill = germline))+ylim(c(0, 7500)) + ggtitle("Alpha values (rescaled for BRCA+ and BRCA-WT)")
plt3
#####
#Construct fake beta
beta = as.data.table(signatures_nmfLasso_germline$best_configuration$beta)
peaks = c("T[C>A]A", "T[C>A]C", "T[C>A]G" , "T[C>A]T","T[C>G]A", "T[C>G]C", "T[C>G]G" , "T[C>G]T")
beta[5,which(colnames(beta) %in% peaks), with=F]

beta_1 = beta[5,]
beta_2 = beta[5,]

beta_1[,c("T[C>A]A", "T[C>A]C", "T[C>A]G" , "T[C>A]T","T[C>G]A", "T[C>G]C", "T[C>G]G" , "T[C>G]T")]=0
beta_1 = beta_1/sum(beta_1)
beta_2[,setdiff(colnames(beta_2), c("T[C>A]A", "T[C>A]C", "T[C>A]G" , "T[C>A]T","T[C>G]A", "T[C>G]C", "T[C>G]G" , "T[C>G]T"))]=0
beta_2 = beta_2/sum(beta_2)

beta_new = rbind(beta[1:4,], beta_1, beta_2, beta[6,])
plotSignatures(as.matrix(beta_new))

#Calculate new alphas
alpha_new = matrix(0, nrow=560, ncol=7)
for(j in 1:560) {
  alpha_new[j,] = nnls(t(beta_new),as.vector(patients[j,]))$x
}
alpha_new = as.data.table(alpha_new)
colnames(alpha_new) = c("background_signature", paste0("S", 1:(ncol(alpha_new)-1)))
alpha_new$patient = rownames(patients)

#Merge germline data
alpha_new = merge(alpha_new, brca_status[, c("Sample", "Gene")], by.x = "patient", by.y = "Sample")
alpha_new[, germline := ifelse(Gene %in% c("BRCA1", "BRCA2"), "BRCA+", "BRCA-WT"), by = 1:nrow(alpha_new)]
alpha_new = melt(alpha_new, id.vars = c(1, 9, 10), variable.name = "signature")

#Normalize alpha to sum to 1 for each patient
alpha_new[, norm:=value/sum(value), by = patient]

#Adjust alpha to equalize means for BRCA and non-BRCA
avg_mutations = alpha_new[, sum(value), by=.(patient, germline)][, mean(V1), by = germline]
adjustment_factor = avg_mutations[germline=="BRCA+", V1]/avg_mutations[germline=="BRCA-WT", V1]
alpha_new[, adjusted:=ifelse(germline=="BRCA+", value/adjustment_factor, value)]

#Plot without normalization
plt1 = ggplot(alpha_new) + geom_boxplot(aes(x = signature, y = value, fill = Gene)) +ylim(c(0, 10000)) + ggtitle("Alpha values (raw)")
plt1
#Plot normalized alphas for carriers and controls
plt2= ggplot(alpha_new) +  geom_boxplot(aes(x = signature, y = norm, fill = germline)) + ggtitle("Alpha values (normalized for each patient)")
plt2
#Plot rescaled alpha
plt3 = ggplot(alpha_new) +  geom_boxplot(aes(x = signature, y = adjusted, fill = germline))+ylim(c(0, 7500)) + ggtitle("Alpha values (rescaled for BRCA+ and BRCA-WT)")
plt3

#
beta_genome = beta_new
for(i in 1:nrow(beta_genome)){
  beta_genome[i,] = beta_genome[i,]/genome$freq
}
plotSignatures(as.matrix(beta_genome))

#Plot with 7 signatures
load("data/signatures_nmfLasso_germline_10.RData")
#Choose best alpha from grid search
alpha10 = as.data.table(signatures_nmfLasso_germline$best_configuration$alpha)
colnames(alpha10)[2:ncol(alpha10)] = paste0("S", 1:(ncol(alpha10)-1))
alpha10$patient = rownames(patients)

#Merge germline data
alpha10 = merge(alpha10, brca_status[, c("Sample", "Gene")], by.x = "patient", by.y = "Sample")
alpha10[, germline := ifelse(Gene %in% c("BRCA1", "BRCA2"), "BRCA+", "BRCA-WT"), by = 1:nrow(alpha10)]
alpha10 = melt(alpha10, id.vars = c(1,10, 11), variable.name = "signature")

#Normalize alpha to sum to 1 for each patient
alpha10[, norm:=value/sum(value), by = patient]

#Adjust alpha to equalize means for BRCA and non-BRCA
avg_mutations = alpha10[, sum(value), by=.(patient, germline)][, median(V1), by = germline]
adjustment_factor = avg_mutations[germline=="BRCA+", V1]/avg_mutations[germline=="BRCA-WT", V1]
best_alpha[, adjusted:=ifelse(germline=="BRCA+", value/adjustment_factor, value)]

#Plot without normalization
plt1 = ggplot(alpha10) + geom_boxplot(aes(x = signature, y = value, fill = Gene)) +ylim(c(0, 10000)) + ggtitle("Alpha values (raw)")
plt1
#Plot normalized alphas for carriers and controls
plt2= ggplot(best_alpha) +  geom_boxplot(aes(x = signature, y = norm, fill = germline)) + ggtitle("Alpha values (normalized for each patient)")
plt2
#Plot rescaled alpha
plt3 = ggplot(best_alpha) +  geom_boxplot(aes(x = signature, y = adjusted, fill = germline))+ylim(c(0, 7500)) + ggtitle("Alpha values (rescaled for BRCA+ and BRCA-WT)")
plt3








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
