setwd("/Users/daniele/Documents/Stanford-github/SparseSignatures/")

library(ggplot2)
library(gridExtra)
library(data.table)

load("/Users/daniele/Documents/Stanford-github/SparseSignatures/paper_final_experiments/final_solution_our_germline_nmf_standard_cv_10_lambda_15.RData")
signatures = final_solution_our_germline_nmf_standard_cv_10_lambda_15
load("/Users/daniele/Documents/Stanford-github/SparseSignatures/data/patients.RData")
load("/Users/daniele/Documents/Stanford-github/SparseSignatures/data/clinical.RData")
load("/Users/daniele/Documents/Stanford-github/SparseSignatures/data/brca_status.RData")
load("/Users/daniele/Documents/Stanford-github/SparseSignatures/data/patients.RData")

# get the best alpha
best_alpha = as.data.table(signatures$alpha)
colnames(best_alpha)[2:ncol(best_alpha)] = paste0("S", 1:(ncol(best_alpha)-1))
best_alpha$patient = rownames(patients)

# merge Subtype data
best_alpha = merge(best_alpha, brca_status[, c("Sample", "Gene")], by.x = "patient", by.y = "Sample")
best_alpha[, Subtype := ifelse(Gene %in% c("BRCA1", "BRCA2"), "BRCA+", "BRCA-WT"), by = 1:nrow(best_alpha)]

# get Subtypes

# Subtype 1: BRCA+ (91 patients)
BRCA_positive = gsub("a","",unique(best_alpha$patient[which(best_alpha$Subtype!="BRCA-WT")]))

# Subtype 2: "normal" breast cancers not BRCA+ (250 patients)
PR_ER_pos_HER2_neg = clinical$sample_name[which(clinical$final.PR=="positive"&clinical$final.ER=="positive"&clinical$final.HER2=="negative")]
PR_ER_pos_HER2_neg = PR_ER_pos_HER2_neg[which(PR_ER_pos_HER2_neg%in%gsub("a","",(unique(best_alpha$patient[which(best_alpha$Subtype=="BRCA-WT")]))))]

# Subtype 3: triple-negative breast cancers not BRCA+ (97 patients)
triple_negative = clinical$sample_name[which(clinical$final.PR=="negative"&clinical$final.ER=="negative"&clinical$final.HER2=="negative")]
triple_negative = triple_negative[which(triple_negative%in%gsub("a","",(unique(best_alpha$patient[which(best_alpha$Subtype=="BRCA-WT")]))))]

# Subtype 4: HER2-enriched not BRCA+ (71 patients)
HER2_positive = clinical$sample_name[which(clinical$final.HER2=="positive")]
HER2_positive = HER2_positive[which(HER2_positive%in%gsub("a","",(unique(best_alpha$patient[which(best_alpha$Subtype=="BRCA-WT")]))))]

# Subtype 5: all the other tumors (51 patients)
other_tumors = gsub("a","",unique(best_alpha$patient)[which(!gsub("a","",unique(best_alpha$patient))%in%c(BRCA_positive,PR_ER_pos_HER2_neg,triple_negative,HER2_positive))])

# set the subtypes
best_alpha$Subtype[which(gsub("a","",best_alpha$patient)%in%BRCA_positive)] = "BRCA+"
best_alpha$Subtype[which(gsub("a","",best_alpha$patient)%in%PR_ER_pos_HER2_neg)] = "PR_ER+_HER2-"
best_alpha$Subtype[which(gsub("a","",best_alpha$patient)%in%triple_negative)] = "TN"
best_alpha$Subtype[which(gsub("a","",best_alpha$patient)%in%HER2_positive)] = "HER2+"
best_alpha$Subtype[which(gsub("a","",best_alpha$patient)%in%other_tumors)] = "Other"

# split BRCA1 and BRCA2 (of the 91, BRCA1 are 52 and BRCA2 are 39)
for(i in 1:length(best_alpha$Subtype)) {
    if(best_alpha$Subtype[i]=="BRCA+") {
        best_alpha$Subtype[i] = best_alpha$Gene[i]
    }
}

best_alpha = melt(best_alpha, id.vars = c(1, 8, 9), variable.name = "signature")

# normalize alpha to sum to 1 for each patient
best_alpha[, norm:=value/sum(value), by = patient]

# plot without normalization (counts)
plt1 = ggplot(best_alpha) + geom_boxplot(aes(x = signature, y = value, fill = Subtype)) +ylim(c(0, 10000)) + ggtitle("Alpha values (raw)")
plt1

# plot normalized alphas for carriers and controls
plt2= ggplot(best_alpha) +  geom_boxplot(aes(x = signature, y = norm, fill = Subtype)) + ggtitle("Alpha values (normalized for each patient)")
plt2

# plot the number of mutations for the patients in each subtype
patients_counts = best_alpha[which(!duplicated(best_alpha$patient)),]
patients_counts$value = NULL
patients_counts$norm = NULL
patients_counts$signature = NULL
patients_counts = cbind(patients_counts,rep(0,nrow(patients_counts)))
colnames(patients_counts)[4] = "Mutations"
for(i in 1:nrow(patients_counts)) {
    patients_counts$Mutations[i] = sum(patients[patients_counts$patient[i],])
}
plt3 = ggplot(patients_counts) + geom_boxplot(aes(x = Subtype, y = Mutations, fill = Subtype)) + ylim(c(0, 25000)) + ggtitle("Number of point mutations")
plt3
