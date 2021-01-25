#' @name patients 
#' @title point mutations for 560 breast tumors 
#' @description dataset of counts of the point mutations detected in 560 breast tumors published in Nik-Zainal, Serena, et al. (2016). 
#' @docType data 
#' @usage data(patients) 
#' @format counts of the point mutations 
#' @source Nik-Zainal, Serena, et al. "Landscape of somatic mutations in 560 breast cancer whole-genome sequences." Nature 534.7605 (2016): 47. 
#' @return counts of point mutations for 560 tumors and 96 trinucleotides 
NULL

#' @name background 
#' @title germline replication error 
#' @description germline replication error estimated in Rahbari, Raheleh, et al. (2016). 
#' @docType data 
#' @usage data(background) 
#' @format vector of rates 
#' @source Rahbari, Raheleh, et al. "Timing, rates and spectra of human germline mutation." Nature genetics 48.2 (2016): 126. 
#' @return vector of rates for the 96 trinucleotides 
NULL

#' @name background2 
#' @title COSMIC replication error 
#' @description background replication error signature derived from COSMIC SBS5. 
#' @docType data 
#' @usage data(background2) 
#' @format vector of rates 
#' @source COSMIC database (https://cancer.sanger.ac.uk/cosmic/signatures) v3. 
#' @return vector of rates for the 96 trinucleotides 
NULL

#' @name mutation_categories 
#' @title trinucleotides mutation categories 
#' @description 96 trinucleotides mutation categories 
#' @docType data 
#' @usage data(mutation_categories) 
#' @format matrix of 96 trinucleotides mutation categories 
#' @return matrix of 96 trinucleotides mutation categories 
NULL

#' @name ssm560_reduced 
#' @title a reduced version of the point mutations for 560 breast tumors in the format compatible with the import function 
#' @description reduced versione of the dataset of counts of the point mutations detected in 560 breast tumors published in Nik-Zainal, Serena, et al. (2016). 
#' @docType data 
#' @usage data(ssm560_reduced) 
#' @format reduced versione of the counts of the point mutations in the format compatible with the import function 
#' @source Nik-Zainal, Serena, et al. "Landscape of somatic mutations in 560 breast cancer whole-genome sequences." Nature 534.7605 (2016): 47. 
#' @return reduced versione of the counts of point mutations for 560 tumors and 96 trinucleotides in the format compatible with the import function 
NULL

#' @name starting_betas_example 
#' @title example of results obtained with the function starting.betas.estimation on the counts input from Nik-Zainal, Serena, et al. (2016). 
#' @description example of results obtained with the function starting.betas.estimation on the counts input from Nik-Zainal, Serena, et al. (2016). 
#' @docType data 
#' @usage data(starting_betas_example) 
#' @format results obtained with the function starting.betas.estimation on the counts input from Nik-Zainal, Serena, et al. (2016) 
#' @return results obtained with the function starting.betas.estimation on the counts input from Nik-Zainal, Serena, et al. (2016) 
NULL

#' @name lambda_range_example 
#' @title example of results obtained with the function evaluate.lambda.range on the counts input from Nik-Zainal, Serena, et al. (2016). 
#' @description example of results obtained with the function evaluate.lambda.range on the counts input from Nik-Zainal, Serena, et al. (2016). 
#' @docType data 
#' @usage data(lambda_range_example) 
#' @format results obtained with the function evaluate.lambda.range on the counts input from Nik-Zainal, Serena, et al. (2016) 
#' @return results obtained with the function evaluate.lambda.range on the counts input from Nik-Zainal, Serena, et al. (2016) 
NULL

#' @name cv_example 
#' @title example of results obtained with the function nmf.LassoCV on the counts input from Nik-Zainal, Serena, et al. (2016). 
#' @description example of results obtained with the function nmf.LassoCV on the counts input from Nik-Zainal, Serena, et al. (2016). 
#' @docType data 
#' @usage data(cv_example) 
#' @format results obtained with the function nmf.LassoCV on the counts input from Nik-Zainal, Serena, et al. (2016) 
#' @return results obtained with the function nmf.LassoCV on the counts input from Nik-Zainal, Serena, et al. (2016) 
NULL

#' @name nmf_LassoK_example 
#' @title example of results obtained with the function nmf.LassoK on the counts input from Nik-Zainal, Serena, et al. (2016). 
#' @description example of results obtained with the function nmf.LassoK on the counts input from Nik-Zainal, Serena, et al. (2016). 
#' @docType data 
#' @usage data(nmf_LassoK_example) 
#' @format results obtained with the function nmf.LassoK on the counts input from Nik-Zainal, Serena, et al. (2016) 
#' @return results obtained with the function nmf.LassoK on the counts input from Nik-Zainal, Serena, et al. (2016) 
NULL

#' @name imported_data 
#' @title example of imported data from Nik-Zainal, Serena, et al. (2016). 
#' @description example of imported data from Nik-Zainal, Serena, et al. (2016). 
#' @docType data 
#' @usage data(imported_data) 
#' @format results obtained with the function import.counts.data on the data from Nik-Zainal, Serena, et al. (2016) 
#' @return results obtained with the function import.counts.data on the data from Nik-Zainal, Serena, et al. (2016) 
NULL
