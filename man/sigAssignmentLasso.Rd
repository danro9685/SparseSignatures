% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signatures.assignment.lasso.R
\name{sigAssignmentLasso}
\alias{sigAssignmentLasso}
\title{sigAssignmentLasso}
\usage{
sigAssignmentLasso(
  x,
  beta,
  normalize_counts = TRUE,
  lambda_rate_alpha = 0.05,
  max_iterations_lasso = 10000,
  seed = NULL,
  verbose = TRUE
)
}
\arguments{
\item{x}{count matrix for a set of n patients and 96 trinucleotides.}

\item{beta}{beta to be fixed during the estimation of alpha.}

\item{normalize_counts}{if true, the input count matrix x is normalize such that the patients have the same number of mutation.}

\item{lambda_rate_alpha}{value of LASSO to be used for alpha between 0 and 1. This value should be greater than 0. 1 is the value of LASSO 
that would shrink all the exposure values to 0 within one step. The higher lambda_rate_alpha is, the sparser are the resulting exposure values, 
but too large values may result in a reduced fit of the observed counts.}

\item{max_iterations_lasso}{Number of maximum iterations to be performed during the sparsification via Lasso.}

\item{seed}{Seed for reproducibility.}

\item{verbose}{boolean; Shall I print all messages?}
}
\value{
A list with the discovered signatures and their assignment to patients. It includes 2 elements: 
             alpha: matrix of the assigned exposure values
             beta: matrix of the discovered signatures
}
\description{
Perform the assignment of somatic mutational signatures to patients given a set of observed counts x and signatures beta.
}
\examples{
data(patients)
data(starting_betas_example)
beta = starting_betas_example[["5_signatures","Value"]]
res = sigAssignmentLasso(x=patients[1:100,],beta=beta,lambda_rate_alpha=0.05,seed=12345)

}
