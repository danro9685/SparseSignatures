SparseSignatures
================

| Branch | Status |
| --- | --- |
| master | [![R-CMD-check-bioc](https://github.com/danro9685/SparseSignatures/actions/workflows/check-bioc.yml/badge.svg?branch=master)](https://github.com/danro9685/SparseSignatures/actions/workflows/check-bioc.yml) |
| development | [![R-CMD-check-bioc](https://github.com/danro9685/SparseSignatures/actions/workflows/check-bioc.yml/badge.svg?branch=development)](https://github.com/danro9685/SparseSignatures/actions/workflows/check-bioc.yml) |

Point mutations occurring in a genome can be divided into 96 categories based on the base being mutated, the base it is mutated into and its two flanking bases. Therefore, for any patient, it is possible to represent all the point mutations occurring in that patient's tumor as a vector of length 96, where each element represents the count of mutations for a given category in the patient. 

A mutational signature represents the pattern of mutations produced by a mutagen or mutagenic process inside the cell. Each signature can also be represented by a vector of length 96, where each element represents the probability that this particular mutagenic process generates a mutation of the 96 above mentioned categories. In this R package, we provide a set of functions to extract and visualize the mutational signatures that best explain the mutation counts of a large number of patients. 
