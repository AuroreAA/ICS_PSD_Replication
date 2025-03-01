# ICS_PSD_Replication

## About ICS
Invariant Coordinate Selection (ICS) is a powerful unsupervised multivariate method designed to identify the structure of multivariate datasets on a subspace. It relies on the joint diagonalization of two scatter matrices and is particularly relevant as a dimension reduction tool prior to clustering or outlier detection.

More information can be found in our article:

Archimbaud, A. (2024). Generalized implementation of invariant coordinate selection with positive semi-definite scatter matrices. arXiv preprint arXiv:2409.02258.

## Reproduce results
This repository provides a collection of [R](https://CRAN.R-project.org/) 
scripts to reproduce all examples, simulations and figures in our manuscript 
and the supplementary report.  

The easiest way to reproduce the results is to clone this repository with 
[RStudio](https://rstudio.com/products/rstudio/download/).  Running the 
scripts within the resulting RStudio project ensures that there are no issues 
with file paths for storing or reading results, or for producing files 
containing plots.  In addition, the RStudio project uses 
[renv](https://rstudio.github.io/renv/) to make sure that the correct 
versions of all required R packages are used.  After opening the RStudio 
project for the fist time, please type `renv::restore()` on the R command 
line to retrieve the correct versions of all required packages.

