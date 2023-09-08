# BayTetra
BayTetra: A Bayesian Semiparametric Approach for Testing Trajectory Differences

BayTetra addresses an important task in the field of biomedical applications, testing differences in longitudinal trajectories among distinct groups of populations. The package offers a Bayesian semiparametric approach for modeling multivariate longitudinal data. It accounts for correlations
among different responses and employs B-splines, along with spike-and-slab priors on the spline coefficients, for flexible and parsimonious trajectory estimation. The package is inspired by the research paper ''BayTetra — A Bayesian
Semiparametric Approach for Testing Trajectory Differences'' by Wei, J & Gao, Q & Xu, Y (2023).

## Steps to install the BayTetra
### List of dependent packages
dep_packages <- c(
    "Rcpp", "RcppArmadillo", "MCMCpack", "MASS", "splines", "dplyr", "tmvtnorm", "Matrix",
    "rstan", "mvtnorm", "truncnorm", "pracma", "loo"
)

### Check if the packages are already installed and install them if they are not
new.packages <- dep_packages[!(dep_packages %in% installed.packages()[,"Package"])]


if(length(new.packages)) {
    install.packages(new.packages)
}

### Now you can proceed to install the BayTetra package using BayTetra_0.1.0.tar.gz
install.packages("BayTetra_0.1.0.zip", repos = NULL)

## Example_code.R provides an example of how to quickly implements the BayTetra

## BayTetra.pdf is the Roxygen file introduced the details of BayTetra package

