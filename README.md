# BayTetra
BayTetra: A Bayesian Semiparametric Approach for Testing Trajectory Differences

BayTetra addresses an important task in the field of biomedical applications, testing differences in longitudinal trajectories among distinct groups of populations. The package offers a Bayesian semiparametric approach for modeling multivariate longitudinal data. It accounts for correlations
among different responses and employs B-splines, along with penalties on smoothness of the spline coefficients, for flexible and parsimonious trajectory estimation. The package is inspired by the research paper ''BayTetra — A Bayesian
Semiparametric Approach for Testing Trajectory Differences'' by Jin, W & Gao, Q & Xu, Y (2023).

## Steps to install the BayTetra
## 1.List of dependent packages
dep_packages <- c(
    "Rcpp", "RcppArmadillo", "MCMCpack", "MASS", "splines", "Matrix",
    "rstan", "mvtnorm", "truncnorm", "pracma", "GIGrvg"
)

## 2.Check if these necessary packages are already installed and install them if they are not
new.packages <- dep_packages[!(dep_packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) {
    install.packages(new.packages)
}

## 3.Now you can proceed to install the BayTetra package using BayTetra_0.1.0.tar.gz (using local file)
install.packages("BayTetra_0.1.0.tar.gz", repos = NULL), which is a source file, and the binary file name is "BayTetra_0.1.0_R_x86_64-pc-linux-gnu.tar.gz"

## 4.Example_code.R provides an example of how to quickly implements the BayTetra

## 5.if you want to use the ex_data inside the package(after installation), using:
library(BayTetra)
data("ex_data")


## 6.BayTetra.pdf is the Roxygen file introduced the details of BayTetra package

