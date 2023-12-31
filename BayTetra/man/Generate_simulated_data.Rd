% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generate_simulated_data.R
\name{Generate_simulated_data}
\alias{Generate_simulated_data}
\title{Generate Simulated Data}
\usage{
Generate_simulated_data()
}
\value{
Longitudinal data with following variables:
\itemize{
\item \strong{ID}: Identity for individuals
\item \strong{VISIT}: Individuals' visit index
\item \strong{time}: Time variable
\item \strong{cov1}: Covariate
\item \strong{Group}: Group memberships for individuals
\item \strong{R1}: Response variable 1
\item \strong{R2}: Response variable 2
}
}
\description{
This function generates a simulated dataset used in Scenario #1 of the paper "BayTetra: A Bayesian Semiparametric Approach for Testing Trajectory Differences"
}
\examples{
\dontrun{
ex_data = Generate_simulated_data()
head(ex_data)
}

}
