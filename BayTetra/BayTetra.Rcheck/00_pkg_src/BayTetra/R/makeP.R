#' Make Penalty matrix
#'
#' This function initializes the parameters that need to be estimated in MCMC.
#'
#' @param dim Dimension.
#' @param degree Degree of differentiation.
#' @param epsilon Small number to ensure non-singularity.
#' @return Penalty matrix.
#' @noRd
makeP <- function(dim, degree, epsilon=1e-3){
  D <- diff(diag(dim), differences=degree)
  return(t(D)%*%D + diag(dim)*epsilon)
}





