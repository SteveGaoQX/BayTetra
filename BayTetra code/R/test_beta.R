#' Test Hypothesis for Beta and Estimate its Value
#'
#' This function performs a hypothesis test on the input beta samples
#' It checks the 95% confidence interval (calculated by posterior samples) will cover zero.
#' If not, the value is considered significant. The function also returns
#' an estimation of beta, calculated as the mean value of the input beta values.
#'
#' @param beta A 4-dimensional array of beta values. The dimensions are
#'   assumed to be (num_samples, K, Q, L).
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{beta_test_result}: A 3-dimensional logical array with dimensions (K, Q, L).
#'     Each value indicates whether the corresponding beta value is significant (TRUE) or not (FALSE).
#'     \item \code{beta_estimation}: A 3-dimensional array of the mean beta values with dimensions (K, Q, L).
#'   }
#'
#'
#' @importFrom stats quantile
#' @noRd

test_beta <- function(beta) {
  # hypothesis test and get the estimation
  # Check that the array is of the correct dimensions
  if (length(dim(beta)) != 4) {
    stop("Beta array should have dimensions (num_samples, K, Q, L)")
  }

  # Get the dimensions
  num_samples <- dim(beta)[1]
  K <- dim(beta)[2]
  Q <- dim(beta)[3]
  L <- dim(beta)[4]

  # Pre-allocate result arrays
  beta_test_result <- array(dim = c(K, Q, L))
  beta_estimation <- array(dim = c(K, Q, L))

  # Calculate mean and standard deviation
  beta_mean <- apply(beta, c(2,3,4), mean)

  # Calculate the 95% confidence interval
  ci_lower <- apply(beta, c(2,3,4), function(x) quantile(x, 0.025))
  ci_upper <- apply(beta, c(2,3,4), function(x) quantile(x, 0.975))

  # Test hypothesis
  beta_test_result <- !(ci_lower <= 0 & ci_upper >= 0)

  # Store mean as estimation
  beta_estimation <- beta_mean

  list(beta_test_result = beta_test_result, beta_estimation = beta_estimation)
}
