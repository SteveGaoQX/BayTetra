#' Estimate Alpha from Posterior Alpha Samples
#'
#' This function calculates the mean value of the posterior alpha samples
#' across the number of samples for each combination of Q and S dimensions.
#'
#' @param post_alpha_samples A 3-dimensional array containing the posterior alpha samples.
#'   Dimensions should be (num_samples, Q, S).
#'
#' @return A matrix containing the mean values for each combination of Q and S.
#'
#' @noRd

estimate_alpha <- function(post_alpha_samples) {
  num_samples <- dim(post_alpha_samples)[1]
  Q <- dim(post_alpha_samples)[2]
  S <- dim(post_alpha_samples)[3]

  # Initialize estimation_alpha matrix
  estimation_alpha <- matrix(0, nrow = Q, ncol = S)

  # Loop through each q in Q and s in S
  for(q in 1:Q) {
    for(s in 1:S) {
      # Calculate mean of post_alpha_samples[ , q, s]
      estimation_alpha[q, s] <- mean(post_alpha_samples[ , q, s])
    }
  }

  return(estimation_alpha)
}
