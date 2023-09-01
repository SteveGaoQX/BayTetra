#' Adjust Beta Samples and Estimate Based on Signal
#'
#' This function adjusts the beta samples and the beta estimate based on the provided
#' signal matrices. If `signal_kq` is 0 for any (k, q) pair, it sets the samples and
#' estimate (from the second to the last dimension) for that pair to 0. Similarly,
#' if `signal_kq0` is 0, it sets the samples (intercept) and estimate for the first dimension
#' of that pair to 0.
#'
#' @param post_beta_samples A 4-dimensional array containing the original beta samples.
#'   Dimensions should be (num_samples, K, Q, L).
#' @param beta_estimate A 3-dimensional array containing the beta estimate.
#'   Dimensions should be (K, Q, L).
#' @param signal_kq A matrix indicating signals for each k, q pair. The output from extract_signal function
#' @param signal_kq0 A matrix indicating signals for the first dimension for each k, q pair.The output from extract_signal function
#'
#' @return A list containing the adjusted `new_beta_samples` and `new_beta_estimate`.
#'
#' @examples
#' \dontrun{
#' beta_samples <- array(runif(8000), dim = c(100, 2, 2, 2))
#' beta_est <- array(runif(12), dim = c(2, 2, 2))
#' sig_kq <- matrix(c(0, 1, 1, 0), nrow = 2)
#' sig_kq0 <- matrix(c(1, 0, 0, 1), nrow = 2)
#' results <- adjust_beta_based_on_signal(beta_samples, beta_est, sig_kq, sig_kq0)
#' }
#' @noRd

update_beta_based_on_signal <- function(post_beta_samples, beta_estimate,
                                        signal_kq, signal_kq0) {
  # make a copy of the original beta samples
  new_beta_samples <- post_beta_samples
  new_beta_estimate <- beta_estimate

  # get the dimensions
  K <- dim(signal_kq)[1]
  Q <- dim(signal_kq)[2]
  L <- dim(new_beta_samples)[4]

  # loop over K and Q
  for (k in 1:K) {
    for (q in 1:Q) {
      # if signal_kq[k,q] is 0, set new_beta_samples[, k, q, 2:L] to 0
      if (signal_kq[k,q] == 0) {
        new_beta_samples[, k, q, 2:L] <- 0
        new_beta_estimate[k, q, 2:L] <- 0
      }

      # if signal_kq0[k,q] is 0, set new_beta_samples[, k, q, 1] to 0
      if (signal_kq0[k, q] == 0) {
        new_beta_samples[, k, q, 1] <- 0
        new_beta_estimate[k, q, 1] <- 0
      }
    }
  }

  return(list(new_beta_samples = new_beta_samples, new_beta_estimate = new_beta_estimate))
}
