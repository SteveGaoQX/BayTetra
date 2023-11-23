#' Extract Signal Information from Gamma Samples
#'
#' This function extracts and processes the signal information from gamma samples.
#' It calculates the probabilities based on the gamma samples and generates matrices
#' for the real and processed signal values.
#'
#' @param post_gamma_kq_samples posterior gamma_kq samples, in dim(num_samples, K,Q)
#' @param post_gamma_kq0_samples posterior gamma_kq0 samples indicating the signal for intercept. dim(num_samples, K,Q)
#' @param K The number of groups.
#' @param Q The number of response.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{signal_kq_real}: The probability of gamma_kq = 1, based on `post_gamma_kq_samples`.
#'     \item \code{signal_kq0_real}: The probability of gamma_kq0 = 1, based on `post_gamma_kq0_samples`.
#'     \item \code{signal_kq}: A matrix of processed signal values based on `post_gamma_kq_samples` (0 or 1).
#'     \item \code{signal_kq0}: A matrix of processed signal values based on `post_gamma_kq0_samples` (0 or 1).
#'   }
#'
#' @noRd

extract_signal <- function(post_gamma_kq_samples, post_gamma_kq0_samples, K, Q) {
  K = dim(post_gamma_kq_samples)[2]
  Q = dim(post_gamma_kq_samples)[3]
  signal_kq <- matrix(nrow = K, ncol = Q)
  signal_kq0 <- matrix(nrow = K, ncol = Q)
  signal_kq_real <- matrix(nrow = K, ncol = Q)
  signal_kq0_real <- matrix(nrow = K, ncol = Q)

  for (q in 1:Q) {
    for (k in 1:K) {
      prob_kq <- mean(post_gamma_kq_samples[ , k,q] == 1)
      prob_kq0 <- mean(post_gamma_kq0_samples[ , k,q] == 1)
      signal_kq_real[k,q] = prob_kq
      signal_kq0_real[k,q] = prob_kq0
      signal_kq[k, q] <- ifelse(prob_kq > 0.5, 1, 0)
      signal_kq0[k, q] <- ifelse(prob_kq0 > 0.5, 1, 0)
    }
  }

  return(list(signal_kq_real = signal_kq_real,
              signal_kq0_real = signal_kq0_real,
              signal_kq = signal_kq,
              signal_kq0 = signal_kq0))
}
