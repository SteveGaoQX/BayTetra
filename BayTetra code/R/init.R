
#' MCMC initialization
#'
#' This function initializes the parameters that need to be estimated in MCMC.
#'
#' @param Q Number of response.
#' @param S Number of dimensions in the alpha vector.
#' @param K Number of groups.
#' @param L length of vector beta.
#' @param I Number of subjects.
#' @param a_rho Shape hyperparameter for the beta distribution used to initialize rho (default is 2). refer to paper section 4.1
#' @param b_rho Rate hyperparameter for the beta distribution used to initialize rho (default is 2).
#' @param a_nu Shape hyperparameter for the gamma distribution used to initialize nu_kq and nu_kq0 (default is 5).
#' @param b_nu Rate hyperparameter for the gamma distribution used to initialize nu_kq and nu_kq0 (default is 25).
#' @param h_1 First parameter for the inverse gamma distribution used to initialize sigma2 (default is 1).
#' @param h_2 Second parameter for the inverse gamma distribution used to initialize sigma2 (default is 1).
#'
#' @return A list containing the initialized values of the parameters alpha, m_kq, xi_kq, rho, rho_0, nu_kq, nu_kq0, gamma_kq, gamma_kq0, eta_kq, beta_kq0, beta_wo_intcp, beta_whole, Sigma_omega, omega, and sigma2.
#'
#'
#' @useDynLib BayTetra, .registration = TRUE
#' @importFrom stats rbeta rgamma rnorm
#' @importFrom MCMCpack rinvgamma
#' @importFrom MASS mvrnorm
#' @importFrom Rcpp sourceCpp
#'
#' @noRd

init <- function( Q,S,K,L,I,
                   a_rho = 2,
                   b_rho = 2,
                   a_nu = 5,
                   b_nu = 25,
    h_1 = 1, h_2 =1){

  # MCMC initialization
  # initialize parameters need to be estimated in MCMC

  alpha_init <- rmvn_rcpp(Q, rep(0, S), diag(1, S))
  # Initialize m_kq
  m_kq <- array(0, dim=c(K, Q, L-1))
  for (k in 1:K){
    for (q in 1:Q){
      for (l in 1:(L-1)){
        m_kq[k, q, l] <- sample(c(-1, 1), 1)  # Bernoulli variables take only values 1 or -1
      }
    }
  }
  # Initialize xi_kq
  xi_kq <- array(0, dim=c(K, Q, L-1))
  for (k in 1:K){
    for (q in 1:Q){
      m_initial <- m_kq[k, q, ]
      xi_kq[k, q, ] <- MASS::mvrnorm(1, mu = m_initial, Sigma = diag(L-1))
    }
  }
  # Initialize rho
  # a_rho <- 2  # Assign a reasonable value
  # b_rho <- 2  # Assign a reasonable value
  rho <- rbeta(1, a_rho, b_rho)
  rho_0 = rbeta(1, a_rho, b_rho)


  # Define a_nu, b_nu
  # a_nu <- 5  # Assign a reasonable value
  # b_nu <- 25  # Assign a reasonable value
  #
  # Initialize nu_kq and nu_kq0
  nu_kq <- array(0, dim=c(K, Q))
  nu_kq0 <- array(0, dim=c(K, Q))
  for (k in 1:K){
    for (q in 1:Q){
      nu_kq[k, q] <- rgamma(1, shape = a_nu, rate = b_nu)
      nu_kq0[k, q] <- rgamma(1, shape = a_nu, rate = b_nu)
    }
  }

  # Initialize gamma_kq and gamma_kq0
  gamma_kq <- array(0, dim=c(K, Q))
  gamma_kq0 <- array(0, dim=c(K, Q))
  for (k in 1:K){
    for (q in 1:Q){
      gamma_kq[k, q] <- sample(c(nu_kq[k, q], 1), 1, prob = c(rho, 1-rho))
      gamma_kq0[k, q] <- sample(c(nu_kq0[k, q], 1), 1, prob = c(rho, 1-rho))
    }
  }

  # Initialize eta_kq and beta_kq0
  eta_kq <- array(0, dim=c(K, Q))
  beta_kq0 <- array(0, dim=c(K, Q))
  for (k in 1:K){
    for (q in 1:Q){
      eta_kq[k, q] <- rnorm(1, mean = 0, sd = sqrt(gamma_kq[k, q]*nu_kq[k, q]))
      beta_kq0[k, q] <- rnorm(1, mean = 0, sd = sqrt(gamma_kq0[k, q]*nu_kq0[k, q]))
    }
  }

  # Initialize beta_wo_intcp
  beta_wo_intcp <- array(0, dim=c(K, Q, L-1))
  for (k in 1:K){
    for (q in 1:Q){
      beta_wo_intcp[k, q, ] <- eta_kq[k, q] * xi_kq[k, q, ]
    }
  }
  # Initialize beta_whole
  beta_whole <- array(0, dim=c(K, Q, L))
  for (k in 1:K){
    for (q in 1:Q){
      beta_whole[k, q, 1] <- beta_kq0[k, q]
      beta_whole[k, q, 2:L] <- beta_wo_intcp[k, q, ]
    }
  }
  Sigma_omega_init <- diag(1, Q)
  omega_init <- matrix(NA, nrow=I, ncol=Q)
  for (i in 1:I){
    omega_init[i,] <- rmvn_rcpp(1, rep(0, Q), Sigma_omega_init)
  }

  sigma2_init <- rinvgamma(Q, h_1, h_2)

  init_list <- list(alpha=alpha_init,m_kq = m_kq, xi_kq = xi_kq,
                    rho = rho, rho_0 = rho_0,nu_kq =nu_kq,nu_kq0 = nu_kq0,
                    gamma_kq = gamma_kq, gamma_kq0 =gamma_kq0,
                    eta_kq = eta_kq, beta_kq0 = beta_kq0,
                    beta_wo_intcp = beta_wo_intcp, beta_whole = beta_whole,
                    Sigma_omega=Sigma_omega_init, omega=omega_init, sigma2=sigma2_init)
  return(init_list)
}





#' @noRd
init_Q1 <- function(Q,S,K,L,I,
                    a_rho = 2,
                    b_rho = 2,
                    a_nu = 5,
                    b_nu = 25,
                    h_1 = 1, h_2 =1){

  alpha_init <- rmvn_rcpp(1, rep(0, S), diag(1, S))
  # Initialize m_kq
  m_kq <- array(0, dim=c(K, L-1))
  for (k in 1:K){
    for (l in 1:(L-1)){
      m_kq[k, l] <- sample(c(-1, 1), 1)  # Bernoulli variables take only values 1 or -1
    }
  }
  # Initialize xi_kq
  xi_kq <- array(0, dim=c(K, L-1))
  for (k in 1:K){
    m_initial <- m_kq[k,]
    xi_kq[k, ] <- MASS::mvrnorm(1, mu = m_initial, Sigma = diag(L-1))
  }

  rho <- rbeta(1, a_rho, b_rho)
  rho_0 = rbeta(1, a_rho, b_rho)


  # Initialize nu_kq and nu_kq0
  nu_kq <- rep(0,K)
  nu_kq0 <- rep(0,K)
  for (k in 1:K){
    nu_kq[k] <- rgamma(1, shape = a_nu, rate = b_nu)
    nu_kq0[k] <- rgamma(1, shape = a_nu, rate = b_nu)
  }

  # Initialize gamma_kq and gamma_kq0
  gamma_kq <- rep(0,K)
  gamma_kq0 <- rep(0,K)
  for (k in 1:K){
    gamma_kq[k] <- sample(c(nu_kq[k], 1), 1, prob = c(rho, 1-rho))
    gamma_kq0[k] <- sample(c(nu_kq0[k], 1), 1, prob = c(rho, 1-rho))

  }

  # Initialize eta_kq and beta_kq0
  eta_kq <- rep(0,K)
  beta_kq0 <- array(0, K)
  for (k in 1:K){
    eta_kq[k] <- rnorm(1, mean = 0, sd = sqrt(gamma_kq[k]*nu_kq[k]))
    beta_kq0[k] <- rnorm(1, mean = 0, sd = sqrt(gamma_kq0[k]*nu_kq0[k]))
  }

  # Initialize beta_wo_intcp
  beta_wo_intcp <- array(0, dim=c(K, L-1))
  for (k in 1:K){

    beta_wo_intcp[k,  ] <- eta_kq[k] * xi_kq[k, ]

  }
  # Initialize beta_whole
  beta_whole <- array(0, dim=c(K,  L))
  for (k in 1:K){
    beta_whole[k,  1] <- beta_kq0[k]
    beta_whole[k,  2:L] <- beta_wo_intcp[k, ]
  }

  sigma2_init <- rinvgamma(1, h_1, h_2)

  init_list <- list(alpha=alpha_init,m_kq = m_kq, xi_kq = xi_kq,
                    rho = rho, rho_0 = rho_0,nu_kq =nu_kq,nu_kq0 = nu_kq0,
                    gamma_kq = gamma_kq, gamma_kq0 =gamma_kq0,
                    eta_kq = eta_kq, beta_kq0 = beta_kq0,
                    beta_wo_intcp = beta_wo_intcp, beta_whole = beta_whole,
                    sigma2=sigma2_init)
  return(init_list)
}







