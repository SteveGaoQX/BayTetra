#' Standardize the eta and xi
#'
#' This function updates (standardize) the values of parameters xi, and eta (without intercept), refer to the MCMC
#' process in the supplements
#'
#' @param eta_update A numeric array representing the updated eta values.
#' @param xi_update A numeric array representing the updated xi values.
#'
#' @return A list containing the following updated parameters:
#' \describe{
#'  \item{current_beta_wo_intcp}{An array of the updated beta values without intercept.}
#'  \item{xi_update_std}{An array of the updated and standardized xi values.}
#'  \item{eta_update_std}{An array of the updated and standardized eta values.}
#' }
#' @noRd

update_beta_wo_eta_xi_std <- function(eta_update, xi_update) {
  # get dimensions
  K = dim(xi_update)[1]
  Q = dim(xi_update)[2]
  L = dim(xi_update)[3] + 1

  # Initialize output arrays
  current_beta_wo_intcp = array(0,dim = c(K,Q,L-1))
  xi_update_std = array(0, dim = c(K,Q,L-1))
  eta_update_std = array(0, dim = c(K,Q))

  # Start the loops
  for (k in 1:K) {
    for (q in 1:Q) {
      current_beta_wo_intcp[k,q,] = eta_update[k,q]*xi_update[k,q, ]
      xi_kq_bar = mean(abs(xi_update[k,q,])) # xi_update length should be L-1
      xi_update_std[k,q,] = xi_update[k,q,]/xi_kq_bar
      eta_update_std[k,q] = eta_update[k,q]*xi_kq_bar
    }
  }

  # return list with the updated variables
  return(list(current_beta_wo_intcp = current_beta_wo_intcp,
              xi_update_std = xi_update_std,
              eta_update_std = eta_update_std))
}


#' Update Gamma Parameter
#'
#' This function updates the gamma parameter values used in the likelihood estimation of a multivariate normal distribution model.
#'
#' @param eta A numeric matrix representing the eta values.
#' @param rho A numeric scalar representing the rho value.
#' @param nu_0 A numeric scalar representing the initial nu value.
#' @param nu A numeric matrix representing the current nu values.
#' @param intercept Logical, default is FALSE. If TRUE, one dimension of K is subtracted in the loop.
#'
#' @return A matrix of updated gamma values.
#' @noRd
update_gamma_kq <- function(eta, rho,nu_0,nu, intercept = FALSE){
  # update gamma
  # args: eta (K,Q),rho(scalar),nu (K*Q), nu_0 (scalar)
  # returns: update_gamma_kq
  K = dim(eta)[1]
  Q = dim(eta)[2]
  update_gamma <- array(0, dim=c(K, Q))
  if(intercept == TRUE){

    K = K-1

  }
  # update gamma, k = 1...K
  for (k in 1:K) {
    for (q in 1:Q) {
      exp_num = (1 - nu_0)*(eta[k,q]^2)/(2*nu_0*nu[k,q])
      if (exp_num>500){
        update_gamma[k, q] = 1
      }else{
        # Compute the new gamma value
        gamma_value <- (sqrt(nu_0) * rho / (1 - rho)) * exp(((1 - nu_0) * eta[k, q]^2) / (2 * nu_0 * nu[k, q]))
        # Compute probabilities A and B
        p0 <- 1 / (1 + gamma_value)
        p1 <- 1 - p0
        # Sample gamma[k, q] from a Bernoulli distribution with probability p1
        update_gamma[k, q] <- sample(c(nu_0, 1),size=1, prob = c(p0,p1))
      }
    }
  }

  return(update_gamma)
}



#' @importFrom stats rgamma
NULL
#' Update Nu Parameter
#'
#' This function updates the nu parameter values
#'
#' @param a_nu A numeric scalar representing the 'a' parameter for the Gamma distribution.
#' @param b_nu A numeric scalar representing the 'b' parameter for the Gamma distribution.
#' @param eta A numeric matrix representing the eta values.
#' @param gamma_kq A numeric matrix representing the gamma values.
#' @param intercept Logical, default is FALSE. If TRUE, one dimension of K is skiped in the loop.
#'
#' @return A matrix of updated nu values. dim(K,Q)
#' @noRd
update_nu_kq = function(a_nu, b_nu, eta, gamma_kq,intercept = FALSE){
  K = dim(gamma_kq)[1]
  Q = dim(gamma_kq)[2]
  nu_update = array(0, dim = c(K,Q))
  if (intercept == TRUE){
    K = K-1
  }
  for (k in 1:K) {
    for (q in 1:Q) {
      shape_updated <- a_nu + 0.5
      scale_updated <- b_nu + (eta[k,q]^2)/(2*gamma_kq[k,q])

      # Sample from the Gamma distribution
      gamma_sample <- rgamma(1, shape = shape_updated, rate = scale_updated)
      nu_update[k,q] = 1/(gamma_sample)
    }
  }

  return(nu_update)
}

#' @importFrom stats rbeta
NULL
#' Update Rho Parameter
#'
#' This function updates the rho parameter values used in the likelihood estimation of a multivariate normal distribution model.
#'
#' @param a_rho hyperparameter for the Beta distribution.
#' @param b_rho hyperparameter for the Beta distribution.
#' @param gamma_kq A numeric matrix representing the gamma values.
#' @param intercept Logical, default is FALSE. If TRUE, one dimension of K is skiped in the loop.
#'
#' @return A numeric value, the updated rho.
#' @noRd
update_rho_kq <- function(a_rho, b_rho, gamma_kq, intercept = FALSE) {
  K = dim(gamma_kq)[1]
  # Adjust K if intercept is TRUE
  if(intercept) {
    K <- K - 1
  }
  gamma_kq <- gamma_kq[1:K, ]

  shape1_updated <- a_rho + sum(gamma_kq == 1)
  shape2_updated <- b_rho + sum(gamma_kq != 1)
  rho_update<- rbeta(1, shape1_updated, shape2_updated)
  # Return the updated rho matrix
  return(rho_update)
}


#' Update m_kql Matrix
#'
#' This function updates the values of the m_kql
#'
#' @param xi A 3-dimensional array representing the xi values.
#'
#' @return A 3-dimensional array, the updated 'm' matrix dim(K,Q,L-1).
#' @noRd
update_mkql = function(xi){
  K = dim(xi)[1]
  Q = dim(xi)[2]
  L = dim(xi)[3] # actually L-1 in formula
  m_update = array(0, dim = c(K,Q,L))
  for (k in 1:K) {
    for (q in 1:Q) {
      for (l in 1:L) {
        pos_prob = 1/(1+exp(-2*xi[k,q,l]))
        neg_prob = 1- pos_prob
        m_update[k, q,l] <- sample(c(-1,1),size = 1, prob = c(neg_prob,pos_prob))
      }
    }
  }
  return(m_update)
}

#' Update the Whole Beta Matrix
#'
#' This function updates the values of the entire beta matrix by combining the beta intercepts with the beta matrix excluding the intercept.
#'
#' @param beta_kq0 A matrix representing the beta intercepts.
#' @param beta_wo_intcp A 3-dimensional array representing the beta values excluding the intercepts.
#'
#' @return A 3-dimensional array, the updated 'beta' matrix dim(K,Q,L).
#' @noRd
update_whole_beta = function(beta_kq0, beta_wo_intcp){
  K = dim(beta_kq0)[1]
  Q = dim(beta_wo_intcp)[2]
  L = dim(beta_wo_intcp)[3]+1
  beta_whole_new = array(0, c(K,Q,L))
  for (k in 1:K) {
    beta_kq0_col <- matrix(beta_kq0[k,], ncol = 1)  # convert to a column matrix
    beta_whole_new[k, , ] <- cbind(beta_kq0_col, beta_wo_intcp[k, , ])  # combine the two matrices
  }
  return(beta_whole_new)
}





#' @noRd
update_beta_wo_eta_xi_std_Q1 <- function(eta_update, xi_update) {
  # get dimensions
  K = dim(xi_update)[1]
  L = dim(xi_update)[2] + 1

  # Initialize output arrays
  current_beta_wo_intcp = array(0,dim = c(K,L-1))
  xi_update_std = array(0, dim = c(K,L-1))
  eta_update_std =rep(0,K)

  # Start the loops
  for (k in 1:K) {

    current_beta_wo_intcp[k,] = eta_update[k]*xi_update[k, ]
    xi_kq_bar = mean(abs(xi_update[k,])) # xi_update length should be L-1
    xi_update_std[k,] = xi_update[k,]/xi_kq_bar
    eta_update_std[k] = eta_update[k]*xi_kq_bar

  }

  # return list with the updated variables
  return(list(current_beta_wo_intcp = current_beta_wo_intcp,
              xi_update_std = xi_update_std,
              eta_update_std = eta_update_std))
}



#' @noRd
update_gamma_kq_Q1 <- function(eta, rho,nu_0,nu, intercept = FALSE){
  # update gamma
  # args: eta (K,Q),rho(scalar),nu (K*Q), nu_0 (scalar)
  # returns: update_gamma_kq
  K = length(eta)
  update_gamma <- rep(0,K)
  if(intercept == TRUE){
    K = K - 1
  }
  # update gamma, k = 1...K
  for (k in 1:K) {

    exp_num = (1 - nu_0)*(eta[k]^2)/(2*nu_0*nu[k])
    if (exp_num>500){
      update_gamma[k] = 1
    }else{
      # Compute the new gamma value
      gamma_value <- (sqrt(nu_0) * rho / (1 - rho)) * exp(exp_num)
      # Compute probabilities A and B
      p1 <- 1 / (1 + gamma_value)
      p0 <- 1 - p1
      # Sample gamma[k, q] from a Bernoulli distribution with probability p1
      update_gamma[k] <- sample(c(nu_0, 1),size=1, prob = c(p1,p0))
    }

  }
  return(update_gamma)
}


#' @noRd
update_nu_kq_Q1 = function(a_nu, b_nu, eta, gamma_kq,intercept = FALSE){
  K = length(gamma_kq)
  nu_update = rep(0,K)
  if (intercept == TRUE){
    K = K - 1
  }
  for (k in 1:K) {
    shape_updated <- a_nu + 0.5
    scale_updated <- b_nu + (eta[k]^2)/(2*gamma_kq[k])

    # Sample from the Gamma distribution
    gamma_sample <- rgamma(1, shape = shape_updated, rate = scale_updated)
    nu_update[k] = 1/(gamma_sample)
  }
  return(nu_update)
}

#' @noRd
update_rho_kq_Q1 <- function(a_rho, b_rho, gamma_kq, intercept = FALSE) {
  K = length(gamma_kq)
  # Adjust K if intercept is TRUE
  if(intercept) {
    K <- K - 1
  }
  # Define the indicator functions delta_1 and delta_nu_0. Replace conditions accordingly.
  delta_1 <- function(x) ifelse(x == 1, 1, 0)     # replace condition accordingly
  delta_nu_0 <- function(x) ifelse(x != 1, 1, 0)  # replace condition accordingly
  # Only use the first K rows of gamma_kq
  gamma_kq <- gamma_kq[1:K]
  shape1_updated <- a_rho + sum(gamma_kq == 1)
  #print(sum(apply(gamma_kq, c(1,2), delta_1)))
  shape2_updated <- b_rho + sum(gamma_kq != 1)
  rho_update<- rbeta(1, shape1_updated, shape2_updated)
  # Return the updated rho matrix
  return(rho_update)
}

#' @noRd
update_mkql_Q1 = function(xi){
  K = dim(xi)[1]
  L = dim(xi)[2] # actually L-1 in formula
  m_update = array(0, dim = c(K,L))
  for (k in 1:K) {
    for (l in 1:L) {
      pos_prob = 1/(1+exp(-2*xi[k,l]))
      neg_prob = 1- pos_prob
      m_update[k, l] <- sample(c(-1,1),size = 1, prob = c(neg_prob,pos_prob))
    }
  }
  return(m_update)
}


#' @noRd
update_whole_beta_Q1 = function(beta_kq0, beta_wo_intcp){
  K = dim(beta_kq0)[1]
  L = dim(beta_wo_intcp)[2]+1
  beta_whole_new = array(0, c(K,L))
  beta_whole_new = cbind(beta_kq0, beta_wo_intcp)
  return(beta_whole_new)
}



