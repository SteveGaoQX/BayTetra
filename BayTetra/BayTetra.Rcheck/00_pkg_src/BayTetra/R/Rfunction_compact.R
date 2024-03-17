


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
#' 
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
#' @return A matrix of updated nu values. dim(K,Q)
#' @noRd

update_nu_kq = function(a_nu, b_nu, beta_kq0,intercept = FALSE){
  K = dim(beta_kq0)[1]
  Q = dim(beta_kq0)[2]
  nu_update = array(0, dim = c(K,Q))
  if (intercept == TRUE){
    K = K-1
  }
  for (k in 1:K) {
    for (q in 1:Q) {
      shape_updated <- a_nu + 0.5
      scale_updated <- b_nu + (beta_kq0[k,q]^2)/2
      
      # Sample from the Gamma distribution
      gamma_sample <- rgamma(1, shape = shape_updated, rate = scale_updated)
      nu_update[k,q] = 1/(gamma_sample)
    }
  }
  
  return(nu_update)
}


#' @importFrom stats rgamma
#' @importFrom GIGrvg rgig
NULL
#' Update eta Parameter
#' @noRd
#' 
update_tau_kq_gamma = function(alpha_g,beta_g, beta, K_mat, intercept = FALSE){
  K = dim(beta)[1]
  Q = dim(beta)[2]
  L = dim(beta)[3]
  
  tau_update = array(0, dim = c(K,Q))
  # if (intercept == TRUE){
  #   K = K-1
  # }
  for (k in 1:K) {
    for (q in 1:Q) {
      beta_kq_vec = beta[k,q,]
      mid = t(beta_kq_vec)%*%K_mat%*%beta_kq_vec
      
      lambda <- alpha_g - 0.5*(L-2)
      chi <- mid[1]
      psi <- 2*beta_g
      
      # Sample from the Gamma distribution
      tau_update[k,q] = rgig(n = 1, lambda = lambda, chi = chi, psi = psi)
      
    }
  }
  
  return(tau_update)
}




