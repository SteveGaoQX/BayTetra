#' Generate Simulated Data
#'
#' This function generates a simulated dataset used in Scenario #1 of the paper "BayTetra: A Bayesian Semiparametric Approach for Testing Trajectory Differences"
#'
#'
#' @return Longitudinal data with following variables:
#'   - \strong{ID}: Identity for individuals
#'   - \strong{VISIT}: Individuals' visit index
#'   - \strong{time}: Time variable
#'   - \strong{cov1}: Covariate
#'   - \strong{Group}: Group memberships for individuals
#'   - \strong{R1}: Response variable 1
#'   - \strong{R2}: Response variable 2
#'
#' @examples
#' \dontrun{
#' ex_data = Generate_simulated_data()
#' head(ex_data)
#' }
#'
#' @useDynLib BayTetra, .registration = TRUE
#' @import truncnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom stats runif
#' @importFrom splines bs
#' @importFrom pracma householder
#' @importFrom stats predict quantile
#' @export
#'
#'
#'
Generate_simulated_data <- function() {
  I <- 1000
  Q <- 2
  J_max <- 10
  K <- 3
  # n = I
  # Create storage for data

  y <- array(NA, c(I, Q, J_max))
  data_index <- array(1, c(I, Q, J_max))
  J_max <- 10 # the maximum number of observations for all the subjects
  J <- sample(3:(J_max-3), I, replace = TRUE) # the maximum number of observations for i-th subject
  # J_max <- 10 # the maximum number of observations for all the subjects
  # J <- sample(5:J_max, I, replace = TRUE) # the maximum number of observations for i-th subject
  data_index <- array(0, dim=c(I, Q, J_max)) # binary indicator of observed data (i.e., missing = 0 and observed = 1)
  for (i in 1:I){
    number_of_data <- J[i] # number of observations for i-th subject at q-th response
    data_index_iq <- sort(sample(1:J_max, number_of_data, replace = FALSE))
    for (q in 1:Q){
      data_index[i,q,data_index_iq] <- rep(1, number_of_data)
    }
  }
  data_index = array(as.integer(data_index), dim = dim(data_index))
  t = Generate_t_cpp(I,Q,J_max,0.2,-2,2, data_index)
  t_min <- min(t,na.rm = TRUE)
  t_max <- max(t,na.rm = TRUE)



  # hist(t)
  S <- 2 # dimension of covariates Z_{iqj}
  Z <- array(NA, dim=c(I, Q, J_max, S))
  # Generate values outside the for loop
  matrix_ones <- matrix(1, nrow=I, ncol=J_max)
  matrix_rmvn <- rmvn_rcpp(I, rep(0, J_max), diag(1, J_max))
  # Fill the array Z with the same values for each q
  for (q in 1:Q) {
    # intercept
    Z[, q, , 1] <- matrix_ones
    # time-invariant binary covariate
    Z[, q, , 2] <- matrix_rmvn
  }
  K <- 3 # number of subject classes
  g <- rep(1:K, length.out = I)
  df = 9
  # Here g should follow the index by Rcpp start from 0
  num_intervals <- df - 3 + 1
  # Calculate the quantiles
  quantiles <- quantile(t, probs = seq(0, 1, 1/num_intervals), na.rm = TRUE)
  # Exclude the first and last quantiles
  middle_quantiles <- quantiles[-c(1, num_intervals + 1)]
  # print(middle_quantiles)
  spline_basis <- bs(seq(t_min, t_max, length.out=1000), knots = middle_quantiles, intercept = TRUE)


  L <- dim(spline_basis)[2]  # degree of freedoms for spline expansion with intercept
  B <- array(NA, dim=c(I, Q, J_max, L)) # time spline basis matrix 
  for (i in 1:I){
    for (q in 1:Q){
      data_index_iq <- which(data_index[i,q,] == 1) # index of data for i-th subject at q-th response
      B[i,q,data_index_iq,] <- cbind( predict(spline_basis, t[i,q,data_index_iq]))
    }
  }
  #
  Z_qr_all = list()
  # Initialize B_after_qr_allq
  B_after_qr_allq <- array(NA, dim=c(I, Q, J_max, L-1))
  
  for (q in 1:Q){
    B_q <- NULL
    record_counter <- 1
    record_rows <- list()
    
    for(i in 1:I){
      data_index_iq <- which(data_index[i, q,] == 1) # Indices of data for i-th subject at q-th response
      B_iq <- B[i,q,data_index_iq,] # Extract data for current i and q
      #B_iq dim is J_iq*(L-1)
      if(length(data_index_iq) > 0){
        B_q <- rbind(B_q, B_iq) # Combine row-wise
        # Store the start and end row indices for each record
        record_rows[[i]] <- c(record_counter, record_counter + length(data_index_iq) - 1)
        record_counter <- record_counter + length(data_index_iq)
      } else {
        record_rows[[i]] <- c(NA, NA)
      }
    }
    
    # Compute B_q^T %*% 1
    trans1 <- t(B_q) %*% rep(1, nrow(B_q))
    
    # Perform Householder QR decomposition
    QR <- householder(trans1)
    Q_matrix <- QR$Q
    R_matrix <- QR$R
    
    # Form Z_qr from 2 to L columns of Q
    Z_qr <- Q_matrix[, -1]
    Z_qr_all[[q]] = Z_qr
    # Compute B_after_qr
    B_after_qr <- B_q %*% Z_qr
    
    # For each i, assign corresponding rows to B_after_qr_complete
    for(i in 1:I){
      #B_after_qr_complete <- array(NA, dim=c(n, J_max, L-1))
      if(!is.na(record_rows[[i]][1])){
        start_row <- record_rows[[i]][1]
        end_row <- record_rows[[i]][2]
        data_index_iq = which(data_index[i,q,] == 1)
        B_after_qr_allq[i,q, data_index_iq , ] = B_after_qr[start_row:end_row, ]
      }
    }
  }
  
  
  B = B_after_qr_allq
  L = dim(B)[4]
  
  P_beta <- makeP(dim=L, degree=2) 
  
  
  omega <- matrix(NA, nrow=I, ncol=Q) # response correlation terms
  Sigma_omega <- diag(1, Q) # correlation matrix
  Sigma_omega[1,2] <- Sigma_omega[2,1] <- 0.75
  for (i in 1:I){
    omega[i,] <- rmvn_rcpp(1, rep(0, Q), Sigma_omega)
  }
  lambda_vec = rep(2,Q)
  # define truth tau_vec
  tau_vec = rep(1,Q)
  # Generate the theta_iq
  theta_iq = array(0, dim = c(I,Q,J_max));
  for (i in 1:I) {
    for (q in 1:Q) {
      lambda_q = lambda_vec[q]
      tau_q = tau_vec[q]
      data_index_iq = which(data_index[i,q,] == 1)
      J_iq = length(data_index_iq)
      t_iq = t[i, q, data_index_iq]
      # Use vectorized operations to create the covariance matrix
      diffs = outer(t_iq, t_iq, "-")
      cov_mat = exp(- (diffs / lambda_q)^2) * tau_q
      theta_cur_iq = rmvn_rcpp(1,rep(0,J_iq),cov_mat);
      theta_iq[i,q, data_index_iq] = theta_cur_iq;
    }
  }
  
  sigma2 <- rep(0.5, Q) # variance of i.i.d Gaussian error
  beta <- array(NA, dim=c(K, Q, L))
  beta[1,1, ] = c(-0.07142442,-0.9023911,-1.6142174,-1.1723674,0.2139180,
                  1.19694463,0.6818716,-1.050855,-2.055200)
  beta[1,2, ] = c(-1.35287953,-1.4701541,-0.6735524,0.8495789,0.9193538,
                  0.03852217,0.6572662,1.492269,2.560624)
  beta[2,1, ] = 0
  beta[2,2, ] = 0
  beta[3,1, ] = c(-1.8003245,0.8970589,2.951867,3.166909,1.819671,
                  0.36735707,0.7319628,2.219103,3.484856)
  beta[3,2, ] = c(-0.3446178,1.0306363,2.302487,2.572165,1.649515,
                  0.05354716,-1.5059283,-1.963917,-2.155285)
  #
  alpha <- array(NA, dim=c(Q,S))
  alpha[,1] = c(2.36,2.74)
  alpha[,2] = c(2.12,2.13)
  # beta_truth = beta
  # alpha_truth = alpha 
  # Generate the y: simulated data
  y <- array(NA, dim = c(I, Q, J_max))  # Initialize y with same dimensions as t
  for (i in 1:I) {
    k = g[i];
    for (q in 1:Q) {
      data_index_iq <- which(data_index[i, q, ] == 1);  # Get the indices where data_index == 1
      J_iq = length(data_index_iq);
      if (length(data_index_iq) > 0) {  # Check if there are any such indices
        t_values <- t[i,q, data_index_iq];  # Get corresponding t values
        if (k == K){
          # if i-th subject in baseline class
          mean <- Z[i,q,data_index_iq,]%*%alpha[q,] + B[i,q,data_index_iq,]%*%beta[K,q,] +
            omega[i,q] + theta_iq[i,q,data_index_iq];
        }else{
          # if i-th subject not in baseline class
          mean <- Z[i,q,data_index_iq,]%*%alpha[q,] + B[i,q,data_index_iq,]%*%beta[K,q,] +
            B[i,q,data_index_iq,]%*%beta[k,q,] + omega[i,q] + theta_iq[i,q,data_index_iq];
        }
        var <- diag(sigma2[q], length(data_index_iq));
        y[i, q, data_index_iq] <- rmvn_rcpp(1, mean, var);
      }
    }
  }
  y[is.na(y)] = 0
  list_of_dfs <- list()
  # Loop over i
  for (i in 1:dim(y)[1]) {
    # Initialize the response matrix for each i, assuming max(data_indices) is the maximum possible index
    data_indices <- which(data_index[i, 1, ] == 1) # or any known maximum
    response_q <- matrix(NA, nrow =length(data_indices) , ncol =Q )
    # Loop over q (assuming q = 1, 2, 3)
    for (q in 1:Q) {
      # Get data indices where data_index[i, q, ] is 1
      data_indices <- which(data_index[i, q, ] == 1)
      # If data_indices is empty, continue to the next q
      if(length(data_indices) == 0) next
      # Get the response for these data indices
      response_q[, q] <- y[i, q, data_indices]
    }
    # If all rows in response_q are NA, continue to the next i
    if(all(is.na(response_q))) next
    # Create the data frame temp_df_i
    temp_df_i <- data.frame(ID = rep(i, length(data_indices)),
                            VISIT = data_indices,
                            time = t[i,1, data_indices],
                            cov1 = Z[i, q, data_indices, 2],
                            Group = rep(g[i], length(data_indices)))
    # Add the response_q columns to temp_df_i
    colnames(response_q) <- paste0("R", 1:Q)
    temp_df_i <- cbind(temp_df_i, response_q)
    # Add temp_df_i to the list of data frames
    list_of_dfs[[length(list_of_dfs) + 1]] <- temp_df_i
  }
  # Combine all data frames in the list into a single data frame
  final_df <- do.call(rbind, list_of_dfs)

  # Return the simulated data
  return(final_df)
}





