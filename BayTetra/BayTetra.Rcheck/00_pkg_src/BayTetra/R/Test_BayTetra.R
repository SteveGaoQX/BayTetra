#' Hypothesis Testing for BayTetra
#'
#' This function implements hypothesis test based on BayTetra posterior samples.
#'
#' @param object An object of class "Post_BayTetra" containing MCMC posterior samples from "mcmc_BayTetra".
#' @param v_rsp A character vector of response variables.
#'
#' @return Return an object of class "Test_BayTetra" containing three elements:
#' \itemize{
#'   \item{BayTetra_summary}: {A dataframe that summarizes the testing information for \eqn{\widetilde{\beta}_{k q 0}} and \eqn{\widetilde{ \boldsymbol{\beta } }_{k q}^{-}}. }
#'   \item{pairwise_significance}: {A matrix indicating pairwise hypothesis testing
#'    results. For each pair of groups, their longitudinal trajectories are significantly
#'    different if the corresponding responses present in the table.}
#'   \item{diff_among_all_grps}: {A named logical vector indicating if there exists difference among all groups for each response.}
#' }
#'
#' @examples
#' \dontrun{
#' mcmc_result = mcmc_BayTetra(data = ex_data,
#'                             v_rsp = paste("R", 1:2,sep = ""),
#'                             v_covs = "cov1",
#'                             v_grp = "Group",
#'                             v_time = "time",
#'                             df = 10)
#' test_result <- Test_BayTetra(mcmc_result,
#'                            v_rsp = paste("R", 1:2,sep = ""))
#' }
#' @importFrom stats quantile sd
#'
#' @export


Test_BayTetra = function(object,
                         v_rsp){
  # v_rsp = paste("R", 1:2,sep = "")
  # mcmc_posterior = mcmc_result

  if(length(v_rsp)>1){

    post_alpha_samples = object$pos.alpha
    post_beta_samples = object$pos.beta
    post_beta_kq0_samples = object$pos.beta_kq0
    post_Sigma_omega_samples = object$pos.Sigma_omega
    post_sigma_samples = object$pos.sigma2


    num_samples <- dim(post_beta_samples)[1]
    K <- dim(post_beta_samples)[2]
    Q <- dim(post_beta_samples)[3]
    L <- dim(post_beta_samples)[4]



    beta_result = test_beta(post_beta_samples)
    beta_estimate = beta_result$beta_estimation
    beta_test = beta_result$beta_test_result
    
    beta_kq0_result = test_beta_kq0(post_beta_kq0_samples)
    beta_kq0_estimate = beta_kq0_result$beta_estimation
    beta_kq0_test = beta_kq0_result$beta_test_result
    

    # Initialize the beta_upper and beta_lower matrices
    beta_upper <- array(0, dim = c(K, Q, L))
    beta_lower <- array(0, dim = c(K, Q, L))
    beta_std = array(0, dim =c(K,Q,L))
    # Calculate the 90% and 10% quantiles
    for(k in 1:K){
      for(q in 1:Q){
        for(l in 1:L){
          beta_upper[k, q, l] <- quantile(post_beta_samples[, k, q, l], 0.975)
          beta_lower[k, q, l] <- quantile(post_beta_samples[, k, q, l], 0.025)
          beta_std[k,q,l] = sd(post_beta_samples[,k,q,l])
        }
      }
    }
    
    # Initialize the beta_upper and beta_lower matrices
    beta_kq0_upper <- array(0, dim = c(K, Q))
    beta_kq0_lower <- array(0, dim = c(K, Q))
    beta_kq0_std = array(0, dim =c(K,Q))
    # Calculate the 90% and 10% quantiles
    for(k in 1:K){
      for(q in 1:Q){
        beta_kq0_upper[k, q] <- quantile(post_beta_kq0_samples[,k,q], 0.975)
        beta_kq0_lower[k, q] <- quantile(post_beta_kq0_samples[,k,q], 0.025)
        beta_kq0_std[k,q] = sd(post_beta_kq0_samples[,k,q])
      }
    }

    # Initialize an empty data frame
    summary_table <- data.frame(
      row_name = character(),
      Posterior_Mean = character(),
      Std = character(),
      Lower_Bound = character(),
      Upper_Bound = character(),
      stringsAsFactors = FALSE
    )

    # Loop through k and q to fill the table
    for (q in 1:Q) {
      for (k in 1:K) {
        # Construct row names
        row_name_beta_kq <- paste0("beta_", k, q)
        row_name_beta_kq0 <- paste0("beta_", k, q, "_0")

        # Get current values
        beta_current <- beta_estimate[k, q, ]
        beta_current_std <- beta_std[k, q, ]
        beta_current_upper <- beta_upper[k, q, ]
        beta_current_lower <- beta_lower[k, q, ]
        
        beta_kq0_current <- beta_kq0_estimate[k, q]
        beta_kq0_current_std <- beta_kq0_std[k, q]
        beta_kq0_current_upper <- beta_kq0_upper[k, q]
        beta_kq0_current_lower <- beta_kq0_lower[k, q]


        # Calculate Mean, std, lower and upper bound for beta_kq and beta_kq0
        mean_beta_kq <- paste("(", paste(round(beta_current, 2), collapse = ", "), ")", sep = "")
        mean_beta_kq0 <- as.character(round(beta_kq0_current, 2))

        std_beta_kq <- paste("(", paste(round(beta_current_std, 2), collapse = ", "), ")", sep = "")
        std_beta_kq0 <- as.character(round(beta_kq0_current_std, 2))

        upper_beta_kq <- paste("(", paste(round(beta_current_upper, 2), collapse = ", "), ")", sep = "")
        upper_beta_kq0 <- as.character(round(beta_kq0_current_upper, 2))

        lower_beta_kq <- paste("(", paste(round(beta_current_lower, 2), collapse = ", "), ")", sep = "")
        lower_beta_kq0 <- as.character(round(beta_kq0_current_lower, 2))
        
        
        # mean_beta_kq <- paste("(", paste(round(beta_current, 2), collapse = ", "), ")", sep = "")
        # mean_beta_kq0 <- as.character(round(beta_kq0_estimate, 2))
        # std_beta_kq <- paste("(", paste(round(beta_current_std[-1], 2), collapse = ", "), ")", sep = "")
        # std_beta_kq0 <- as.character(round(beta_current_std[1], 2))
        # upper_beta_kq <- paste("(", paste(round(beta_current_upper[-1], 2), collapse = ", "), ")", sep = "")
        # upper_beta_kq0 <- as.character(round(beta_current_upper[1], 2))
        # lower_beta_kq <- paste("(", paste(round(beta_current_lower[-1], 2), collapse = ", "), ")", sep = "")
        # lower_beta_kq0 <- as.character(round(beta_current_lower[1], 2))

        # Append to the data frame
        summary_table <- rbind(
          summary_table,
          data.frame(
            row_name = row_name_beta_kq,
            Posterior_Mean = mean_beta_kq,
            Std = std_beta_kq,
            Lower_Bound = lower_beta_kq,
            Upper_Bound = upper_beta_kq,
            stringsAsFactors = FALSE
          ),
          data.frame(
            row_name = row_name_beta_kq0,
            Posterior_Mean = mean_beta_kq0,
            Std = std_beta_kq0,
            Lower_Bound = lower_beta_kq0,
            Upper_Bound = upper_beta_kq0,
            stringsAsFactors = FALSE
          )
        )
      }
    }

    # Show the final table
    rownames(summary_table) <- summary_table$row_name

    # Remove the "row_name" column as it's now redundant
    summary_table$row_name <- NULL



    # make a copy of the original beta samples
    new_beta_samples = post_beta_samples
    new_beta_estimate = beta_estimate


    col_group_names <- row_group_names <- paste0("group", 1:K)
    hypo_matrix <- matrix(NA, ncol=K, nrow=K)
    colnames(hypo_matrix) <- col_group_names
    rownames(hypo_matrix) <- row_group_names

    responses = v_rsp

    tests = list()
    # Loop over the Q groups
    for (q in 1:Q) {
      current_response = responses[q]
      record_matrix <- matrix(NA, ncol=K, nrow=K)
      colnames(record_matrix) <- col_group_names
      rownames(record_matrix) <- row_group_names
      # Loop over the groups
      for (k1 in 1:(K-1)) {
        for (k2 in (k1+1):K) {  # Only start k2 from k1+1 to avoid duplicates

          if (k1 != k2) {  # Skip comparison of the same group
            if (k2 == K) {
              diff_temp = new_beta_samples[, k1, q, ]
            } else {
              diff_temp <- new_beta_samples[, k1, q, ] - new_beta_samples[, k2, q, ]
            }

            # Calculate the confidence intervals
            ci_lower <- sapply(1:L, function(l) quantile(diff_temp[, l], 0.025))
            ci_upper <- sapply(1:L, function(l) quantile(diff_temp[, l], 0.975))

            # Perform the hypothesis test
            test_temp <- !(0 >= ci_lower & 0 <= ci_upper)

            # Update the hypothesis matrix based on the test result
            if (any(test_temp == TRUE)) {
              if (is.na(hypo_matrix[k1, k2]) ){
                hypo_matrix[k1, k2] <- current_response
              } else {
                hypo_matrix[k1, k2] <- paste(hypo_matrix[k1, k2], current_response, sep=", ")
              }
              record_matrix[k1, k2] = 1

            }
          }
        }
      }

      tests[[q]] = record_matrix
    }

    # Post-process each hypo_matrix
    for (i in 1:K) {
      for (j in 1:K) {
        if (i == j || j < i) {
          # Set diagonal and lower triangular to NA
          hypo_matrix[i, j] <- NA
        } else {
          if (is.na(hypo_matrix[i, j])) {
            hypo_matrix[i, j] <- "None"
          } else {
            hypo_matrix[i, j] <- paste0("(", hypo_matrix[i, j], ")")
          }
        }
      }
    }


    result <- c(list(BayTetra_summary = summary_table), list(pairwise_significance = hypo_matrix))

    # Initialize an empty named vector
    diff_all <- logical(0)
    names(diff_all) <- character(0)

    # Loop over the Q groups
    for (q in 1:Q) {
      # Retrieve the matrix for the q-th group
      current_matrix <- tests[[q]]

      # Retrieve the response name
      current_response <- responses[q]

      # Check if any value in the matrix is 1
      if (any(current_matrix == 1)) {
        diff_all[current_response] <- TRUE
      } else {
        diff_all[current_response] <- FALSE
      }
    }

    result$diff_among_all_grps <- diff_all


    class(result) = "Test_BayTetra"
    return(result)

  }else{
      stop(print("Q =1 Case is not supported"))
  }





}





