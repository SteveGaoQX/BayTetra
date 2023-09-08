#' Hypothesis Testing for BayTetra
#'
#' This function implements hypothesis test based on BayTetra posterior samples.
#'
#' @param object An object of class "Post_BayTetra" containing MCMC posterior samples from "mcmc_BayTetra".
#' @param v_rsp A character vector of response variables.
#'
#' @return Return an object of class "Test_BayTetra" containing three elements:
#' \itemize{
#'   \item{BayTetra_summary}: {A dataframe that summarizes the testing information for \eqn{\widetilde{\beta}_{k q 0}} and \eqn{\widetilde{ \boldsymbol{\beta } }_{k q}^{-}}. The 'Pr(Signal)' column represents the confidence level \eqn{P(\gamma_{kq} = 1)}.}
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
#'                             df = 5)
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
    post_Sigma_omega_samples = object$pos.Sigma_omega
    post_sigma_samples = object$pos.sigma2
    post_gamma_kq_samples = object$pos.gamma_kq
    post_gamma_kq0_samples = object$pos.gamma_kq0

    num_samples <- dim(post_beta_samples)[1]
    K <- dim(post_beta_samples)[2]
    Q <- dim(post_beta_samples)[3]
    L <- dim(post_beta_samples)[4]



    beta_result = test_beta(post_beta_samples)
    beta_estimate = beta_result$beta_estimation
    beta_test = beta_result$beta_test_result


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

    # Initialize an empty data frame
    summary_table <- data.frame(
      row_name = character(),
      Posterior_Mean = character(),
      Std = character(),
      Lower_Bound = character(),
      Upper_Bound = character(),
      Pr_Signal = numeric(),
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

        signal_current_kq <- signal_kq_real[k, q]
        signal_current_kq0 <- signal_kq0_real[k, q]

        # Calculate Mean, std, lower and upper bound for beta_kq and beta_kq0
        mean_beta_kq <- paste("(", paste(round(beta_current[-1], 2), collapse = ", "), ")", sep = "")
        mean_beta_kq0 <- as.character(round(beta_current[1], 2))

        std_beta_kq <- paste("(", paste(round(beta_current_std[-1], 2), collapse = ", "), ")", sep = "")
        std_beta_kq0 <- as.character(round(beta_current_std[1], 2))

        upper_beta_kq <- paste("(", paste(round(beta_current_upper[-1], 2), collapse = ", "), ")", sep = "")
        upper_beta_kq0 <- as.character(round(beta_current_upper[1], 2))

        lower_beta_kq <- paste("(", paste(round(beta_current_lower[-1], 2), collapse = ", "), ")", sep = "")
        lower_beta_kq0 <- as.character(round(beta_current_lower[1], 2))

        # Append to the data frame
        summary_table <- rbind(
          summary_table,
          data.frame(
            row_name = row_name_beta_kq,
            Posterior_Mean = mean_beta_kq,
            Std = std_beta_kq,
            Lower_Bound = lower_beta_kq,
            Upper_Bound = upper_beta_kq,
            Pr_Signal = signal_current_kq,
            stringsAsFactors = FALSE
          ),
          data.frame(
            row_name = row_name_beta_kq0,
            Posterior_Mean = mean_beta_kq0,
            Std = std_beta_kq0,
            Lower_Bound = lower_beta_kq0,
            Upper_Bound = upper_beta_kq0,
            Pr_Signal = signal_current_kq0,
            stringsAsFactors = FALSE
          )
        )
      }
    }

    # Show the final table
    rownames(summary_table) <- summary_table$row_name

    # Remove the "row_name" column as it's now redundant
    summary_table$row_name <- NULL
    colnames(summary_table)[colnames(summary_table) == "Pr_Signal"] <- "Pr(Signal)"



    # make a copy of the original beta samples
    new_beta_samples = post_beta_samples
    new_beta_estimate = beta_estimate


    # loop over K and Q
    for (k in 1:K) {
      for (q in 1:Q) {
        # if signal_kq[k,q] is 0, set new_beta_samples[, k, q, 2:L] to 0
        if (signal_kq[k,q] == 0) {
          new_beta_samples[, k, q, 2:L] <- 0
          new_beta_estimate[k,q,2:L] = 0
        }

        # if signal_kq0[k,q] is 0, set new_beta_samples[, k, q, 1] to 0
        if (signal_kq0[k,q] == 0) {
          new_beta_samples[, k, q, 1] <- 0
          new_beta_estimate[k,q,1] = 0
        }
      }
    }



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
    post_alpha_samples = object$pos.alpha
    post_beta_samples = object$pos.beta
    post_sigma_samples = object$pos.sigma2
    post_gamma_kq_samples = object$pos.gamma_kq
    post_gamma_kq0_samples = object$pos.gamma_kq0

    num_samples <- dim(post_beta_samples)[1]
    K <- dim(post_beta_samples)[2]
    L <- dim(post_beta_samples)[3]


    beta_result = test_beta_Q1(post_beta_samples)
    beta_estimate = beta_result$beta_estimation
    beta_test = beta_result$beta_test_result



    signal_kq <- rep(0,K)
    signal_kq0 <- rep(0,K)
    signal_kq_real <- rep(0,K)
    signal_kq0_real <- rep(0,K)

    for (k in 1:K) {
      prob_kq <- mean(post_gamma_kq_samples[ , k] == 1)
      prob_kq0 <- mean(post_gamma_kq0_samples[ , k] == 1)
      signal_kq_real[k] = prob_kq
      signal_kq0_real[k] = prob_kq0
      signal_kq[k] <- ifelse(prob_kq > 0.5, 1, 0)
      signal_kq0[k] <- ifelse(prob_kq0 > 0.5, 1, 0)
    }

    # Initialize the beta_upper and beta_lower matrices
    beta_upper <- array(0, dim = c(K, L))
    beta_lower <- array(0, dim = c(K, L))
    beta_std = array(0, dim =c(K, L))
    # Calculate the 90% and 10% quantiles
    for(k in 1:K){

      for(l in 1:L){
        beta_upper[k, l] <- quantile(post_beta_samples[, k, l], 0.975)
        beta_lower[k, l] <- quantile(post_beta_samples[, k, l], 0.025)
        beta_std[k,l] = sd(post_beta_samples[,k,l])
      }

    }

    # Initialize an empty data frame
    summary_table <- data.frame(
      row_name = character(),
      Posterior_Mean = character(),
      Std = character(),
      Lower_Bound = character(),
      Upper_Bound = character(),
      Pr_Signal = numeric(),
      stringsAsFactors = FALSE
    )

    # Loop through k and q to fill the table

    for (k in 1:K) {
      # Construct row names
      row_name_beta_kq <- paste0("beta_", k)
      row_name_beta_kq0 <- paste0("beta_", k, "_0")

      # Get current values
      beta_current <- beta_estimate[k, ]
      beta_current_std <- beta_std[k, ]
      beta_current_upper <- beta_upper[k, ]
      beta_current_lower <- beta_lower[k, ]

      signal_current_kq <- signal_kq_real[k]
      signal_current_kq0 <- signal_kq0_real[k]

      # Calculate Mean, std, lower and upper bound for beta_kq and beta_kq0
      mean_beta_kq <- paste("(", paste(round(beta_current[-1], 2), collapse = ", "), ")", sep = "")
      mean_beta_kq0 <- as.character(round(beta_current[1], 2))

      std_beta_kq <- paste("(", paste(round(beta_current_std[-1], 2), collapse = ", "), ")", sep = "")
      std_beta_kq0 <- as.character(round(beta_current_std[1], 2))

      upper_beta_kq <- paste("(", paste(round(beta_current_upper[-1], 2), collapse = ", "), ")", sep = "")
      upper_beta_kq0 <- as.character(round(beta_current_upper[1], 2))

      lower_beta_kq <- paste("(", paste(round(beta_current_lower[-1], 2), collapse = ", "), ")", sep = "")
      lower_beta_kq0 <- as.character(round(beta_current_lower[1], 2))

      # Append to the data frame
      summary_table <- rbind(
        summary_table,
        data.frame(
          row_name = row_name_beta_kq,
          Posterior_Mean = mean_beta_kq,
          Std = std_beta_kq,
          Lower_Bound = lower_beta_kq,
          Upper_Bound = upper_beta_kq,
          Pr_Signal = signal_current_kq,
          stringsAsFactors = FALSE
        ),
        data.frame(
          row_name = row_name_beta_kq0,
          Posterior_Mean = mean_beta_kq0,
          Std = std_beta_kq0,
          Lower_Bound = lower_beta_kq0,
          Upper_Bound = upper_beta_kq0,
          Pr_Signal = signal_current_kq0,
          stringsAsFactors = FALSE
        )
      )
    }


    # Show the final table
    rownames(summary_table) <- summary_table$row_name

    # Remove the "row_name" column as it's now redundant
    summary_table$row_name <- NULL
    colnames(summary_table)[colnames(summary_table) == "Pr_Signal"] <- "Pr(Signal)"


    # second part group difference

    new_beta_samples = post_beta_samples
    new_beta_estimate = beta_estimate
    # get the dimensions
    K <- dim(new_beta_samples)[2]
    L <- dim(new_beta_samples)[3]

    # loop over K and Q
    for (k in 1:K) {

      # if signal_kq[k,q] is 0, set new_beta_samples[, k, q, 2:L] to 0
      if (signal_kq[k] == 0) {
        new_beta_samples[, k,  2:L] <- 0
        new_beta_estimate[k,2:L] = 0
      }

      # if signal_kq0[k,q] is 0, set new_beta_samples[, k, q, 1] to 0
      if (signal_kq0[k] == 0) {
        new_beta_samples[, k, 1] <- 0
        new_beta_estimate[k,1] = 0
      }

    }

    # start to create the test matrix
    col_group_names <- row_group_names <- paste0("group", 1:K)
    hypo_matrix <- matrix(NA, ncol=K, nrow=K)
    colnames(hypo_matrix) <- col_group_names
    rownames(hypo_matrix) <- row_group_names

    responses = v_rsp

    tests = NULL
    # Loop over the Q groups

    current_response = responses
    record_matrix <- matrix(NA, ncol=K, nrow=K)
    colnames(record_matrix) <- col_group_names
    rownames(record_matrix) <- row_group_names
    # Loop over the groups
    for (k1 in 1:(K-1)) {
      for (k2 in (k1+1):K) {  # Only start k2 from k1+1 to avoid duplicates

        if (k1 != k2) {  # Skip comparison of the same group
          if (k2 == K) {
            diff_temp = new_beta_samples[, k1,  ]
          } else {
            diff_temp <- new_beta_samples[, k1,  ] - new_beta_samples[, k2,  ]
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

    tests = record_matrix


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

    diff_all <- logical(0)
    names(diff_all) <- character(0)

    if (any(tests == 1)) {
      diff_all[current_response] <- TRUE
    } else {
      diff_all[current_response] <- FALSE
    }


    result$diff_among_all_grps <- diff_all

    class(result) = "Test_BayTetra"
    return(result)

  }





}





