#' Model Selection for BayTetra
#'
#' Function implements the model selection of BayTetra for the degree of freedom in B-splines.
#'
#'
#' @param df_min Minimum value for degrees of freedom (df). Must be 3 or greater. See `mcmc_BayTetra` documentation.
#' @param df_max Maximum value for degrees of freedom (df).
#' @param data See `mcmc_BayTetra` documentation.
#' @param v_rsp See `mcmc_BayTetra` documentation.
#' @param v_covs See `mcmc_BayTetra` documentation.
#' @param v_grp See `mcmc_BayTetra` documentation.
#' @param v_time See `mcmc_BayTetra` documentation.
#' @param mcmc See `mcmc_BayTetra` documentation.
#' @param prior See `mcmc_BayTetra` documentation.
#' @param display_process See `mcmc_BayTetra` documentation.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item{`selection_elpd`: A list containing a named vector of elpd values,
#'    which elements' names are the degrees of freedom. It also contains the optimal degree of freedom.}
#'   \item{`best_mcmc_result`: The best model's MCMC posterior samples (model with the highest elpd).}
#' }
#'
#' @examples
#' \dontrun{
#' selection_result = Model_selection_BayTetra(
#'                          df_min = 4, df_max = 6,
#'                          data = ex_data,
#'                          v_rsp = paste("R", 1:2, sep = ""),
#'                          v_covs = "cov1",
#'                          v_grp = "Group",
#'                          v_time = "time"
#'                          )
#' }
#'
#' @seealso [mcmc_BayTetra()]
#' @importFrom utils modifyList
#' @importFrom loo loo
#' @export


Model_selection_BayTetra = function(df_min,df_max,
                                  data,
                                  v_rsp,
                                  v_covs,
                                  v_grp,
                                  v_time,
                                  mcmc = list(),
                                  prior = list(),
                                  display_process = TRUE
                                  ) {
  Response_vec = v_rsp;
  total_columns = length(v_covs);
  S = as.integer(total_columns + 1) # add 1 column for intercept


  responses = v_rsp
  J_max <- max(data$VISIT, na.rm = TRUE)
  I <- length(unique(data$ID))
  Q <- length(responses)

  Time_name = v_time


  data_index <- array(0, dim = c(I, Q, J_max))

  indice_problem <- list()  # Initialize list to store problematic indices

  for(q in 1:Q) {
    target <- responses[q]

    for(i in unique(data$ID)) {
      # Find rows with the current ID and non-NA target
      rows_with_values <- data[data$ID == i & !is.na(data[[target]]), ]

      # Extract the VISIT column from these rows
      visits_with_values <- rows_with_values$VISIT

      # Update data_index
      data_index[i, q, visits_with_values] <- 1

      if (all(data_index[i, q, ] == 0)) {
        # If all elements are 0, add the indices to the list
        indice_problem <- c(indice_problem, list(list(i = i, q = q)))
      }
    }
  }

  # Check if there are problematic indices
  if (length(indice_problem) > 0) {
    # Stop the execution and report the problematic indices
    stop("Found problematic indices: ", toString(indice_problem))
  }



  # Initialize array to store t
  t <- array(NA, dim = c(I, Q, J_max))

  # Loop over all IDs
  for(i in 1:I) {
    # Filter the data for the current ID
    data_i <- data[data$ID == i,]

    for(q in 1:Q) {
      # Target response column for this q
      target <- responses[q]

      # Filter data_i based on non-NA values for this specific q
      data_i_q <- data_i[!is.na(data_i[[target]]), ]

      # Get the VISIT times and corresponding AGE (or Time_name in this case)
      visit_times_q <- data_i_q$VISIT
      age_at_visit_q <- data_i_q[[Time_name]]

      # Update the t array
      t[i, q, visit_times_q] <- age_at_visit_q
    }
  }

  t_min = min(t,na.rm = TRUE)
  t_max = max(t,na.rm = TRUE)




  y <- array(0, dim = c(I, Q, J_max))

  # Loop through each ID (i)
  for (i in unique(data$ID)) {
    # Filter data for the current ID
    data_i <- data[data$ID == i, ]
    # Loop through each target variable (q)
    for (q in 1:Q) {
      data_index_iq = which(data_index[i,q,] == 1)
      target_values <- data_i[[responses[q]]]
      target_values <- target_values[data_index_iq]
      # Store these values in the array
      y[i, q, data_index_iq] <- target_values
    }
  }

  y[is.na(y)] <- 0


  total_columns = length(v_covs)
  S = as.integer(total_columns + 1) # add 1 column for intercept
  Z <- array(NA, dim = c(I, Q, J_max, S))
  for(i in 1:I){
    # Filter the data for the current ID
    # data_i <- subset(data, ID == i)
    data_i = data[data[["ID"]] == i, ]
    for(q in 1:Q){
      # Extract the indices of the non-zero elements in data_index[i,q,]
      data_index_iq <- which(data_index[i,q,] == 1)
      # First dimension of Z for intercept
      Z[i,q,,1] <- matrix(1, nrow=1, ncol=J_max)
      for(j in data_index_iq){
        # Select the row corresponding to the visit and the covariate columns
        selected_data <- data_i[data_i$VISIT == j, v_covs]
        selected_data_vector <- as.vector(t(selected_data))
        # Store this row in Z
        Z[i,q,j,2:S] <- selected_data_vector
      }
    }
  }

  first_visits <- do.call(rbind, lapply(split(data, data$ID), function(df) df[which.min(df$VISIT), ]))
  g = first_visits[[v_grp]]






  default_simulation_params <- list(
    Nit = 4000,
    burn_in = 2000,
    thin_factor = 10
  )

  # Update the default values with user-provided values
  final_simulation_params <- modifyList(default_simulation_params, mcmc)

  # # Extract variables for use
  Nit <- final_simulation_params$Nit
  burn_in <- final_simulation_params$burn_in
  thin_factor <- final_simulation_params$thin_factor
  #

  # Default hyperparameters
  default_hyper_params <- list(
    mu_alpha = rep(0, S),
    V_alpha = diag(100, S),
    nu_0 = 2.5e-5,
    a_nu = 5,
    b_nu = 25,
    a_rho = 1,
    b_rho = 1,
    h_1 = 1,
    h_2 = 1
  )

  # Update the default values with user-provided values
  final_hyper_params <- modifyList(default_hyper_params, prior)


  elpd_list = list()
  if(df_min < 4){
    stop("df_min should be 4 or greater.")
  }
  highest_elpd <- -Inf        # Initialize variable to track highest elpd
  best_mcmc_result <- NULL    # Initialize variable to store best mcmc_result

  for (df in df_min:df_max) {
    print(paste0("Current df:", df))


    num_intervals <- df - 3
    quantiles <- quantile(t, probs = seq(0, 1, 1/num_intervals), na.rm = TRUE)
    middle_quantiles <- quantiles[-c(1, num_intervals + 1)]

    if(df < 3) {
      stop("df must be at least 3.")
    }
    if(df == 3) {
      suppressWarnings({
        spline_basis <- bs(seq(t_min, t_max, length.out=1000), degree = 2, intercept = TRUE)
      })
    } else if(df == 4){
      suppressWarnings({
        spline_basis <- bs(seq(t_min, t_max, length.out=1000), df = 3, intercept = TRUE)
      })
    } else {
      spline_basis <- bs(seq(t_min, t_max, length.out=1000), knots = middle_quantiles, intercept = TRUE)
    }
    L <- dim(spline_basis)[2]  # degree of freedoms for spline expansion with intercept
    B <- array(NA, dim=c(I, Q, J_max, L)) # time spline basis matrix
    for (i in 1:I){
      for (q in 1:Q){
        data_index_iq <- which(data_index[i,q,] == 1) # index of data for i-th subject at q-th response
        B[i,q,data_index_iq,] <- cbind( predict(spline_basis, t[i,q,data_index_iq]))
      }
    }

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
    B_new = array(NA, dim=c(I, Q, J_max, L))
    for (i in 1:I){
      for (q in 1:Q){
        data_index_iq <- which(data_index[i,q,] == 1) # index of data for i-th subject at q-th response
        if (length(data_index_iq) == 1) {
          B_new[i,q,data_index_iq,] = c(1, B_after_qr_allq[i,q, data_index_iq , ])
        } else {
          B_new[i,q,data_index_iq,] = cbind(rep(1,length(data_index_iq)), B_after_qr_allq[i,q, data_index_iq , ])
        }
      }
    }
    B = B_new
    B[is.na(B)] <- 0




    mcmc_result <- mcmc_BayTetra(data = data,
                              v_rsp = Response_vec,
                              v_covs = v_covs,
                              v_grp = v_grp,
                              v_time = v_time,
                              df = df,
                              prior = final_hyper_params,
                              mcmc = final_simulation_params,
                              display_process = display_process)


    post_alpha_samples = mcmc_result$pos.alpha
    post_beta_samples = mcmc_result$pos.beta
    post_Sigma_omega_samples = mcmc_result$pos.Sigma_omega
    post_sigma_samples = mcmc_result$pos.sigma2
    post_gamma_kq_samples = mcmc_result$pos.gamma_kq
    post_gamma_kq0_samples = mcmc_result$pos.gamma_kq0
    # B = mcmc_result$B
    # y = mcmc_result$y
    # g = mcmc_result$g
    # Z = mcmc_result$Z
    # data_index = mcmc_result$data_index

    beta_result = test_beta(post_beta_samples)
    beta_estimate = beta_result$beta_estimation
    beta_test = beta_result$beta_test_result


    signal_all = extract_signal(post_gamma_kq_samples,
                                post_gamma_kq0_samples)
    signal_kq_real = signal_all$signal_kq_real
    signal_kq0_real = signal_all$signal_kq0_real
    signal_kq = signal_all$signal_kq
    signal_kq0 = signal_all$signal_kq0

    beta_update_result = update_beta_based_on_signal(post_beta_samples,beta_estimate,
                                                     signal_kq,signal_kq0)
    new_beta_samples = beta_update_result$new_beta_samples
    new_beta_estimate = beta_update_result$new_beta_estimate
    # estimation_alpha = estimate_alpha(post_alpha_samples)

    post_sample_length = (Nit-burn_in)/thin_factor


    I <- dim(data_index)[1]
    Q <- dim(data_index)[2]
    J_max = dim(y)[3]
    K = dim(post_beta_samples)[2]
    L = dim(post_beta_samples)[4]


    total_records <- dim(data)[1]


    # Initialize model likelihood matrix
    model_llh_matrix = array(0, dim = c(total_records, post_sample_length))

    for (c in 1:post_sample_length) {
      # c = 1
      current_Sigma_omega = post_Sigma_omega_samples[c, , ]
      current_alpha = post_alpha_samples[c, , ]
      current_beta = new_beta_samples[c, , , ]
      current_sigma2 = post_sigma_samples[c, ]

      # Initialize a vector to store the log likelihood values
      sum_ll = vector("numeric", length = total_records)

      # Loop counters for filling up the sum_ll vector
      counter = 1

      for (i in 1:I) {
        k = g[i]
        data_index_iq = which(data_index[i,1,] == 1)
        for (j in data_index_iq) {
          # Initialize vectors for y_ij and mean
          y_ij = vector("numeric", length = Q)
          mean = vector("numeric", length = Q)

          for (q in 1:Q) {
            y_ij[q] = y[i,q,j]
            if (k == K) {
              mean[q] = t(Z[i,q,j, ])%*%current_alpha[q,] +
                t(B[i,q,j, ])%*%current_beta[K,q,]
            } else {
              mean[q] = t(Z[i,q,j, ])%*%current_alpha[q,] +
                t(B[i,q,j, ])%*%current_beta[K,q,] +
                t(B[i,q,j, ])%*%current_beta[k,q,]
            }
          }

          # Calculate sigma1, sigma2, and sigma_all
          sigma1 = diag(current_sigma2)
          sigma2 = current_Sigma_omega
          sigma_all = sigma1 + sigma2
          # print(sigma_all)
          # Compute the log likelihood and store in sum_ll vector
          sum_ll[counter] = dmvn_rcpp(y_ij, mean, sigma_all, logd = TRUE)

          # Increment counter
          counter = counter + 1
        }
      }

      model_llh_matrix[, c] <- sum_ll
    }

    model_llh_matrix = t(model_llh_matrix)

    suppressWarnings({
      loo_result = loo(model_llh_matrix)
    })
    current_elpd <- loo_result$estimates[1, 1]
    elpd_list = append(elpd_list, current_elpd)

    if (current_elpd > highest_elpd) {
      highest_elpd = current_elpd
      highest_df <- df
      best_mcmc_result <- mcmc_result
    }

  }
  selection_elpd <- unlist(lapply(elpd_list, function(x) x[[1]]))
  names(selection_elpd) <- paste0("df", df_min:df_max)

  selection_info <- list(
    elpd = selection_elpd,
    max_elpd_info = paste("The maximum elpd is reached at df =", highest_df)
  )

  selection_elpd = selection_info
  # class(selection_elpd) <- "Select_BayTetra"

  # # Attach the highest_df as an attribute
  # attributes(selection_elpd) <- list(highest_df = highest_df)

  return(list(selection_elpd,best_mcmc_result) )
}




