#' Posterior inference for BayTetra
#'
#' Draw posterior samples of the parameters of interest from BayTetra
#'
#'
#'
#' @details
#' The model of the BayTetra is:
#' \deqn{y_{i q j}=\boldsymbol{Z}_{i j}^{\mathrm{T}} \boldsymbol{\alpha}_q+\sum_{l=1}^{L-1} \widetilde{\beta}_{K q l} \widetilde{\boldsymbol{B}}_l\left(t_{i q j}\right)+\sum_{k=1}^{K-1} \mathbb{I}\left(g_i=k\right)\left(\widetilde{\beta}_{k q 0}+\sum_{l=1}^{L-1} \widetilde{\beta}_{k q l} \widetilde{\boldsymbol{B}}_l\left(t_{i q j}\right)\right)+\omega_{i q}+ \theta_{iqj}+ \epsilon_{i q j},}
#'
#' \deqn{ \boldsymbol{\omega}_i=\left(\omega_{i 1}, \ldots, \omega_{i Q}\right) \sim \mathcal{N}\left(\mathbf{0}, \Sigma_\omega\right),\boldsymbol{\theta}_{i q}=\left(\theta_{i q 1}, \ldots, \theta_{i q, J_i}\right) \sim \mathcal{N}\left(\mathbf{0}, \boldsymbol{\Sigma}_{\boldsymbol{\theta}_{i q}}\right), \epsilon_{i q j} \sim \mathcal{N}\left(0, \sigma_q^2\right),}
#'
#'
#'
#' where \eqn{\widetilde{\boldsymbol{B}}_l\left(t_{i q j}\right)} denote the \eqn{l}-th basis function for the \eqn{L-1} dimensional
#' cubic B-spline expansion at time \eqn{t_{iqj}}, where \eqn{\boldsymbol{\Sigma}_{\boldsymbol{\omega}}} is a correlation matrix,
#' and \eqn{\boldsymbol{\Sigma}_{\boldsymbol{\theta}_{i q}}} is a \eqn{J_i \times J_i} squared exponential covariance matrix whose
#' \eqn{\left(j, j^{\prime}\right)}-th entry is \eqn{\tau_q^2 \exp \left\{-\left(\frac{t_{i q j}-t_{i q j^{\prime}}}{\lambda_q}\right)^2\right\}}.
#'
#'
#' We set \eqn{\widetilde{\beta}_{K q 0}} to 0 for identifiability and denote
#' \eqn{\widetilde{\boldsymbol{\beta}}_{k q}=\left(\widetilde{\beta}_{k q 0}, \widetilde{\beta}_{k q 1}, \ldots, \widetilde{\beta}_{k q, L-1}\right)^{\mathrm{T}}=\left(\widetilde{\beta}_{k q 0},\left(\widetilde{ \boldsymbol{\beta } }_{k q}^{-}\right)^{\mathrm{T}}\right)^{\mathrm{T}} .}
#'
#' We assign priors:
#' \deqn{\widetilde{\boldsymbol{\beta}}_{k q}^{-}=\eta_{k q} \boldsymbol{\xi}_{k q},}
#' \deqn{\eta_{k q} \sim \mathcal{N}\left(0, \gamma_{k q} \nu_{k q}\right),}
#' \deqn{\gamma_{k q} \sim \rho \delta_1\left(\gamma_{k q}\right)+(1-\rho) \delta_{\nu_0}\left(\gamma_{k q}\right),}
#' \deqn{\nu_{k q} \sim \operatorname{Inverse-Gamma}\left(a_\nu, b_\nu\right),}
#' \deqn{\rho \sim \text{Beta}(a_{\rho},b_{\rho}), \xi_{kql} \sim \mathcal{N}(m_{kql}, 1),}
#' \deqn{m_{kql} \sim \frac{1}{2}\delta_1(m_{kql}) + \frac{1}{2}\delta_{-1}(m_{kql}), l=1,2,\dots,L-1.}
#'
#'
#' For the intercept \eqn{\widetilde{\beta}_{k q 0}}, we assume its prior:
#' \deqn{\widetilde{\beta}_{k q 0} \sim \mathcal{N}\left(0, \gamma_{k q 0} \nu_{k q 0}\right),}
#' \deqn{\gamma_{k q 0} \sim \rho_0 \delta_1\left(\gamma_{k q 0}\right)+\left(1-\rho_0\right) \delta_{\nu_0}\left(\gamma_{k q 0}\right).}
#'
#'
#' The prior of other parameters are:
#' \deqn{\boldsymbol{\alpha}_q \sim \mathcal{N}\left(0, \Sigma_\alpha\right), p\left(\Sigma_\omega\right) \propto 1 ,}
#' \deqn{\tau_q^2 \sim \text { Inverse-Gamma }\left(a_\tau, b_\tau\right) ,}
#' \deqn{\lambda_q \sim \text { Inverse-Gamma }\left(a_\lambda, b_\lambda\right) ,}
#' \deqn{\sigma_q^2 \sim \text { Inverse-Gamma }\left(a_\sigma, b_\sigma\right) .}
#'
#'
#'
#' @param data longitudinal data with ID, VISIT, Group, and Covariates, Responses, Time.
#' @param v_rsp Column names corresponding to responses.
#' @param v_covs Column names corresponding to covariates.
#' @param v_grp Column name corresponding to group memberships.
#' @param v_time Column name corresponding to time.
#' @param df This parameter specifies the degree of freedom of B-spline and is used to select
#' the number of interior knots. Default value is 4 and minimum value is 3.
#' \itemize{
#'   \item {df = 3: Function uses a degree 2 B-spline with 0 interior knots.}
#'   \item {df = 4: Function uses a degree 3 B-spline with 0 interior knots.}
#'   \item {df >= 5: Function uses a degree 3 B-spline with (df - 4) interior knots.}
#' }
#' @param prior A list giving the prior information.
#'  - \code{mu_alpha}: The mean in normal prior for \eqn{\alpha_{q}}. Default value is a zero vector.
#'  - \code{V_alpha}: The covariance matrix in normal prior for \eqn{\alpha_{q}}. Default value is 100 * \eqn{I} where \eqn{I} is the identity matrix.
#'  - \code{nu_0}: The hyperparameter \eqn{\nu_0} in prior for \eqn{\gamma_{kq}} and \eqn{\gamma_{kq0}}. Default value is 2.5e-5.
#'  - \code{a_nu}: The hyperparameter \eqn{a_{\nu}} in prior for \eqn{\nu_{kq}} and \eqn{\nu_{kq0}}. Default value is 5.
#'  - \code{b_nu}: The hyperparameter \eqn{b_{\nu}} in prior for \eqn{\nu_{kq}} and \eqn{\nu_{kq0}}. Default value is 25.
#'  - \code{a_rho}: The hyperparameter \eqn{a_{\rho}} in prior for \eqn{\rho_0}. Default value is 1.
#'  - \code{b_rho}: The hyperparameter \eqn{b_{\rho}} in prior for \eqn{\rho_0}. Default value is 1.
#'  - \code{a_tau}: The hyperparameter \eqn{a_{\tau}} in prior for \eqn{\tau_q^2}. Default value is 1.
#'  - \code{b_tau}: The hyperparameter \eqn{b_{\tau}} in prior for \eqn{\tau_q^2}. Default value is 1.
#'  - \code{a_lamb}: The hyperparameter \eqn{a_{\lambda}} in prior for \eqn{\lambda_q}. Default value is 1.
#'  - \code{b_lamb}: The hyperparameter \eqn{b_{\lambda}} in prior for \eqn{\lambda_q}. Default value is 1.
#'  - \code{h_1}: The hyperparameter \eqn{a_{\sigma}} in prior for \eqn{\sigma_q^2}. Default value is 1.
#'  - \code{h_2}: The hyperparameter \eqn{b_{\sigma}} in prior for \eqn{\sigma_q^2}. Default value is 1.
#' @param mcmc A list giving the MCMC parameters.
#'  - \code{Nit}: The number of iterations for the MCMC chain. Default is 4000.
#'  - \code{burn_in}: The number of burn-in samples in the MCMC chain. Default is 2000.
#'  - \code{thin_factor}: The thinning factor for the chain. Default is 10.
#' @param display_process A bool value; if TRUE, progress will be displayed every 1000 iteration by default.
#'
#' @return An object of class 'Post_BayTetra' containing posterior samples:
#'   \itemize{
#'     \item \code{pos.alpha}: Posterior samples for \eqn{\alpha_{q}}.
#'     \item \code{pos.beta}: Posterior samples for \eqn{\widetilde{\boldsymbol{\beta}}_{k q}}.
#'     \item \code{pos.gamma_kq}: Posterior samples for \eqn{\gamma_{kq}}.
#'     \item \code{pos.gamma_kq0}: Posterior samples for \eqn{\gamma_{kq0}}.
#'     \item \code{pos.Sigma_omega}: Posterior samples for \eqn{\Sigma_{\omega}}.
#'     \item \code{pos.tau_q}: Posterior samples for \eqn{\tau^2_q}.
#'     \item \code{pos.lambda_q}: Posterior samples for \eqn{\lambda^2_q}.
#'     \item \code{pos.sigma2}: Posterior samples for \eqn{\sigma^2_q}.
#'    }

#'
#' @examples
#' \dontrun{
#' mcmc_result = mcmc_BayTetra(data = ex_data,
#'                             v_rsp = paste("R", 1:2,sep = ""),
#'                             v_covs = "cov1",
#'                             v_grp = "Group",
#'                             v_time = "time",
#'                             df = 5)
#' }
#'
#'
#'
#' @importFrom MCMCpack rinvgamma
#' @importFrom MASS mvrnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom splines bs
#' @importFrom pracma householder
#' @importFrom stats predict quantile rbeta rgamma rnorm
#' @importFrom utils modifyList

#' @export

mcmc_BayTetra <- function(data,
                          v_rsp,
                          v_covs,
                          v_grp,
                          v_time,
                          df ,
                          prior = list(),
                          mcmc = list(),
                          display_process = TRUE) {


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
      # Store these values in the array
      y[i, q, data_index_iq] <- target_values
    }
  }

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




  K = max(data[[v_grp]],na.rm = TRUE)

  default_simulation_params <- list(
    Nit = 4000,
    burn_in = 2000,
    thin_factor = 10
  )

  # Update the default values with user-provided values
  final_simulation_params <- modifyList(default_simulation_params, mcmc)

  # Extract variables for use
  Nit <- final_simulation_params$Nit
  burn_in <- final_simulation_params$burn_in
  thin_factor <- final_simulation_params$thin_factor


  # Default hyperparameters
  default_hyper_params <- list(
    mu_alpha = rep(0, S),
    V_alpha = diag(100, S),
    nu_0 = 2.5e-5,
    a_nu = 5,
    b_nu = 25,
    a_rho = 1,
    b_rho = 1,
    a_lamb = 1,
    b_lamb = 1,
    a_tau = 1,
    b_tau = 1,
    h_1 = 1,
    h_2 = 1
  )

  # Update the default values with user-provided values
  final_hyper_params <- modifyList(default_hyper_params, prior)

  # Extract variables for use
  mu_alpha <- final_hyper_params$mu_alpha
  V_alpha <- final_hyper_params$V_alpha
  nu_0 <- final_hyper_params$nu_0
  a_nu <- final_hyper_params$a_nu
  b_nu <- final_hyper_params$b_nu
  a_rho <- final_hyper_params$a_rho
  b_rho <- final_hyper_params$b_rho
  a_lamb <- final_hyper_params$a_lamb
  b_lamb <- final_hyper_params$b_lamb
  a_tau <- final_hyper_params$a_tau
  b_tau <- final_hyper_params$b_tau
  h_1 <- final_hyper_params$h_1
  h_2 <- final_hyper_params$h_2


  # t_min = min(t,na.rm = TRUE)
  # t_max = max(t,na.rm = TRUE)
  I = dim(y)[1]
  Q = dim(y)[2]
  J_max = dim(y)[3]
  S = dim(Z)[4]
  K = max(g, na.rm = TRUE)

  num_intervals <- df - 3
  quantiles <- quantile(t, probs = seq(0, 1, 1/num_intervals), na.rm = TRUE)
  middle_quantiles <- quantiles[-c(1, num_intervals + 1)]

  # print(paste("t_min",t_min,"t_max",t_max))
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


  if(Q > 1){
    V_alpha_inv <- chol2inv(chol(V_alpha))
    V_alpha_inv_mu_alpha <- V_alpha_inv %*% mu_alpha

    Z_sum <- array(0, dim=c(Q, S, S))


    for (q in 1:Q){
      for (i in 1:I){
        k <- g[i]
        data_index_iq <- which(data_index[i,q,] == 1)

        if(length(data_index_iq)>1){
          Z_sum[q,,] <- Z_sum[q,,] + t(Z[i,q,data_index_iq,])%*%Z[i,q,data_index_iq,]
        } else if(length(data_index_iq)==1) {
          Z_iq <- as.matrix(Z[i, q, data_index_iq, ])
          Z_sum[q,,] <- Z_sum[q,,] + Z_iq %*% t(Z_iq)
        }

      }
    }


    B[is.na(B)] <- 0
    Z[is.na(Z)] <- 0
    y[is.na(y)] <- 0

    B_cpp <- vector("list", I)
    for (i in 1:I) {
      B_cpp[[i]] <- B[i, , , ]
    }

    Z_cpp <- vector("list", I)
    for (i in 1:I) {
      Z_cpp[[i]] <- Z[i, , , ]
    }

    g_cpp = g-1


    # Initialize MCMC storage object
    mcmc <- list()

    # Initialize alpha, beta, omega, Sigma_omega, sigma2
    mcmc$alpha <- array(NA, dim=c(Nit, Q, S))
    # mcmc$omega <- array(NA, dim=c(I, Q))
    mcmc$Sigma_omega <- array(NA, dim=c(Nit, Q, Q))
    mcmc$sigma2 <- array(NA, dim=c(Nit, Q))
    # Store additional parameters
    mcmc$m_kq <- array(NA, dim=c(Nit, K, Q, L-1))
    mcmc$xi_kq <- array(NA, dim=c(Nit, K, Q, L-1))
    mcmc$rho <- array(NA, dim=c(Nit, 1))
    mcmc$rho_0 <- array(NA,dim = c(Nit,1))
    mcmc$nu_kq <- array(NA, dim=c(Nit, K, Q))
    mcmc$nu_kq0 <- array(NA, dim=c(Nit, K, Q))
    mcmc$gamma_kq <- array(NA, dim=c(Nit, K, Q))
    mcmc$gamma_kq0 <- array(NA, dim=c(Nit, K, Q))
    mcmc$eta_kq <- array(NA, dim=c(Nit, K, Q))
    mcmc$beta_kq0 <- array(NA, dim=c(Nit, K, Q))
    mcmc$beta_wo_intcp <- array(NA, dim=c(Nit, K, Q, L-1))
    mcmc$beta_whole <- array(NA,dim = c(Nit, K,Q, L))
    mcmc$tau_vec = array(NA,dim = c(Nit,Q))
    mcmc$lambda_vec = array(NA,dim = c(Nit,Q))

    initial <- init(Q,S,K,L,I,J_max,data_index,t)


    # Store initial values in the first row of mcmc
    mcmc$alpha[1, , ] <- initial$alpha
    # mcmc$omega[ , ] <- initial$omega
    current_omega <- initial$omega
    mcmc$Sigma_omega[1, , ] <- initial$Sigma_omega
    mcmc$sigma2[1, ] <- initial$sigma2
    mcmc$m_kq[1, , , ] <- initial$m_kq
    mcmc$xi_kq[1, , , ] <- initial$xi_kq
    mcmc$rho[1] <- initial$rho
    mcmc$rho_0[1] <- initial$rho_0
    mcmc$nu_kq[1, , ] <- initial$nu_kq
    mcmc$nu_kq0[1, , ] <- initial$nu_kq0
    mcmc$gamma_kq[1, , ] <- initial$gamma_kq
    mcmc$gamma_kq0[1, , ] <- initial$gamma_kq0
    mcmc$eta_kq[1, , ] <- initial$eta_kq
    mcmc$beta_kq0[1, , ] <- initial$beta_kq0
    mcmc$beta_wo_intcp[1, , , ] <- initial$beta_wo_intcp
    mcmc$beta_whole[1, , , ] <- initial$beta_whole
    mcmc$tau_vec[1, ] = initial$tau_vec
    mcmc$lambda_vec[1,] = initial$lambda_vec
    current_theta_iq = initial$theta_iq

    for (nit in 2:Nit) {
      if (display_process == TRUE && nit %% 1000 == 0) {
        print(paste0("current iteration ", nit))
      }
      mcmc$alpha[nit,,] <- update_alpha_cpp(mcmc$beta_whole[nit-1,,,], current_omega,current_theta_iq,
                                            c(mcmc$sigma2[nit-1,]),data_index,
                                            y, B_cpp, Z_cpp, g_cpp, Z_sum,  V_alpha_inv,
                                            V_alpha_inv_mu_alpha)


      mcmc$eta_kq[nit, , ] = update_eta_kq_cpp(mcmc$alpha[nit,,], mcmc$beta_wo_intcp[nit-1, , , ],
                                               current_omega,current_theta_iq,
                                               mcmc$sigma2[nit-1, ], mcmc$beta_kq0[nit-1, , ], mcmc$xi_kq[nit-1, , ,],
                                               mcmc$gamma_kq[nit-1, , ], mcmc$nu_kq[nit-1, , ],
                                               y, Z_cpp, B_cpp,g_cpp,data_index
      )


      mcmc$xi_kq[nit, , , ] = update_xi_kq_cpp(mcmc$alpha[nit, ,], mcmc$beta_wo_intcp[nit-1, , , ],
                                               current_omega,current_theta_iq,
                                               c(mcmc$sigma2[nit-1,]), mcmc$beta_kq0[nit-1, , ],
                                               mcmc$eta_kq[nit, ,], mcmc$m_kq[nit-1, , , ],
                                               g_cpp, data_index,y,Z_cpp,B_cpp)

      update_beta_wo_eta_xi_std_list = update_beta_wo_eta_xi_std(mcmc$eta_kq[nit, , ], mcmc$xi_kq[nit, , , ])
      mcmc$beta_wo_intcp[nit, , , ] = update_beta_wo_eta_xi_std_list$current_beta_wo_intcp
      mcmc$eta_kq[nit, , ] = update_beta_wo_eta_xi_std_list$eta_update_std
      mcmc$xi_kq[nit, , ,] = update_beta_wo_eta_xi_std_list$xi_update_std

      mcmc$gamma_kq[nit, , ] = update_gamma_kq(mcmc$eta_kq[nit, , ],mcmc$rho[nit-1],nu_0, mcmc$nu_kq[nit-1, , ])
      mcmc$nu_kq[nit, ,] = update_nu_kq(a_nu,b_nu,mcmc$eta_kq[nit, ,],mcmc$gamma_kq[nit, ,])
      mcmc$rho[nit] = update_rho_kq(a_rho,b_rho,mcmc$gamma_kq[nit, ,])
      mcmc$m_kq[nit, ,,] = update_mkql(mcmc$xi_kq[nit, ,,])

      # update intercept
      mcmc$beta_kq0[nit, ,] = update_beta_kq0_cpp(mcmc$alpha[nit,,], mcmc$beta_wo_intcp[nit, , ,],
                                                  current_omega,current_theta_iq,
                                                  mcmc$sigma2[nit-1,],
                                                  mcmc$gamma_kq0[nit-1,,],mcmc$nu_kq0[nit-1,,],
                                                  y,Z_cpp,B_cpp,g_cpp,data_index)
      # update the hyper parameters for intercept
      mcmc$gamma_kq0[nit, ,] = update_gamma_kq(mcmc$beta_kq0[nit, ,], mcmc$rho_0[nit-1],nu_0,
                                               mcmc$nu_kq0[nit-1,,],intercept = TRUE)
      mcmc$nu_kq0[nit, ,] = update_nu_kq(a_nu,b_nu,mcmc$beta_kq0[nit, ,],mcmc$gamma_kq0[nit, ,],
                                         intercept = TRUE)
      mcmc$rho_0[nit] = update_rho_kq(a_rho,b_rho,mcmc$gamma_kq0[nit, ,],intercept = TRUE)
      # update a complete beta
      mcmc$beta_whole[nit, , , ] = update_whole_beta(mcmc$beta_kq0[nit, ,],mcmc$beta_wo_intcp[nit, , ,])

      # update omega
      # mcmc$omega[nit,,] = omega
      current_omega <- update_omega_cpp(mcmc$alpha[nit,,], mcmc$beta_whole[nit, , , ],
                                        current_omega,current_theta_iq,
                                        c(mcmc$sigma2[nit-1,]), mcmc$Sigma_omega[nit-1,,],
                                        y, Z_cpp, B_cpp, g_cpp,
                                        data_index)
      # update Sigma_omega
      # mcmc$Sigma_omega[nit,,] = Sigma_omega
      mcmc$Sigma_omega[nit,,] <- update_Sigma_omega_cpp(current_omega, mcmc$Sigma_omega[nit-1,,])


      # update Sigma_omega
      current_theta_iq = update_theta_iq_cpp(mcmc$alpha[nit,,], mcmc$beta_whole[nit, , , ],
                                             current_omega,c(mcmc$sigma2[nit-1,]),
                                             data_index, y, B_cpp, Z_cpp,g_cpp,
                                             t, mcmc$tau_vec[nit-1,],mcmc$lambda_vec[nit-1,])

      # update tau_q

      mcmc$tau_vec[nit, ] = update_tau_q_cpp(current_theta_iq, a_tau, b_tau,
                                             data_index, t, mcmc$lambda_vec[nit-1,])

      # update lambda_q
      mcmc$lambda_vec[nit, ] = update_lambda_q_cpp( mcmc$lambda_vec[nit-1,], a_lamb, b_lamb,
                                                    current_theta_iq, mcmc$tau_vec[nit, ],
                                                    data_index, t, 0.05)


      # update sigma2
      # mcmc$sigma2[nit,] = sigma2
      mcmc$sigma2[nit,] <- update_sigma2_cpp(mcmc$alpha[nit,,], mcmc$beta_whole[nit, , , ],
                                             current_omega, current_theta_iq,
                                             y, Z_cpp, B_cpp, data_index, g_cpp,h_1, h_2)
    }




    # mcmc=mcmc_result
    total_iterations = Nit
    # Calculate the number of posterior samples
    num_samples <- total_iterations - burn_in
    # Generate the index list
    index_list <- seq(burn_in+1, total_iterations, by = thin_factor)


    post_alpha_samples = mcmc$alpha[index_list, , ]
    post_beta_samples = mcmc$beta_whole[index_list, , , ]
    post_Sigma_omega_samples = mcmc$Sigma_omega[index_list, , ]
    post_sigma_samples = mcmc$sigma2[index_list, ]
    post_gamma_kq_samples = mcmc$gamma_kq[index_list, , ]
    post_gamma_kq0_samples = mcmc$gamma_kq0[index_list, , ]
    post_tau_samples = mcmc$tau_vec[index_list,  ]
    post_lambda = mcmc$lambda_vec[index_list , ]



    post_samples <- list(
      pos.alpha = post_alpha_samples,
      pos.beta = post_beta_samples,
      pos.Sigma_omega = post_Sigma_omega_samples,
      pos.tau_q = post_tau_samples,
      pos.lambda_q = post_lambda,
      pos.sigma2 = post_sigma_samples,
      pos.gamma_kq = post_gamma_kq_samples,
      pos.gamma_kq0 = post_gamma_kq0_samples

    )
    class(post_samples) = "Post_BayTetra"
    # Return only the posterior samples
    return(post_samples)

  }else{

    # start to collapse the dim

    y_cs = y[,1,]
    Z_cs = Z[,1, , ]
    B_cs = B[,1, , ]
    data_index_cs = data_index[,1,]



    V_alpha_inv <- chol2inv(chol(V_alpha));
    V_alpha_inv_mu_alpha <- V_alpha_inv %*% mu_alpha;


    # pre-calculated data for saving running time
    Z_sum <- array(0, dim=c(S, S))


    for (i in 1:I){
      k <- g[i]
      data_index_iq <- which(data_index_cs[i,] == 1)

      if(length(data_index_iq)>1){
        Z_sum <- Z_sum + t(Z_cs[i,data_index_iq,])%*%Z_cs[i,data_index_iq,]

      } else if(length(data_index_iq)==1) {
        Z_iq <- as.matrix(Z_cs[i, data_index_iq, ])
        Z_sum <- Z_sum + Z_iq %*% t(Z_iq)

      }
    }

    y = y_cs
    Z = Z_cs
    B = B_cs
    data_index = data_index_cs


    B[is.na(B)] <- 0
    Z[is.na(Z)] <- 0
    y[is.na(y)] <- 0


    g_cpp = g-1 # to match the index in Rcpp, so that g[i] from 0 to K-1


    # Initialize MCMC storage object
    mcmc <- list()

    # Initialize alpha, beta, omega, Sigma_omega, sigma2
    mcmc$alpha <- array(NA, dim=c(Nit,  S))
    mcmc$sigma2 <- array(NA, dim=c(Nit))
    # Store additional parameters
    mcmc$m_kq <- array(NA, dim=c(Nit, K, L-1))
    mcmc$xi_kq <- array(NA, dim=c(Nit, K, L-1))
    mcmc$rho <- array(NA, dim=c(Nit))
    mcmc$rho_0 <- array(NA,dim = c(Nit))
    mcmc$nu_kq <- array(NA, dim=c(Nit, K))
    mcmc$nu_kq0 <- array(NA, dim=c(Nit, K))
    mcmc$gamma_kq <- array(NA, dim=c(Nit, K))
    mcmc$gamma_kq0 <- array(NA, dim=c(Nit, K))
    mcmc$eta_kq <- array(NA, dim=c(Nit, K))
    mcmc$beta_kq0 <- array(NA, dim=c(Nit, K))
    mcmc$beta_wo_intcp <- array(NA, dim=c(Nit, K,  L-1))
    mcmc$beta_whole <- array(NA,dim = c(Nit, K, L))

    initial <- init_Q1(Q,S,K,L,I)


    # Store initial values in the first row of mcmc
    mcmc$alpha[1, ] <- initial$alpha
    mcmc$sigma2[1] <- initial$sigma2
    mcmc$m_kq[1, ,  ] <- initial$m_kq
    mcmc$xi_kq[1, , ] <- initial$xi_kq
    mcmc$rho[1] <- initial$rho
    mcmc$rho_0[1] <- initial$rho
    mcmc$nu_kq[1, ] <- initial$nu_kq
    mcmc$nu_kq0[1, ] <- initial$nu_kq0
    mcmc$gamma_kq[1, ] <- initial$gamma_kq
    mcmc$gamma_kq0[1,  ] <- initial$gamma_kq
    mcmc$eta_kq[1, ] <- initial$eta_kq
    mcmc$beta_kq0[1, ] <- initial$beta_kq0
    mcmc$beta_wo_intcp[1, , ] <- initial$beta_wo_intcp
    mcmc$beta_whole[1, , ] <- initial$beta_whole


    for (nit in 2:Nit) {
      if (display_process == TRUE && nit %% 1000 == 0) {
        print(paste0("current iteration ", nit))
      }
      mcmc$alpha[nit,] <- update_alpha_Q1_cpp(mcmc$beta_whole[nit-1,,],mcmc$sigma2[nit-1],
                                              data_index, y, B, Z, g_cpp, Z_sum,  V_alpha_inv,
                                              V_alpha_inv_mu_alpha)
      # update eta_kq
      mcmc$eta_kq[nit, ] = update_eta_kq_Q1_cpp(mcmc$alpha[nit,], mcmc$beta_wo_intcp[nit-1, , ],
                                                mcmc$sigma2[nit-1],
                                                mcmc$beta_kq0[nit-1, ], mcmc$xi_kq[nit-1, , ],
                                                mcmc$gamma_kq[nit-1, ], mcmc$nu_kq[nit-1,  ],
                                                y, Z, B,g_cpp,data_index)

      mcmc$xi_kq[nit, , ] = update_xi_kq_Q1_cpp(mcmc$alpha[nit, ], mcmc$beta_wo_intcp[nit-1, , ],
                                                mcmc$sigma2[nit-1], mcmc$beta_kq0[nit-1, ],
                                                mcmc$eta_kq[nit, ], mcmc$m_kq[nit-1, ,  ],
                                                g_cpp, data_index,y,Z,B)

      # standardize these eta_kq and xi_kq
      update_beta_wo_eta_xi_std_list = update_beta_wo_eta_xi_std_Q1(mcmc$eta_kq[nit, ], mcmc$xi_kq[nit, , ])
      mcmc$beta_wo_intcp[nit, , ] = update_beta_wo_eta_xi_std_list$current_beta_wo_intcp
      mcmc$eta_kq[nit,  ] = update_beta_wo_eta_xi_std_list$eta_update_std
      mcmc$xi_kq[nit, , ] = update_beta_wo_eta_xi_std_list$xi_update_std

      # update hyper parameter for non-intercept
      mcmc$gamma_kq[nit,  ] = update_gamma_kq_Q1(mcmc$eta_kq[nit,  ],mcmc$rho[nit-1],nu_0, mcmc$nu_kq[nit-1,  ])

      mcmc$nu_kq[nit, ] = update_nu_kq_Q1(a_nu,b_nu,mcmc$eta_kq[nit, ],mcmc$gamma_kq[nit, ])
      mcmc$rho[nit] = update_rho_kq_Q1(a_rho,b_rho,mcmc$gamma_kq[nit, ])

      mcmc$m_kq[nit, ,] = update_mkql_Q1(mcmc$xi_kq[nit, ,])


      # update intercept
      mcmc$beta_kq0[nit, ] = update_beta_kq0_Q1_cpp(mcmc$alpha[nit,], mcmc$beta_wo_intcp[nit, , ],
                                                    mcmc$sigma2[nit-1], mcmc$gamma_kq0[nit-1,],
                                                    mcmc$nu_kq0[nit-1,],
                                                    y,Z,B,g_cpp,data_index)

      # update the hyper parameters for intercept
      mcmc$gamma_kq0[nit, ] = update_gamma_kq_Q1(mcmc$beta_kq0[nit, ], mcmc$rho_0[nit-1],nu_0,
                                                 mcmc$nu_kq0[nit-1,],intercept = TRUE)

      mcmc$nu_kq0[nit, ] = update_nu_kq_Q1(a_nu,b_nu,mcmc$beta_kq0[nit, ],mcmc$gamma_kq0[nit, ],
                                           intercept = TRUE)
      mcmc$rho_0[nit] = update_rho_kq_Q1(a_rho,b_rho,mcmc$gamma_kq0[nit, ],intercept = TRUE)
      # update a complete beta
      mcmc$beta_whole[nit, , ] = update_whole_beta_Q1(mcmc$beta_kq0[nit, ],mcmc$beta_wo_intcp[nit, , ])



      # update sigma2
      mcmc$sigma2[nit] <- update_sigma2_Q1_cpp(mcmc$alpha[nit,], mcmc$beta_whole[nit, , ],
                                               y,Z,B, data_index , g_cpp,h_1,h_2)

    }




    # mcmc=mcmc_result
    total_iterations = Nit
    # Calculate the number of posterior samples
    num_samples <- total_iterations - burn_in
    # Generate the index list
    index_list <- seq(burn_in+1, total_iterations, by = thin_factor)




    post_alpha_samples = mcmc$alpha[index_list,  ]
    post_beta_samples = mcmc$beta_whole[index_list,  , ]
    post_sigma_samples = mcmc$sigma2[index_list]
    post_gamma_kq_samples = mcmc$gamma_kq[index_list, ]
    post_gamma_kq0_samples = mcmc$gamma_kq0[index_list, ]


    post_samples <- list(
      pos.alpha = post_alpha_samples,
      pos.beta = post_beta_samples,
      pos.sigma2 = post_sigma_samples,
      pos.gamma_kq = post_gamma_kq_samples,
      pos.gamma_kq0 = post_gamma_kq0_samples
    )
    class(post_samples) = "Post_BayTetra"
    # Return only the posterior samples
    return(post_samples)



  }





}























