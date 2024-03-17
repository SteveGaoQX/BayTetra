//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>


using namespace std;
using namespace arma;
using namespace Rcpp;


const double log2pi = log(2.0*M_PI);
//' @noRd

 // [[Rcpp::export]]
 double dmvn_rcpp(const arma::rowvec& x,const arma::rowvec& mean,const arma::mat& sigma, bool logd = false){
   // calculate density of multivariate normal distribution
   // args: x: row vector data
   //      mean: row vector mean, sigma: covariance matrix
   //      logd: true for taking log
   // returns: out: pdf (or log pdf) of multivariate normal distribution

   int xdim = x.size();
   arma::mat rooti = trans(inv(trimatu(chol(sigma))));
   double rootisum = sum(log(rooti.diag()));
   double constants = -(static_cast<double>(xdim)/2.0)*log2pi;

   arma::vec z = rooti*trans(x-mean);
   double out = constants-0.5*sum(z%z)+rootisum;

   if (logd == false){ out = exp(out); }
   return(out);
 }



//' @noRd
 // [[Rcpp::export]]
 arma::cube Generate_t_cpp(int I, int Q, int J_max,
                           double min_dist, double lower, double upper,
                           const arma::cube data_index) {
   arma::cube t(I, Q, J_max, arma::fill::none);

   //double range = upper - lower - (min_dist * (J_max - 1));


   for (int i = 0; i < I; ++i) {
     double q = 0;
     arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1);
     int J_iq = data_index_iq.n_elem;

     // cout << "J_iq" << J_iq << endl;
     arma::vec t_values = arma::vec(J_max, fill::none);
     double interval_size = (upper - lower - (min_dist * (J_iq - 1)))/J_iq;
     // cout << "interval_size" << interval_size << endl;

     for (int j = 0; j < J_iq; ++j) {
       // Define the bounds of the interval
       double interval_lower = lower + j * (interval_size + min_dist);
       double interval_upper = interval_lower + interval_size;

       // cout << "j th (in cpp):" << j << " interval :" << interval_lower << "&" << interval_upper << endl;
       // Sample from the interval
       t_values(data_index_iq(j)) = R::runif(interval_lower, interval_upper);
     }
     // Assign the adjusted values to the cube
     //t.tube(i, q) = t_values;
     for (int q = 0; q < Q; ++q) {
       // Assign the adjusted values to the cube
       t.tube(i, q) = t_values;
     }
   }

   return t;
 }




//' @noRd

 // [[Rcpp::export]]
 arma::mat rmvn_rcpp(const int n,const arma::vec& mean,const arma::mat& sigma){
   // randomly generate samples from multivariate normal distribution
   // args: n: number of data
   //      mean: row vector mean, sigma: covariance matrix
   // returns: out: random samples from multivariate normal distribution

   int k = sigma.n_cols; // dimension of the multivariate normal distribution
   arma::mat z = randn(n, k);
   arma::mat out = repmat(mean,1,n).t()+z*chol(sigma);
   return(out);
 }


//' @noRd
// [[Rcpp::export]]
double rinvgamma_rcpp(double a, double b){

  // generate random samples from inverse-gamma distribution
  // args: inverse-gamma(a, b)
  // returns: random sample from inverse-gamma distribution

  return(1/R::rgamma(a, 1/b));
}

//' @noRd
// [[Rcpp::export]]
double dinv_gamma_rcpp(double x, double alpha, double beta, bool logd = false) {
  if (x <= 0) return logd ? R_NegInf : 0.0;

  double log_pdf = alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - beta / x;
  return logd ? log_pdf : exp(log_pdf);
}


//' @noRd
// [[Rcpp::export]]
arma::mat update_alpha_cpp(const arma::cube& beta, const arma::mat& beta_kq0, 
                           const arma::mat& omega, const arma::cube& theta_iq, const arma::vec& sigma2,
                           const arma::cube& y,const arma::cube& data_index,
                           const arma::field<arma::cube>& B, const arma::field<arma::cube>& Z, 
                           const arma::vec& g_cpp, const arma::cube& Z_sum,
                           const arma::mat& V_alpha_inv, const arma::vec& V_alpha_inv_mu_alpha){
  // Update the alpha
  // args: 1: beta is K*Q*L at t-1 time
  //       2: beta_kq0 is K*Q dim 
  //       3: beta_kqt is K*Q
  //       4: omega is the I*Q mat at t-1 time
  //       5: theta_iq dim I*Q*J_max
  //       6: sigma2 is the length Q col vect at t-1 time
  //       7: data_index dim I*Q*J_max
  //       8: y dim I*Q*J_max
  //       9: t dim I*Q*J_max
  //       10: B dim I*Q*J_max*L (L is the length after shrinkage, here we do not add intercept on B)
  //       11: Z dim I*Q*J_max*S
  //       12: g vec length I
  //       13: Z_sum dim Q*S*S (?)
  
  // returns: the updates alpha at dim Q*S
  
  
  // Get dimensions from the data
  int K = beta.n_rows;
  int Q = beta.n_cols;
  int S = Z_sum.n_cols; // Z_sum in dim Q*S*S
  int I = g_cpp.n_elem;
  int J_max = data_index.n_slices; // data_index in dim I*Q*J_max
  //cout << "K:" << K << "  Q: "<< Q <<endl;
  arma::mat alpha_update(Q, S, arma::fill::none);
  arma::mat V_n(S, S);
  arma::vec mu_n(S);
  
  for (int q = 0; q < Q; ++q){
    
    arma::vec Zy_sum(S, arma::fill::zeros);
    
    for (int i = 0; i < I; i++){
      
      int k = g_cpp(i);
      //cout<< "i:  "<< i <<" k: "<< k<<endl;
      arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1);
      
      arma::mat B_iq_raw = B(i).tube(q,0,q, J_max-1 );
      arma::mat B_iq = B_iq_raw.rows(data_index_iq);
      arma::vec beta_kq = beta.tube(k, q);
      arma::vec beta_Kq = beta.tube(K-1, q);
      
      // arma::vec t_vec = t.tube(i, q);
      // arma::vec t_selected = t_vec.elem(data_index_iq);
      // 
      // arma::vec t_beta_kqt = beta_kqt(k,q)*t_selected;
      // arma::vec t_beta_Kqt = u(K-1,q)*t_selected;
      
      arma::vec beta_kq0_vec(data_index_iq.size(), arma::fill::zeros);
      beta_kq0_vec.fill(beta_kq0(k, q));
      
      arma::vec omega_vec(data_index_iq.size(), arma::fill::zeros);
      omega_vec.fill(omega(i, q));
      
      
      // extract theta_iq
      arma::vec cur_theta_iq = arma::vectorise(theta_iq.tube(i,q));
      arma::vec select_cur_theta_iq = cur_theta_iq.elem(data_index_iq);
      
      arma::vec y_vec = arma::vectorise(y.tube(i, q));
      arma::vec y_selected = y_vec.elem(data_index_iq);
      arma::vec y_tilde_iq;
      if (k == K-1){ // baseline
        y_tilde_iq = y_selected - B_iq * beta_Kq - omega_vec - select_cur_theta_iq;
      } else { // not baseline
        y_tilde_iq = y_selected - B_iq * beta_Kq - beta_kq0_vec - B_iq * beta_kq -
          omega_vec - select_cur_theta_iq;
      }
      arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1);
      arma::mat Z_selected = Z_iq_raw.rows(data_index_iq);
      Zy_sum += Z_selected.t() * y_tilde_iq;
    }
    arma::mat Z_sum_q = Z_sum.tube(q,0, q,S-1);
    V_n = arma::inv(Z_sum_q / sigma2(q) + V_alpha_inv);
    mu_n = V_n * (Zy_sum / sigma2(q) + V_alpha_inv_mu_alpha);
    alpha_update.row(q) = rmvn_rcpp(1, mu_n, V_n).row(0);
  }
  return alpha_update;
}


//' @noRd
 // [[Rcpp::export]]
 arma::mat update_beta_kq0_cpp(const arma::mat& alpha, const arma::cube& beta, 
                               const arma::mat& omega, const arma::cube& theta_iq, const arma::vec& sigma2,
                               const arma::mat& nu_kq0,
                               const arma::cube& y, const arma::cube& data_index, 
                               const arma::field<arma::cube>& B, const arma::field<arma::cube>& Z,const arma::vec& g_cpp){
   
   int K = beta.n_rows;
   int Q = beta.n_cols;
   //int L = B(0).n_slices;
   //int n = y.n_rows;
   int J_max = B(0).n_cols;
   arma::mat beta_kq0_update(K, Q, arma::fill::zeros);
   
   for (int k = 0; k < K-1; ++k){
     for (int q = 0; q < Q; ++q){
       double Sigma_beta_kq0_sum = 0;
       double mu_beta_kq0_sum = 0;
       arma::uvec index_k = find(g_cpp == k);
       
       for (auto i : index_k){
         arma::uvec data_index_iq = find(data_index.tube(i,q) == 1);
         //Rcpp::Rcout << "data_index_iq: " << data_index_iq << std::endl;
         int J_iq = data_index_iq.size();
         
         arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
         arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
         arma::vec alpha_q = alpha.row(q).t();
         
         // take B[i,q,index_index_iq,1:L], i.e do not count intercept part (we do not have intercept 1 on B)
         arma::mat B_iq_raw = B(i).tube(q,0,q, J_max-1);
         arma::mat B_iq = B_iq_raw.rows(data_index_iq);
         
         arma::vec beta_Kq = beta.tube(K-1, q); 
         arma::vec beta_kq = beta.tube(k, q);
         
         // take vec_t*beta_kqt or vec_t*beta_Kqt
         // arma::vec t_vec = t.tube(i, q);
         // arma::vec t_selected = t_vec.elem(data_index_iq);
         // 
         // arma::vec t_beta_kqt = beta_kqt(k,q)*t_selected;
         // arma::vec t_beta_Kqt = beta_kqt(K-1,q)*t_selected;
         
         
         // take omega_i
         arma::vec omega_vec(J_iq, arma::fill::zeros);
         omega_vec.fill(omega(i, q));
         
         // extract theta_iq
         arma::vec cur_theta_iq = arma::vectorise(theta_iq.tube(i,q));
         arma::vec select_cur_theta_iq = cur_theta_iq.elem(data_index_iq);
         
         arma::vec y_vec = arma::vectorise(y.tube(i, q));
         arma::vec y_selected = y_vec.elem(data_index_iq);
         
         arma::vec y_tilde_iq;
         
         y_tilde_iq = y_selected - Z_iq * alpha_q - 
           B_iq* beta_Kq - B_iq*beta_kq - omega_vec - select_cur_theta_iq;
         
         mu_beta_kq0_sum = mu_beta_kq0_sum + sum(y_tilde_iq);
         
         Sigma_beta_kq0_sum = Sigma_beta_kq0_sum + J_iq;
         
       }
       
       double sigma_beta_kq0 = 1/ (Sigma_beta_kq0_sum/sigma2(q) + 1/(nu_kq0(k,q)));
       double mu_beta_kq0 = sigma_beta_kq0*mu_beta_kq0_sum/sigma2(q);
       beta_kq0_update(k,q) = Rcpp::rnorm(1, mu_beta_kq0, sqrt(sigma_beta_kq0))[0];
       
     }
   }
   return beta_kq0_update;
 }




//' @noRd
//[[Rcpp::export]]
arma::cube update_beta_kq_cpp(const arma::mat& alpha, const arma::cube& beta, const arma::mat& beta_kq0,
                              const arma::mat& omega, const arma::cube& theta_iq,const arma::vec& sigma2, 
                              const arma::cube& y, const arma::cube& data_index,
                              const arma::mat& K_mat, const arma::mat& tau_kq2,
                              const arma::field<arma::cube>& B, const arma::field<arma::cube>& Z,const arma::vec& g_cpp){
  //cout<< "begin update xi in cpp"<< endl;
  int K = beta.n_rows;
  int Q = beta.n_cols;
  int L = beta.n_slices;
  int I = y.n_rows;
  int J_max = y.n_slices;
  
  arma::cube beta_kq_update(K, Q, L, arma::fill::zeros);
  
  
  // update xi_kq, k = 1,2,...,K-1
  for (int k = 0; k < K-1; k++){
    for (int q = 0; q < Q; q++){
      
      // Indexes for subject in k-th class
      arma::uvec index_k = arma::find(g_cpp == k);
      // Initialize sum values
      arma::mat Sigma_beta_sum(L, L, arma::fill::zeros);
      arma::vec mu_beta_sum(L, arma::fill::zeros);
      for (auto i : index_k){
        
        arma::uvec data_index_iq = arma::find(data_index.tube(i,q) == 1);
        int J_iq = data_index_iq.size();
        
        arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
        arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
        arma::vec alpha_q = alpha.row(q).t();
        
        // take B[i,q,index_index_iq,1:L], i.e do not count intercept part (we do not have intercept 1 on B)
        arma::mat B_iq_raw = B(i).tube(q,0,q, J_max-1);
        arma::mat B_iq = B_iq_raw.rows(data_index_iq);
        
        arma::vec beta_Kq = beta.tube(K-1, q); // this can be put outside the for loop for i
        
        
        // take beta_kq0
        arma::vec beta_kq0_vec(J_iq, arma::fill::zeros);
        beta_kq0_vec.fill(beta_kq0(k,q));
        
        // take omega_i
        arma::vec omega_vec(J_iq, arma::fill::zeros);
        omega_vec.fill(omega(i, q));
        
        // extract theta_iq
        arma::vec cur_theta_iq = arma::vectorise(theta_iq.tube(i,q));
        arma::vec select_cur_theta_iq = cur_theta_iq.elem(data_index_iq);
        
        arma::vec y_vec = arma::vectorise(y.tube(i, q));
        arma::vec y_selected = y_vec.elem(data_index_iq);
        
        arma::vec y_tilde_iq;
        
        y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq*beta_Kq -
          beta_kq0_vec - omega_vec - select_cur_theta_iq;
        
        arma::mat B_iqj = B_iq; //J*L
        arma::vec B_iqj_wcoeff = trans(B_iqj) * y_tilde_iq;
        
        arma::mat vi = trans(B_iqj) * B_iqj; // L*L
        Sigma_beta_sum = Sigma_beta_sum + vi;
        
        mu_beta_sum = mu_beta_sum + B_iqj_wcoeff;
      }
      //cout<< "current k:" <<k << "q:" <<q <<endl;
      //cout<< "Sigma_xi_sum" << Sigma_xi_sum << endl;
      //cout<< "eta(k,q)" << eta(k,q) <<endl;
      
      // posterior mean and covariance matrix
      arma::mat Identity = arma::eye(L, L);
      arma::mat mat = Sigma_beta_sum/sigma2(q) + K_mat/(tau_kq2(k,q));
      
      // check if the matrix is symmetric positive definite
      if(mat.is_sympd()) {
        arma::mat V_xi = arma::inv_sympd(mat);
        // continue your computation here
      } else {
        arma::mat V_xi = Identity;
        Rcpp::Rcout << "The matrix is not symmetric positive definite: " << std::endl;
        mat.print();
        throw std::range_error("Error: Matrix is not symmetric positive definite.");
      }
      //cout << V_xi << endl;
      arma::mat V_beta = arma::inv_sympd(mat);
      arma::vec mu_beta = V_beta*(mu_beta_sum/sigma2(q));
      //cout<< "m_kq" <<m_kq <<endl;
      
      //Rcpp::Rcout << "this is k: " << k + 1 << " q: " << q << "\n";
      //Rcpp::Rcout << "the mean is: " << mu_xi << "\n";
      //Rcpp::Rcout << "the variance is: " << V_xi << "\n";
      //arma::vec xi_temp_K = arma::mvnrnd(mu_xi, V_xi); // Draw from multivariate normal distribution
      //arma::rowvec xi_temp_K_row = xi_temp_K.t(); // Transpose to a row vector
      
      beta_kq_update.tube(k,q) = rmvn_rcpp(1, mu_beta, V_beta).row(0); // Store row vector
    }
  }
  
  try {
    // update xi_Kq
    for (int q = 0; q < Q; q++){
      // Initialize sum values
      arma::mat Sigma_beta_sum_K(L, L, arma::fill::zeros);
      arma::vec mu_beta_sum_K(L, arma::fill::zeros);
      
      for (int i = 0; i < I; i++){
        int k = g_cpp(i);
        
        arma::uvec data_index_iq = arma::find(data_index.tube(i,q) == 1);
        int J_iq = data_index_iq.size();
        
        arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
        arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
        arma::vec alpha_q = alpha.row(q).t();
        
        // take B[i,q,index_index_iq,1:L], i.e do not count intercept part (we do not have intercept 1 on B)
        arma::mat B_iq_raw = B(i).tube(q,0,q, J_max-1);
        arma::mat B_iq = B_iq_raw.rows(data_index_iq);
        
        arma::vec beta_Kq = beta.tube(K-1, q); // this can be put outside the for loop for i
        arma::vec beta_kq = beta.tube(k,q);
        
        // take vec_t*beta_kqt or vec_t*beta_Kqt
        // arma::vec t_vec = t.tube(i, q);
        // arma::vec t_selected = t_vec.elem(data_index_iq);
        // 
        // arma::vec t_beta_kqt = beta_kqt(k,q)*t_selected;
        // arma::vec t_beta_Kqt = beta_kqt(K-1,q)*t_selected;
        
        // take beta_kq0
        arma::vec beta_kq0_vec(J_iq, arma::fill::zeros);
        beta_kq0_vec.fill(beta_kq0(k,q));
        
        // take omega_i
        arma::vec omega_vec(J_iq, arma::fill::zeros);
        omega_vec.fill(omega(i, q));
        
        // extract theta_iq
        arma::vec cur_theta_iq = arma::vectorise(theta_iq.tube(i,q));
        arma::vec select_cur_theta_iq = cur_theta_iq.elem(data_index_iq);
        
        arma::vec y_vec = arma::vectorise(y.tube(i, q));
        arma::vec y_selected = y_vec.elem(data_index_iq);
        
        arma::vec y_tilde_iq;
        
        if (k == K-1){
          y_tilde_iq = y_selected - Z_iq*alpha_q - omega_vec - select_cur_theta_iq;
        } else{
          y_tilde_iq = y_selected - Z_iq*alpha_q - beta_kq0_vec -
            B_iq*beta_kq - omega_vec - select_cur_theta_iq;
        }
        
        arma::mat B_iqj = B_iq;
        arma::vec B_iqj_wcoeff = trans(B_iqj) * y_tilde_iq;
        
        arma::mat vi = trans(B_iqj) * B_iqj;
        Sigma_beta_sum_K = Sigma_beta_sum_K + vi;
        
        mu_beta_sum_K = mu_beta_sum_K + B_iqj_wcoeff;
      }
      
      
      // arma::mat Identity = arma::eye(L, L);
      arma::mat V_beta = arma::inv_sympd(Sigma_beta_sum_K/sigma2(q) + K_mat/(tau_kq2(K-1,q)));
      arma::vec mu_beta = V_beta * (mu_beta_sum_K/sigma2(q));
      
      beta_kq_update.tube(K-1,q) = rmvn_rcpp(1, mu_beta, V_beta).row(0);
    }
  } catch (const std::runtime_error& e) {
    // Catch a std::runtime_error, which is thrown by Armadillo for errors such as dimension mismatches
    Rcpp::Rcout << "Caught a runtime_error exception: " << e.what() << std::endl;
  } catch (const std::exception& e) {
    // Catch any other standard exceptions
    Rcpp::Rcout << "Caught an exception: " << e.what() << std::endl;
  } catch (...) {
    // Catch all other exceptions
    Rcpp::Rcout << "Caught an unknown exception!" << std::endl;
  }
  
  
  //cout << "finish update xi in cpp " << endl;
  return beta_kq_update;
  //gc();
}





//' @noRd
// [[Rcpp::export]]
double logpost_omega_cpp(int i, const arma::rowvec& omega_i, 
                          const arma::mat& alpha, const arma::cube& beta, const arma::mat& beta_kq0,
                          const arma::cube& theta_iq, const arma::vec& sigma2, const arma::mat& Sigma_omega,
                          const arma::cube& y,const arma::cube& data_index, 
                          const arma::field<arma::cube>& B, const arma::field<arma::cube>& Z,  const arma::vec& g_cpp,
                          int K, int Q, int J_max){
   // calculate the log-posterior for omega_i
   // args: index of i-th subject, i from 0 to I-1 since it call by other Rcpp function
   //    omega_i: omega values for i-subjects, as a rowvector length Q
   double logpost = 0.0;
   
   // calculate log-likelihood
   int k = g_cpp(i) ; // index of class for i-th subject
   //Rcpp::Rcout << "k =  " << k << "\n";
   for (int q = 0; q < Q; ++q){
     //Rcpp::Rcout << "q =  " << q << "\n";
     
     arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1); // index of data for i-th subject at q-th response
     int J_iq = data_index_iq.size();
     
     arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
     arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
     arma::vec alpha_q = alpha.row(q).t();
     
     // take B[i,q,index_index_iq,1:L], i.e do not count intercept part (we do not have intercept 1 on B)
     arma::mat B_iq_raw = B(i).tube(q,0,q, J_max-1);
     arma::mat B_iq = B_iq_raw.rows(data_index_iq);
     
     arma::vec beta_Kq = beta.tube(K-1, q); 
     arma::vec beta_kq = beta.tube(k, q);
     
     // // take vec_t*beta_kqt or vec_t*beta_Kqt
     // arma::vec t_vec = t.tube(i, q);
     // arma::vec t_selected = t_vec.elem(data_index_iq);
     // 
     // arma::vec t_beta_kqt = beta_kqt(k,q)*t_selected;
     // arma::vec t_beta_Kqt = beta_kqt(K-1,q)*t_selected;
     
     // take beta_kq0
     arma::vec beta_kq0_vec(J_iq, arma::fill::zeros);
     beta_kq0_vec.fill(beta_kq0(k,q));
     
     // extract theta_iq
     arma::vec cur_theta_iq = arma::vectorise(theta_iq.tube(i,q));
     arma::vec select_cur_theta_iq = cur_theta_iq.elem(data_index_iq);
     
     arma::vec y_vec = arma::vectorise(y.tube(i, q));
     arma::vec y_selected = y_vec.elem(data_index_iq);
     
     arma::vec y_tilde_iq;
     
     if (k == K-1){ // baseline
       y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq - select_cur_theta_iq;
     } else { // not baseline
       y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq - beta_kq0_vec -
         B_iq * beta_kq - select_cur_theta_iq;
     }
     
     arma::rowvec y_tilde_iq_rowvec = y_tilde_iq.t();
     
     arma::vec mu(J_iq, arma::fill::ones);
     
     mu *= omega_i(q);
     arma::rowvec mu_rowvec = mu.t();
     
     arma::mat sigma_mat = arma::eye<arma::mat>(J_iq, J_iq) * sigma2(q);
     logpost += dmvn_rcpp(y_tilde_iq_rowvec, mu_rowvec, sigma_mat, true);
     
   }
   //arma::rowvec omega_i_rowvec = omega_i.t();
   arma::rowvec zeros = arma::zeros<arma::rowvec>(Q);
   arma::rowvec omega_i_temp = omega_i;
   logpost += dmvn_rcpp(omega_i_temp, zeros, Sigma_omega, true);
   
   return logpost;
 }


//' @noRd
// [[Rcpp::export]]
arma::mat update_omega_cpp(const arma::mat& alpha, const arma::cube& beta, const arma::mat& beta_kq0,
                           const arma::mat& omega,const arma::cube& theta_iq,
                           const arma::vec& sigma2, const arma::mat& Sigma_omega,
                           const arma::cube& y, const arma::cube& data_index, 
                           const arma::field<arma::cube>& B, const arma::field<arma::cube>& Z, 
                           const arma::vec& g_cpp, const double omega_var_update = 0.04){
  
  int K = beta.n_rows;
  int Q = beta.n_cols;
  int J_max = Z(0).n_cols;
  arma::mat omega_update = omega;
  int I = omega.n_rows;
  
  for (int i = 0; i < I; i++){
    
    // propose new omega_i
    arma::rowvec omega_i_row = omega_update.row(i);
    arma::vec omega_i = arma::conv_to<arma::vec>::from(omega_i_row); // convert rowvec to vec
    arma::mat var_mat = arma::eye(Q,Q) * omega_var_update;
    arma::mat omega_i_new = rmvn_rcpp(1, omega_i, var_mat);
    
    arma::rowvec omega_i_new_row = omega_i_new.row(0);
    
    // acceptance ratio
    double ratio = logpost_omega_cpp(i, omega_i_new_row, alpha, beta, beta_kq0, theta_iq, sigma2,
                                     Sigma_omega, y, data_index, B, Z,  g_cpp, K, Q, J_max) -
                                       logpost_omega_cpp(i, omega_i_row, alpha, beta, beta_kq0, theta_iq, sigma2,
                                                         Sigma_omega, y, data_index, B, Z, g_cpp, K, Q, J_max);
    
    // accept or reject
    if (log(R::runif(0, 1)) < ratio){
      omega_update.row(i) = omega_i_new.row(0);
    }
  }
  
  return(omega_update);
}



//' @noRd
// [[Rcpp::export]]
double logpost_Sigma_omega_cpp(const arma::mat& Sigma_omega, const arma::mat& omega){
   
   double logpost = 0.0;
   arma::rowvec zeros = arma::zeros<arma::rowvec>(omega.n_cols);
   int I = omega.n_rows; // Added a semicolon at the end
   
   for (int i = 0; i < I; i++){
     arma::rowvec omega_i = omega.row(i);
     logpost += dmvn_rcpp(omega_i, zeros, Sigma_omega, true);
   }
   
   return logpost;
 }




//' @importFrom Rcpp sourceCpp

//' @noRd
// [[Rcpp::export]]
arma::mat update_Sigma_omega_cpp(const arma::mat& omega, const arma::mat& Sigma_omega, const double Sigma_omega_step = 0.008){
   
   arma::mat Sigma_omega_update = Sigma_omega;
   double eps = 1e-300;
   
   int Q = Sigma_omega.n_rows;
   
   for (int i = 0; i < Q; i++){
     for (int j = 0; j < Q; j++){
       
       if (i < j){
         double lower = std::max(Sigma_omega_update(i, j) - Sigma_omega_step, -1.0 + eps);
         double upper = std::min(Sigma_omega_update(i, j) + Sigma_omega_step, 1.0 - eps);
         
         arma::mat Sigma_omega_new = Sigma_omega_update;
         Sigma_omega_new(i, j) = R::runif(lower, upper);
         Sigma_omega_new(j, i) = Sigma_omega_new(i, j);
         
         if (Sigma_omega_new.is_sympd()){
           double ratio = logpost_Sigma_omega_cpp(Sigma_omega_new, omega) -
             logpost_Sigma_omega_cpp(Sigma_omega_update, omega);
           
           if (std::log(R::runif(0, 1)) < ratio){
             Sigma_omega_update = Sigma_omega_new;
           }
         }
       }
     }
   }
   
   return Sigma_omega_update;
 }




//' @noRd
// [[Rcpp::export]]
arma::cube update_theta_iq_cpp(const arma::mat& alpha, const arma::cube& beta, const arma::mat& beta_kq0, 
                                const arma::mat& omega, const arma::vec& sigma2,
                                const arma::cube& y, const arma::cube& data_index, const arma::cube& t_std,const arma::cube& t_org, 
                                const arma::field<arma::cube>& B, const arma::field<arma::cube>& Z, const arma::vec& g_cpp,
                                const arma::vec& tau_q_vec,const arma::vec& lambda_q_vec ){
   
   // Get dimensions from the data
   int K = beta.n_rows;
   int Q = beta.n_cols;
   int I = y.n_rows;
   int J_max = data_index.n_slices; // data_index in dim I*Q*J_max
   //cout << "K:" << K << "  Q: "<< Q <<endl;
   arma::cube theta_iq_update(I, Q, J_max, arma::fill::zeros);
   
   
   for (int q = 0; q < Q; ++q){
     for (int i = 0; i < I; i++){
       
       // extract index
       int k = g_cpp(i);
       arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1);
       int J_iq = data_index_iq.size();
       
       arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
       arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
       arma::vec alpha_q = alpha.row(q).t();
       
       // take B[i,q,index_index_iq,1:L], i.e do not count intercept part (we do not have intercept 1 on B)
       arma::mat B_iq_raw = B(i).tube(q,0,q, J_max-1);
       arma::mat B_iq = B_iq_raw.rows(data_index_iq);
       
       arma::vec beta_Kq = beta.tube(K-1, q); 
       arma::vec beta_kq = beta.tube(k, q);
       
       // take vec_t*beta_kqt or vec_t*beta_Kqt
       arma::vec t_vec = t_std.tube(i, q);
       arma::vec t_selected = t_vec.elem(data_index_iq);
       
       // arma::vec t_beta_kqt = beta_kqt(k,q)*t_selected;
       // arma::vec t_beta_Kqt = beta_kqt(K-1,q)*t_selected;
       
       // take beta_kq0
       arma::vec beta_kq0_vec(J_iq, arma::fill::zeros);
       beta_kq0_vec.fill(beta_kq0(k,q));
       
       // take omega_i
       arma::vec omega_vec(J_iq, arma::fill::zeros);
       omega_vec.fill(omega(i, q));
       
       arma::vec y_vec = arma::vectorise(y.tube(i, q));
       arma::vec y_selected = y_vec.elem(data_index_iq);
       
       arma::vec y_tilde_iq;
       
       // create the dynamic covariance matrix
       arma::mat V_iq(J_iq, J_iq, arma::fill::zeros);
       arma::vec mu_iq(J_iq, arma::fill::zeros);
       
       if (k == K-1){ // baseline
         y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq - omega_vec;
       } else { // not baseline
         y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq - beta_kq0_vec -
           B_iq * beta_kq - omega_vec;
       }
       
       //create the current cov matrix
       double lambda_q = lambda_q_vec(q);
       double tau_q = tau_q_vec(q);
       
       arma::mat cov_mat(J_iq, J_iq, arma::fill::eye);
       for (int j = 0; j < J_iq; j++) {
         for (int j_prime = j + 1; j_prime < J_iq; j_prime++) {
           double diff = (t_org(i, q, data_index_iq(j)) - t_org(i, q, data_index_iq(j_prime)) )/lambda_q;
           double cov_value = std::exp(-std::pow(diff, 2));
           cov_mat(j, j_prime) = cov_value;
           cov_mat(j_prime, j) = cov_value;
         }
       }
       cov_mat = tau_q* cov_mat;
       
       if(abs(arma::det(cov_mat))<1e-35){
         cov_mat.diag() += 1e-4;
       }
       
       arma::mat V_theta_iq_idx = cov_mat;
       
       
       // get Sigma_q in the formula
       arma::mat Sigma_q = sigma2(q) * eye<mat>(J_iq, J_iq);
       // get V_theta in the formula
       V_iq = arma::inv( arma::inv(Sigma_q) + arma::inv(V_theta_iq_idx) );
       // get mu_theta in the formula
       mu_iq = V_iq*(arma::inv(Sigma_q))*y_tilde_iq;
       // do the update
       arma::vec current_theta_iq = rmvn_rcpp(1, mu_iq, V_iq).row(0).t();
       
       // store current update theta_iq to the cube we return
       arma::vec temp_vec(J_max,fill::zeros);
       temp_vec.elem(data_index_iq) = current_theta_iq;
       theta_iq_update.tube(i, q) = temp_vec;
       
     }
   }
   return theta_iq_update;
 }


//' @noRd
// [[Rcpp::export]]
arma::vec update_tau_q_cpp(const arma::cube& theta_iq, double a_tau, double b_tau,
                           const arma::cube& data_index,
                           const arma::cube& t, const arma::vec& lambda_q_vec ){


  // Get dimensions from the data
  int I = theta_iq.n_rows;
  int Q = theta_iq.n_cols;
  //int J_max = theta_iq.n_slices; // data_index in dim I*Q*J_max

  //cout << "K:" << K << "  Q: "<< Q <<endl;
  arma::vec tau_q_update(Q, arma::fill::zeros);


  for (int q = 0; q < Q; ++q){
    double counter = 0;

    double sum_b_tau = 0;
    for (int i = 0; i < I; i++){
      //if (!(i == 0 && q == 0)) continue;  // add this line

      //cout << "q:" << q+1 << "  i: "<< i+1 <<endl;

      // extract index
      arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1);
      int J_iq = data_index_iq.size();

      counter += J_iq;

      // extract theta_iq
      arma::vec cur_theta_iq = arma::vectorise(theta_iq.tube(i,q));

      arma::vec select_cur_theta_iq = cur_theta_iq.elem(data_index_iq);

      // create tilde_V_theta_iq
      double lambda_q = lambda_q_vec(q);
      arma::mat cov_mat(J_iq, J_iq, arma::fill::eye);
      for (int j = 0; j < J_iq; j++) {
        for (int j_prime = j + 1; j_prime < J_iq; j_prime++) {
          double diff = (t(i, q, data_index_iq(j)) - t(i, q, data_index_iq(j_prime)) )/lambda_q;
          double cov_value = std::exp(-std::pow(diff, 2));
          cov_mat(j, j_prime) = cov_value;
          cov_mat(j_prime, j) = cov_value;
        }
      }
      if(arma::det(cov_mat)<1e-35){
        cov_mat.diag() += 1e-3;
      }
      //cout << "cov_mat:" << cov_mat <<endl;
      //cout << "theta_iq" << select_cur_theta_iq <<endl;

      // calculate temp value
      arma::vec temp = arma::inv(cov_mat)*select_cur_theta_iq;
      sum_b_tau +=  dot(select_cur_theta_iq, temp);
      //cout << "sum_tau" << sum_b_tau;
    }
    double update_a_tau = a_tau + 0.5*counter;
    double update_b_tau = b_tau + 0.5*sum_b_tau;

    tau_q_update(q) = rinvgamma_rcpp(update_a_tau, update_b_tau);
  }
  return tau_q_update;
}




//' @noRd
// [[Rcpp::export]]
double logpost_lambda_normpart_cpp(double q, const arma::cube& theta_iq,
                                    double lambda_q, double tau_q ,
                                    const arma::cube& data_index,
                                    const arma::cube& t){
   double threshold = 1e-35;
   double epsilon = 1e-4;
   double logpost_sum = 0.0;
   int I = theta_iq.n_rows;
   
   
   for (int i = 0; i < I; i++){
     
     //cout << "i:" << i << " q: "<< q <<endl; 
     //cout << " current  lambda_q:" << lambda_q << endl;
     
     // value for evaluation
     arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1);
     int J_iq = data_index_iq.size();
     
     // extract theta_iq
     arma::vec cur_theta_iq = arma::vectorise(theta_iq.tube(i,q));
     arma::vec select_cur_theta_iq = cur_theta_iq.elem(data_index_iq);
     
     //create the current cov matrix
     arma::mat cov_mat(J_iq, J_iq, arma::fill::eye);
     for (int j = 0; j < J_iq; j++) {
       for (int j_prime = j + 1; j_prime < J_iq; j_prime++) {
         double diff = (t(i, q, data_index_iq(j)) - t(i, q, data_index_iq(j_prime)) )/lambda_q;
         double cov_value = std::exp(-std::pow(diff, 2));
         cov_mat(j, j_prime) = cov_value;
         cov_mat(j_prime, j) = cov_value;
       }
     }
     cov_mat = tau_q* cov_mat;
     
     double det_val = arma::det(cov_mat);
     
     // If determinant is less than the threshold, add epsilon to the diagonal
     if(abs(det_val) < threshold) {
       cov_mat.diag() += epsilon;
       // return -1e30;
     }
     //cov_mat.diag() += epsilon;
     
     //cout << " cov_mat :" << cov_mat << endl;
     
     // create the mean
     arma::vec mean_vec(J_iq,arma::fill::zeros);
     
     
     // do the transfer
     arma::rowvec mean_rowvec = mean_vec.t();
     arma::rowvec theta_rowvec = select_cur_theta_iq.t();
     
     
     logpost_sum += dmvn_rcpp(theta_rowvec, mean_rowvec, cov_mat, true);
   }
   
   return logpost_sum;
 }


//' @importFrom Rcpp sourceCpp
//' @noRd
// [[Rcpp::export]]
arma::vec update_lambda_q_cpp(const arma::vec& lambda, double a_lam, double b_lam,
                               const arma::cube& theta_iq, const arma::vec& tau_q_vec,
                               const arma::cube& data_index, const arma::cube& t,
                               const double lambda_step = 0.1){
   
   //int I = theta_iq.n_rows;
   int Q = theta_iq.n_cols;
   
   arma::vec lambda_update = lambda;
   double eps = 1e-10;
   int accepted = 0; // Counter for accepted proposals
   
   for (int q = 0; q < Q; ++q){
     
     double lower = std::max(lambda(q) - lambda_step, 0.01 + eps);
     double upper = std::min(lambda(q) + lambda_step, 100 - eps);
     
     double lambda_q_current = lambda(q);
     double lambda_q_new = R::runif(lower, upper);
     
     double tau_q_current = tau_q_vec(q);
     
     double sum_of_norm_new = logpost_lambda_normpart_cpp(q, theta_iq,
                                                          lambda_q_new,  tau_q_current,
                                                          data_index,t );
     double sum_of_norm_cur = logpost_lambda_normpart_cpp(q, theta_iq,
                                                          lambda_q_current,  tau_q_current,
                                                          data_index,t );
     
     double log_IG_new = dinv_gamma_rcpp(lambda_q_new, a_lam, b_lam, true);
     double log_IG_cur = dinv_gamma_rcpp(lambda_q_current, a_lam, b_lam, true);
     
     
     double Ratio = (log_IG_new + sum_of_norm_new) - (log_IG_cur + sum_of_norm_cur);
     
     if(std::log(R::runif(0, 1)) < Ratio){
       lambda_update(q) = lambda_q_new;
       accepted++;
     }
   }
   
   
   return lambda_update;
 }



//' @noRd
// [[Rcpp::export]]
Rcpp::List update_lambda_q_cpp_wth_accept(const arma::vec& lambda, double a_lam, double b_lam,
                              const arma::cube& theta_iq, const arma::vec& tau_q_vec,
                              const arma::cube& data_index, const arma::cube& t,
                              const double lambda_step = 0.1){

  //int I = theta_iq.n_rows;
  int Q = theta_iq.n_cols;

  arma::vec lambda_update = lambda;
  double eps = 1e-10;
  int accepted = 0; // Counter for accepted proposals

  for (int q = 0; q < Q; ++q){

    double lower = std::max(lambda(q) - lambda_step, 0.01 + eps);
    double upper = std::min(lambda(q) + lambda_step, 100 - eps);

    double lambda_q_current = lambda(q);
    double lambda_q_new = R::runif(lower, upper);

    double tau_q_current = tau_q_vec(q);

    double sum_of_norm_new = logpost_lambda_normpart_cpp(q, theta_iq,
                                                         lambda_q_new,  tau_q_current,
                                                         data_index,t );
    double sum_of_norm_cur = logpost_lambda_normpart_cpp(q, theta_iq,
                                                         lambda_q_current,  tau_q_current,
                                                         data_index,t );

    double log_IG_new = dinv_gamma_rcpp(lambda_q_new, a_lam, b_lam, true);
    double log_IG_cur = dinv_gamma_rcpp(lambda_q_current, a_lam, b_lam, true);


    double Ratio = (log_IG_new + sum_of_norm_new) - (log_IG_cur + sum_of_norm_cur);

    if(std::log(R::runif(0, 1)) < Ratio){
      lambda_update(q) = lambda_q_new;
      accepted++;
    }
  }



  // Return both the updated lambda values and the acceptance count
  return Rcpp::List::create(
    Rcpp::Named("lambda") = lambda_update,
    Rcpp::Named("acceptance_count") = accepted
  );

  //return lambda_update;
}



//' @noRd
// [[Rcpp::export]]
arma::vec update_sigma2_cpp(arma::mat& alpha, const arma::cube& beta, const arma::mat& beta_kq0, 
                             arma::mat& omega, const arma::cube& theta_iq,
                             const arma::cube& y, const arma::cube& data_index, 
                             const arma::field<arma::cube>& B, const arma::field<arma::cube>& Z, 
                             const arma::vec& g_cpp,
                             double h_1, double h_2){
   int K = beta.n_rows;
   int Q = beta.n_cols;
   int I = y.n_rows;
   int J_max = Z(0).n_cols;
   
   arma::vec sigma2_update(Q, arma::fill::zeros);
   
   for (int q = 0; q < Q; ++q){
     
     double y_tilde_sum = 0.0;
     
     arma::mat data_index_q = data_index.tube(0, q, I - 1, q);
     int n_iq = arma::accu(data_index_q);
     
     for (int i = 0; i < I; ++i){
       int k = g_cpp(i); // index of class for i-th subject
       
       arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1); // index of data for i-th subject at q-th response
       int J_iq = data_index_iq.size();
       
       arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
       arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
       arma::vec alpha_q = alpha.row(q).t();
       
       // take B[i,q,index_index_iq,1:L], i.e do not count intercept part (we do not have intercept 1 on B)
       arma::mat B_iq_raw = B(i).tube(q,0,q, J_max-1);
       arma::mat B_iq = B_iq_raw.rows(data_index_iq);
       
       arma::vec beta_Kq = beta.tube(K-1, q); 
       arma::vec beta_kq = beta.tube(k, q);
       
       // take vec_t*beta_kqt or vec_t*beta_Kqt
       // arma::vec t_vec = t.tube(i, q);
       // arma::vec t_selected = t_vec.elem(data_index_iq);
       // 
       // arma::vec t_beta_kqt = beta_kqt(k,q)*t_selected;
       // arma::vec t_beta_Kqt = beta_kqt(K-1,q)*t_selected;
       
       // take beta_kq0
       arma::vec beta_kq0_vec(J_iq, arma::fill::zeros);
       beta_kq0_vec.fill(beta_kq0(k,q));
       
       // take omega_i
       arma::vec omega_vec(J_iq, arma::fill::zeros);
       omega_vec.fill(omega(i, q));
       
       // extract theta_iq
       arma::vec cur_theta_iq = arma::vectorise(theta_iq.tube(i,q));
       arma::vec select_cur_theta_iq = cur_theta_iq.elem(data_index_iq);
       
       arma::vec y_vec = arma::vectorise(y.tube(i, q));
       arma::vec y_selected = y_vec.elem(data_index_iq);
       
       arma::vec y_tilde_iq;
       
       if (k == K-1){ // baseline
         y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq - omega_vec - select_cur_theta_iq;
       } else { // not baseline
         y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq - beta_kq0_vec -
           B_iq * beta_kq - omega_vec - select_cur_theta_iq;
       }
       
       // y_tilde_sum += arma::accu(arma::square(y_tilde_iq));
       y_tilde_sum += dot(y_tilde_iq,y_tilde_iq);
     }
     
     double h_1_star = h_1 + n_iq / 2;
     double h_2_star = h_2 + y_tilde_sum / 2;
     
     sigma2_update(q) = rinvgamma_rcpp(h_1_star, h_2_star);
   }
   
   return sigma2_update;
 }





//' @noRd
 // [[Rcpp::export]]
 arma::vec update_alpha_Q1_cpp(const arma::mat& beta, const double& sigma2,
                               const arma::mat& data_index, const arma::mat& y, const arma::cube& B,
                               const arma::cube& Z, const arma::vec& g, const arma::mat& Z_sum,
                               const arma::mat& V_alpha_inv, const arma::vec& V_alpha_inv_mu_alpha){
   // Update the alpha
   // args: 1: beta is K*L at t-1 time
   //       3: sigma2 is double
   // returns: the updates alpha


   // Get dimensions from the data
   int K = beta.n_rows;
   //int Q = beta.n_cols;
   int S = Z_sum.n_cols; // Z_sum in dim S*S
   int I = g.n_elem;
   //cout<< "K:" << K << "S:" << S << "I:" << I <<endl;
   arma::vec alpha_update(S);
   arma::mat V_n(S, S);
   arma::vec mu_n(S);
   arma::vec Zy_sum(S, arma::fill::zeros);

   for (int i = 0; i < I; i++){
     //cout << "i:" << i <<"\n"<< endl;
     int k = g[i] ;
     //cout <<"k: " << k <<"\n" << endl;
     arma::uvec data_index_iq = arma::find(data_index.row(i) == 1);
     arma::mat B_iq_raw  = B.row_as_mat(i).t();
     //cout<< "B_iq_raw:"<< B_iq_raw<<"\n" << endl;

     arma::mat B_iq = B_iq_raw.rows(data_index_iq);
     //cout<< "B_iq" << B_iq<< "\n" << endl;
     arma::vec beta_kq = arma::vectorise(beta.row(k));
     arma::vec beta_Kq = arma::vectorise(beta.row(K-1));
     //cout<< "beta_kq:"<< beta_kq<<"\n" << endl;
     //cout<< "beta_Kq:"<< beta_Kq<<"\n" << endl;

     arma::vec y_vec = arma::vectorise(y.row(i));
     arma::vec y_selected = y_vec.elem(data_index_iq);
     arma::vec y_tilde_iq;
     if (k == K-1){ // baseline
       y_tilde_iq = y_selected - B_iq * beta_Kq;
     } else { // not baseline
       y_tilde_iq = y_selected - B_iq * beta_Kq -
         B_iq * beta_kq;
     }
     //cout<< "y_tilde_iq" << y_tilde_iq<<endl;
     arma::mat Z_iq_raw = Z.row_as_mat(i).t();
     arma::mat Z_selected = Z_iq_raw.rows(data_index_iq);
     //cout << "Z_selected" << Z_selected << endl;
     Zy_sum += Z_selected.t() * y_tilde_iq;
   }
   //cout <<"Zy_sum" << Zy_sum<<"\n"<< endl;
   V_n = arma::inv(Z_sum / sigma2 + V_alpha_inv);
   mu_n = V_n * (Zy_sum / sigma2 + V_alpha_inv_mu_alpha);
   //cout << "V_n" << V_n << "\n"<< endl;
   //cout << "mu_n" <<mu_n<< "\n" << endl;
   alpha_update = mvnrnd(mu_n,V_n);
   //cout<< "K:" << K << "S:" << S << "I:" << I <<endl;
   return alpha_update;
 }


//' @noRd
 // [[Rcpp::export]]
 arma::vec update_eta_kq_Q1_cpp(const arma::vec& alpha, const arma::mat& beta,
                                const double sigma2, const arma::vec& beta_kq0,
                                const arma::mat& xi, const arma::vec& gamma_kq, const arma::vec& nu_kq,
                                const arma::mat& y, const arma::cube& Z, const arma::cube& B, const arma::vec& g_cpp, const arma::mat& data_index){

   int K = beta.n_rows;
   //cout<< "K:"<< K <<" \n ; "<<endl;
   // here the beta is at dim K*L-1
   int L = beta.n_cols + 1;
   //cout<< "L:" << L <<" \n ; "<<endl;
   int I = y.n_rows;
   int J_max = y.n_cols;
   arma::vec eta_update(K, arma::fill::zeros);
   //cout<< "eta_update:" << eta_update<< endl;
   for (int k = 0; k < K-1; ++k){
     double Sigma_eta_sum = 0;
     double mu_eta_sum = 0;
     arma::uvec index_k = find(g_cpp == k);
     arma::vec beta_Kq = arma::vectorise(beta.row(K-1)); // this can be put outside the for loop for i
     //cout<< "beta_wo_intc K:"<< beta_Kq<<endl;
     for (auto i : index_k){
       //cout<< "current i:" << index_k << endl;

       arma::uvec data_index_iq = find(data_index.row(i) == 1);
       int Iq_j = data_index_iq.size();
       arma::vec y_tilde_iq;

       arma::mat B_iq_raw  = B.subcube(i,0,1,i, J_max-1,L-1);
       arma::mat B_iq = B_iq_raw.rows(data_index_iq);
       //cout<< "B_iq"<< B_iq <<"\n"<<endl;
       arma::mat Z_iq_raw = Z.row_as_mat(i).t();
       arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
       //Rcpp::Rcout << "Z_iq: " << Z_iq<<"\n" << std::endl;

       arma::vec y_vec = arma::vectorise(y.row(i));
       arma::vec y_selected = y_vec.elem(data_index_iq);
       //Rcpp::Rcout << "y_selected: " << y_selected <<" \n "<< std::endl;
       arma::vec alpha_q = alpha;

       //
       arma::vec beta_intcp_vec(Iq_j, arma::fill::zeros);
       beta_intcp_vec.fill(beta_kq0(k));

       y_tilde_iq = y_selected - Z_iq * alpha_q -
         B_iq* beta_Kq - beta_intcp_vec;
       //cout<< "y_tilde_iq:"<< y_tilde_iq<<"\n"<< endl;
       arma::vec B_iqj_wcoeff = trans(B_iq) * y_tilde_iq;
       //cout << "B_iqj_wcoeff : " << B_iqj_wcoeff<<" \n " << std::endl;
       arma::vec xi_kq = xi.row(k).t();
       //cout << "xi_kq: " << xi_kq <<" \n " << std::endl;

       mu_eta_sum = mu_eta_sum + dot(B_iqj_wcoeff, xi_kq);
       //cout << "dot result: " << dot(B_iqj_wcoeff, xi_kq) << std::endl;

       arma::vec vi = B_iq * xi_kq;
       //Rcpp::Rcout << "vi: " << vi << std::endl;

       Sigma_eta_sum = Sigma_eta_sum + dot(vi,vi);
       //Rcpp::Rcout << "Sigma_eta_sum: " << Sigma_eta_sum << std::endl;
     }

     double sigma_eta = 1/ (Sigma_eta_sum/sigma2 + 1/(gamma_kq(k)*nu_kq(k)));
     //sigma_eta = abs(sigma_eta);
     //Rcpp::Rcout << "sigma_eta " << sigma_eta << std::endl;

     double mu_eta = sigma_eta*mu_eta_sum/sigma2;
     eta_update(k) = arma::randn( distr_param(mu_eta,sqrt(sigma_eta)) );
     //cout<< "finish k:" << k << endl;
   }
   // cout<<"update baseline"<< endl;
   // update the baseline
   arma::mat Sigma_eta_sum_K(L-1, L-1, arma::fill::zeros);
   arma::vec mu_eta_sum_K(L-1, arma::fill::zeros);

   for (int i = 0; i < I; ++i){
     //cout<<"i: "<< i<< "\n"<<endl;
     int k = g_cpp(i);

     arma::uvec data_index_iq = find(data_index.row(i) == 1);
     int Iq_j = data_index_iq.size();
     arma::vec y_tilde_iq;

     arma::mat B_iq_raw  = B.subcube(i,0,1,i, J_max-1,L-1);
     arma::mat B_iq = B_iq_raw.rows(data_index_iq);
     //cout<<"B_iq:"<<B_iq <<endl;
     arma::mat Z_iq_raw = Z.row_as_mat(i).t();
     arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
     //Rcpp::Rcout << "Z_iq: " << Z_iq << std::endl;

     arma::vec y_vec = arma::vectorise(y.row(i));
     arma::vec y_selected = y_vec.elem(data_index_iq);
     //Rcpp::Rcout << "y_selected: " << y_selected << std::endl;
     arma::vec alpha_q = alpha;

     if (k == K-1){
       y_tilde_iq = y_selected - Z_iq*alpha_q ;
     } else{
       arma::vec beta_kq = arma::vectorise(beta.row(k));
       arma::vec beta_intcp_vec(Iq_j, arma::fill::zeros);
       beta_intcp_vec.fill(beta_kq0(k));

       y_tilde_iq = y_selected - Z_iq*alpha_q -
         B_iq*beta_kq - beta_intcp_vec;
     }
     //cout<< "y_tilde_iq:" << y_tilde_iq<< endl;
     arma::vec B_iqj_wcoeff = trans(B_iq) * y_tilde_iq;
     //cout<< "B_iqj"<<endl;
     mu_eta_sum_K = mu_eta_sum_K + B_iqj_wcoeff;

     arma::mat vi = trans(B_iq) * B_iq;
     //cout<< "vi :"<< vi<< endl;
     Sigma_eta_sum_K = Sigma_eta_sum_K + vi;
   }
   //cout<< "Sigma"<< Sigma_eta_sum_K << endl;
   arma::vec xi_Kq = xi.row(K-1).t();
   //cout<< "xi_Kq: " << xi_Kq << endl;
   arma::mat Sigma_eta_sum_trans = trans(xi_Kq)*Sigma_eta_sum_K*xi_Kq;
   //Rcpp::Rcout << "Sigma_eta_sum_trans " << Sigma_eta_sum_trans << std::endl;

   double trans_scalar = arma::as_scalar(Sigma_eta_sum_trans);
   //Rcpp::Rcout << "trans_scalar " << trans_scalar << std::endl;

   double sigma_eta = 1/ (trans_scalar/sigma2 + 1/(gamma_kq(K-1)*nu_kq(K-1)));
   sigma_eta = abs(sigma_eta);
   //Rcpp::Rcout << "sigma_eta " << sigma_eta << std::endl;

   double mu_eta = sigma_eta*(dot(mu_eta_sum_K, xi_Kq))/sigma2;
   eta_update(K-1) = arma::randn( distr_param(mu_eta,sqrt(sigma_eta)) );


   return eta_update;
 }


//' @noRd
 //[[Rcpp::export]]
 arma::mat update_xi_kq_Q1_cpp(const arma::vec& alpha, const arma::mat& beta,
                               const double sigma2, const arma::vec& beta_kq0, const arma::vec& eta,
                               const arma::mat& m,
                               const arma::vec& g_cpp, const arma::mat& data_index,
                               const arma::mat& y, const arma::cube& Z, const arma::cube& B){
   // beta K*(L-1)
   int K = beta.n_rows;
   //int Q = beta.n_cols;
   int L = beta.n_cols +1;
   int I = y.n_rows;
   int J_max = y.n_cols;

   arma::mat xi_kq_update(K, L-1, arma::fill::zeros);

   // update xi_kq, k = 1,2,...,K-1
   for (int k = 0; k < K-1; k++){

     arma::vec beta_Kq = beta.row(K-1).t(); // this can be put outside the for loop for i

     // Indexes for subject in k-th class
     arma::uvec index_k = arma::find(g_cpp == k);
     // Initialize sum values
     arma::mat Sigma_xi_sum(L-1, L-1, arma::fill::zeros);
     arma::vec mu_xi_sum(L-1, arma::fill::zeros);
     for (auto i : index_k){
       arma::uvec data_index_iq = find(data_index.row(i) == 1);
       int Iq_j = data_index_iq.size();

       arma::mat B_iq_raw  = B.subcube(i,0,1,i, J_max-1,L-1);
       arma::mat B_iq = B_iq_raw.rows(data_index_iq);
       //cout<< "B_iq"<< B_iq <<"\n"<<endl;
       arma::mat Z_iq_raw = Z.row_as_mat(i).t();
       arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
       //cout << "Z_iq: " << Z_iq<<"\n" << std::endl;

       arma::vec y_vec = arma::vectorise(y.row(i));
       arma::vec y_selected = y_vec.elem(data_index_iq);
       //cout << "y_selected: " << y_selected <<" \n "<< std::endl;
       arma::vec alpha_q = alpha;

       //
       arma::vec beta_intcp_vec(Iq_j, arma::fill::zeros);
       beta_intcp_vec.fill(beta_kq0(k));

       arma::vec y_tilde_iq;
       y_tilde_iq = y_selected - Z_iq * alpha_q -
         B_iq* beta_Kq - beta_intcp_vec;

       //cout<< "y_tilde_iq" << y_tilde_iq<< endl;
       //arma::mat B_iqj = B_iq; //J*(L-1)
       arma::vec B_iqj_wcoeff = trans(B_iq) * y_tilde_iq;

       arma::mat vi = trans(B_iq) * B_iq; //(L-1)*(L-1)
       Sigma_xi_sum = Sigma_xi_sum + vi;

       mu_xi_sum = mu_xi_sum + B_iqj_wcoeff;
     }
     //cout<<"finish loops for i"<< endl;
     // posterior mean and covariance matrix
     arma::mat Identity = arma::eye(L-1, L-1);
     arma::mat V_xi = arma::inv_sympd(Sigma_xi_sum*(pow(eta(k),2)/sigma2) + Identity);
     arma::vec m_kq = m.row(k).t();
     arma::vec mu_xi = V_xi * (eta(k)/sigma2*mu_xi_sum + m_kq);
     //cout<< "mu_xi"<< mu_xi <<endl;
     arma::vec xi_temp_K = arma::mvnrnd(mu_xi, V_xi); // Draw from multivariate normal distribution
     arma::rowvec xi_temp_K_row = xi_temp_K.t(); // Transpose to a row vector
     xi_kq_update.row(k) = xi_temp_K_row; // Store row vector
   }

   //update baseline
   // Initialize sum values
   arma::mat Sigma_xi_sum_K(L-1, L-1, arma::fill::zeros);
   arma::vec mu_xi_sum_K(L-1, arma::fill::zeros);

   for (int i = 0; i < I; i++){
     //cout << "i :" << i<< endl;
     int k = g_cpp(i);
     arma::uvec data_index_iq = find(data_index.row(i) == 1);
     int Iq_j = data_index_iq.size();
     arma::vec y_tilde_iq;

     arma::mat B_iq_raw  = B.subcube(i,0,1,i, J_max-1,L-1);
     arma::mat B_iq = B_iq_raw.rows(data_index_iq);
     //cout<< "B_iq"<< B_iq <<"\n"<<endl;
     arma::mat Z_iq_raw = Z.row_as_mat(i).t();
     arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
     //Rcpp::Rcout << "Z_iq: " << Z_iq<<"\n" << std::endl;

     arma::vec y_vec = arma::vectorise(y.row(i));
     arma::vec y_selected = y_vec.elem(data_index_iq);
     //Rcpp::Rcout << "y_selected: " << y_selected <<" \n "<< std::endl;
     arma::vec alpha_q = alpha;



     if (k == K-1){
       y_tilde_iq = y_selected - Z_iq*alpha_q;
     } else{
       arma::vec beta_kq = beta.row(k).t();
       arma::vec beta_intcp_vec(Iq_j, arma::fill::zeros);
       beta_intcp_vec.fill(beta_kq0(k));
       y_tilde_iq = y_selected - Z_iq*alpha_q -
         B_iq*beta_kq - beta_intcp_vec;
     }
     //cout<< "y_tilde"<< y_tilde_iq<< endl;

     //arma::mat B_iqj = B_iq;
     arma::vec B_iqj_wcoeff = trans(B_iq) * y_tilde_iq;

     arma::mat vi = trans(B_iq) * B_iq;
     Sigma_xi_sum_K = Sigma_xi_sum_K + vi;

     mu_xi_sum_K = mu_xi_sum_K + B_iqj_wcoeff;
   }
   //cout<< "finish loops for baseline"<< endl;
   // posterior mean and covariance matrix
   arma::mat Identity = arma::eye(L-1, L-1);
   arma::mat V_xi = arma::inv_sympd(Sigma_xi_sum_K*(pow(eta(K-1),2)/sigma2) + Identity);
   arma::vec m_Kq = m.row(K-1).t();
   arma::vec mu_xi = V_xi * ((eta(K-1)/sigma2)*mu_xi_sum_K + m_Kq);

   arma::vec xi_temp_K = arma::mvnrnd(mu_xi, V_xi); // Draw from multivariate normal distribution
   arma::rowvec xi_temp_K_row = xi_temp_K.t(); // Transpose to a row vector
   xi_kq_update.row(K-1) = xi_temp_K_row;

   return xi_kq_update;
 }

//' @noRd
 // [[Rcpp::export]]
 arma::vec update_beta_kq0_Q1_cpp(const arma::vec& alpha, const arma::mat& beta,
                                  const double sigma2, const arma::vec& gamma_kq0,
                                  const arma::vec& nu_kq0,
                                  const arma::mat& y, const arma::cube& Z, const arma::cube& B,
                                  const arma::vec& g_cpp, const arma::mat& data_index){

   int K = beta.n_rows;
   int L = B.n_slices;
   //int n = y.n_rows;
   int J_max = B.n_cols;
   arma::vec beta_kq0_update(K, arma::fill::zeros);

   for (int k = 0; k < K-1; ++k){

     double Sigma_beta_kq0_sum = 0;
     double mu_beta_kq0_sum = 0;
     arma::uvec index_k = find(g_cpp == k);
     arma::vec beta_Kq = beta.row(K-1).t(); // this can be put outside the for loop for i
     arma::vec beta_kq = beta.row(k).t(); // this can be put outside the for loop for i

     //Rcpp::Rcout << "k = " << k << ", q = " << q << std::endl;

     for (auto i : index_k){
       arma::uvec data_index_iq = find(data_index.row(i) == 1);
       //Rcpp::Rcout << "data_index_iq: " << data_index_iq << std::endl;
       int Iq_j = data_index_iq.size();
       arma::vec y_tilde_iq;
       arma::mat B_iq_raw = B.subcube(i,0,1,i, J_max-1,L-1);
       // B_iq_raw now is dim J_max*(L-1)
       arma::mat B_iq = B_iq_raw.rows(data_index_iq);
       arma::mat Z_iq_raw = Z.tube(i,0,i, J_max-1 );
       arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
       arma::vec y_vec = arma::vectorise(y.row(i));
       arma::vec y_selected = y_vec.elem(data_index_iq);
       arma::vec alpha_q = alpha;


       y_tilde_iq = y_selected - Z_iq * alpha_q -
         B_iq* beta_Kq - B_iq*beta_kq;

       mu_beta_kq0_sum = mu_beta_kq0_sum + sum(y_tilde_iq);
       //Rcpp::Rcout << "dot result: " << dot(B_iqj_wcoeff, xi_kq) << std::endl;

       Sigma_beta_kq0_sum = Sigma_beta_kq0_sum + Iq_j;
       //Rcpp::Rcout << "Sigma_beta_kq0_sum: " << Sigma_beta_kq0_sum << std::endl;
     }

     double sigma_beta_kq0 = 1/ (Sigma_beta_kq0_sum/sigma2 + 1/(gamma_kq0(k)*nu_kq0(k)));
     //sigma_beta_kq0 = abs(sigma_beta_kq0);
     //Rcpp::Rcout << "sigma_eta " << sigma_eta << std::endl;

     double mu_beta_kq0 = sigma_beta_kq0*mu_beta_kq0_sum/sigma2;
     beta_kq0_update(k) = arma::randn( distr_param(mu_beta_kq0,sqrt(sigma_beta_kq0)) );
     //eta_update(k,q) = R::rnorm(mu_eta, sqrt(sigma_eta));

   }
   return beta_kq0_update;
 }



//' @noRd
 // [[Rcpp::export]]
 double update_sigma2_Q1_cpp(arma::vec& alpha, arma::mat& beta,
                             const arma::mat& y, const arma::cube& Z, const arma::cube& B,
                             const arma::mat& data_index, const arma::vec& g, double h_1, double h_2){
   int K = beta.n_rows;
   //int Q = beta.n_cols;
   int I = g.n_elem;
   //int J_max = y.n_cols;

   double sigma2_update;
   // Print beta dimensions
   //Rcpp::Rcout << "beta dimensions: " << beta.n_rows << " x " << beta.n_cols << " x " << beta.n_slices << "\n";

   // Print K
   //Rcpp::Rcout << "K: " << K << "\n";



   //Rcpp::Rcout << "q : " << q << "\n";

   double y_tilde_sum = 0.0;

   int n_iq = arma::accu(data_index);

   for (int i = 0; i < I; ++i){
     int k = g(i); // index of class for i-th subject
     //Rcpp::Rcout << "i : " << i << "\n";
     //Rcpp::Rcout << "k : " << k << "\n";
     arma::uvec data_index_iq = find(data_index.row(i) == 1);
     //Rcpp::Rcout << "data_index_iq length: " << data_index_iq.n_elem << "\n";

     //int Iq_j = data_index_iq.size();

     arma::mat B_iq_raw  = B.row_as_mat(i).t();
     arma::mat B_iq = B_iq_raw.rows(data_index_iq);
     //cout<< "B_iq"<< B_iq <<"\n"<<endl;
     arma::mat Z_iq_raw = Z.row_as_mat(i).t();
     arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
     //Rcpp::Rcout << "Z_iq: " << Z_iq<<"\n" << std::endl;

     arma::vec y_vec = arma::vectorise(y.row(i));
     arma::vec y_selected = y_vec.elem(data_index_iq);
     //Rcpp::Rcout << "y_selected: " << y_selected <<" \n "<< std::endl;
     arma::vec alpha_q = alpha;


     arma::vec beta_Kq = beta.row(K-1).t();
     //Rcpp::Rcout << "beta_Kq dimensions: " << beta_Kq.n_rows << " x " << beta_Kq.n_cols << "\n";

     arma::vec y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq;
     if (k != K) {
       arma::vec beta_kq = beta.row(k).t();
       //Rcpp::Rcout << "beta_kq dimensions: " << beta_kq.n_rows << " x " << beta_kq.n_cols << "\n";
       y_tilde_iq -= B_iq * beta_kq;
     }

     y_tilde_sum += arma::accu(arma::square(y_tilde_iq));
   }

   double h_1_star = h_1 + n_iq / 2;
   double h_2_star = h_2 + y_tilde_sum / 2;

   sigma2_update = rinvgamma_rcpp(h_1_star, h_2_star);


   return sigma2_update;
 }










