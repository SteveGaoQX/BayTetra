//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
// # include "utils_cpp.h"

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
arma::mat update_alpha_cpp(const arma::cube& beta, const arma::mat& omega, const arma::vec& sigma2,
                           const arma::cube& data_index, const arma::cube& y, const arma::field<arma::cube>& B,
                           const arma::field<arma::cube>& Z, const arma::vec& g, const arma::cube& Z_sum,
                           const arma::mat& V_alpha_inv, const arma::vec& V_alpha_inv_mu_alpha){
  // Update the alpha
  // args: 1: beta is K*Q*L at t-1 time
  //       2: omega is the I*Q mat at t-1 time
  //       3: sigma2 is the length Q col vect at t-1 time
  // returns: the updates alpha


  // Get dimensions from the data
  int K = beta.n_rows;
  int Q = beta.n_cols;
  int S = Z_sum.n_cols; // Z_sum in dim Q*S*S
  int I = g.n_elem;
  int J_max = data_index.n_slices; // data_index in dim I*Q*J_max
  //cout << "K:" << K << "  Q: "<< Q <<endl;
  arma::mat alpha_update(Q, S, arma::fill::none);
  arma::mat V_n(S, S);
  arma::vec mu_n(S);

  for (int q = 0; q < Q; ++q){

    arma::vec Zy_sum(S, arma::fill::zeros);

    for (int i = 0; i < I; i++){

      int k = g(i);
      //cout<< "i:  "<< i <<" k: "<< k<<endl;
      arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1);
      arma::vec omega_vec(data_index_iq.size(), arma::fill::zeros);
      omega_vec.fill(omega(i, q));
      arma::mat B_iq_raw = B(i).tube(q,0,q, J_max-1 );
      //cout<< "B_iq_raw: " << B_iq_raw <<endl;
      arma::mat B_iq = B_iq_raw.rows(data_index_iq);
      arma::vec beta_kq = beta.tube(k, q);
      //cout<< "beta_kq" << beta_kq <<endl;
      arma::vec beta_Kq = beta.tube(K-1, q);
      arma::vec y_vec = arma::vectorise(y.tube(i, q));
      arma::vec y_selected = y_vec.elem(data_index_iq);
      arma::vec y_tilde_iq;
      if (k == K-1){ // baseline
        y_tilde_iq = y_selected - B_iq * beta_kq - omega_vec;
      } else { // not baseline
        y_tilde_iq = y_selected - B_iq * beta_Kq -
          B_iq * beta_kq - omega_vec;
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
arma::mat update_eta_kq_cpp(const arma::mat& alpha, const arma::cube& beta, const arma::mat& omega,
                            const arma::vec& sigma2, const arma::mat& beta_kq0,
                            const arma::cube& xi, const arma::mat& gamma_kq, const arma::mat& nu_kq,
                            const arma::cube& y, const arma::field<arma::cube>& Z,
                            const arma::field<arma::cube>& B, const arma::vec& g_cpp, const arma::cube& data_index){
  //cout << "updating eta in cpp" << endl;
  int K = beta.n_rows;
  int Q = beta.n_cols;
  int L = B(0).n_slices;
  int n = y.n_rows;
  int J_max = B(0).n_cols;
  arma::mat eta_update(K, Q, arma::fill::zeros);
  try {
    for (int k = 0; k < K-1; ++k){
      for (int q = 0; q < Q; ++q){
        double Sigma_eta_sum = 0;
        double mu_eta_sum = 0;
        arma::uvec index_k = find(g_cpp == k);
        arma::vec beta_Kq = beta.tube(K-1, q); // this can be put outside the for loop for i
        //cout<< "Current k:"<< k <<"q:"<< q <<"\n"<<"\n"<<endl;

        for (auto i : index_k){
          //Rcpp::Rcout << "i = " << i << std::endl;

          arma::uvec data_index_iq = find(data_index.tube(i,q) == 1);
          //Rcpp::Rcout << "data_index_iq: " << data_index_iq << std::endl;

          int Iq_j = data_index_iq.size();
          arma::vec y_tilde_iq;

          arma::mat B_iq_raw = B(i).subcube(q,0,1,q, J_max-1,L-1);
          //B_iq_raw now is dim J_max*(L-1)
          arma::mat B_iq = B_iq_raw.rows(data_index_iq);
          //cout<< "B_iq" << B_iq << "\n" << endl;
          arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
          arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
          //Rcpp::Rcout << "Z_iq: " << Z_iq <<"\n" << std::endl;

          arma::vec y_vec = arma::vectorise(y.tube(i, q));
          arma::vec y_selected = y_vec.elem(data_index_iq);
          //Rcpp::Rcout << "y_selected: " << y_selected <<"\n"<< std::endl;
          //cout<< "B_iq* beta_Kq" << B_iq* beta_Kq <<endl;
          arma::vec alpha_q = alpha.row(q).t();
          //cout<< "Z_iq * alpha_q: "<< Z_iq * alpha_q <<"\n" <<endl;
          arma::vec omega_vec(Iq_j, arma::fill::zeros);
          omega_vec.fill(omega(i, q));
          //cout<< "omega_iq" << omega(i,q) <<" \n "<< endl;
          arma::vec beta_intcp_vec(Iq_j, arma::fill::zeros);
          beta_intcp_vec.fill(beta_kq0(k,q));
          //cout<< "beta_kq0(k,q)" << beta_kq0(k,q) <<"\n"<< endl;

          y_tilde_iq = y_selected - Z_iq * alpha_q -
            B_iq* beta_Kq - beta_intcp_vec -
            omega_vec;

          //Rcpp::Rcout << "y_tilde_iq: " << y_tilde_iq <<"\n"<< std::endl;
          arma::mat B_iqj = B_iq;
          //Rcpp::Rcout << "B_iqj: " << B_iqj << std::endl;
          arma::vec B_iqj_wcoeff = trans(B_iqj) * y_tilde_iq;
          //Rcpp::Rcout << "B_iqj_wcoeff: " << B_iqj_wcoeff << std::endl;

          arma::vec xi_kq = xi.tube(k,q);
          //Rcpp::Rcout << "xi_kq: " << xi_kq << std::endl;

          mu_eta_sum = mu_eta_sum + dot(B_iqj_wcoeff, xi_kq);
          //Rcpp::Rcout << "dot result: " << dot(B_iqj_wcoeff, xi_kq) << std::endl;

          arma::vec vi = B_iqj * xi_kq;
          //Rcpp::Rcout << "vi: " << vi << std::endl;

          Sigma_eta_sum = Sigma_eta_sum + dot(vi,vi);
          //Rcpp::Rcout << "Sigma_eta_sum: " << Sigma_eta_sum << std::endl;
        }
        //Rcpp::Rcout << "Check for final " << std::endl;
        //Rcpp::Rcout << "Sigma_eta_sum: " << Sigma_eta_sum << std::endl;
        //Rcpp::Rcout << "mu_eta_sum: " << mu_eta_sum << std::endl;
        //Rcpp::Rcout << "sigma2(q) " << sigma2(q) << std::endl;
        //Rcpp::Rcout << "gamma_kq(k,q) " << gamma_kq(k,q) << std::endl;
        //Rcpp::Rcout << "nu_kq(k,q) " << nu_kq(k,q) << std::endl;

        double sigma_eta = 1/ (Sigma_eta_sum/sigma2(q) + 1/(gamma_kq(k,q)*nu_kq(k,q)));
        //sigma_eta = abs(sigma_eta);
        //Rcpp::Rcout << "sigma_eta " << sigma_eta << std::endl;
        //cout << "mu_eta_sum" << mu_eta_sum << endl;
        double mu_eta = sigma_eta*mu_eta_sum/sigma2(q);
        //cout<< "current k:" << k + 1 << "current q:" << q + 1 << endl;
        //cout << "mean is" << mu_eta <<endl;
        //cout<< "variance is" << sigma_eta << endl;
        //eta_update(k,q) = arma::randn( distr_param(mu_eta,sqrt(sigma_eta)) );
        //eta_update(k,q) = mu_eta;
        //cout<< "eta_updata(k,q)" << eta_update(k,q) << endl;
        eta_update(k,q) = Rcpp::rnorm( 1, mu_eta, sqrt(sigma_eta) )[0];
        //eta_update(k,q) = R::rnorm(mu_eta, sqrt(sigma_eta));
      }
    }
  } catch (const std::runtime_error& e) {
    // Catch a std::runtime_error, which is thrown by Armadillo for errors such as dimension mismatches
    Rcpp::Rcout << "Caught a runtime_error exception on eta_first_loops: " << e.what() << std::endl;
  } catch (const std::exception& e) {
    // Catch any other standard exceptions
    Rcpp::Rcout << "Caught an exception on on eta_first_loops: " << e.what() << std::endl;
  } catch (...) {
    // Catch all other exceptions
    Rcpp::Rcout << "Caught an unknown exception on on eta_first_loops!" << std::endl;
  }



  //cout<< "begin update K for eta" << endl;


  try {
    for (int q = 0; q < Q; ++q){
      double Sigma_eta_sum_K = 0;
      double mu_eta_sum_K = 0;
      for (int i = 0; i < n; ++i){
        int k = g_cpp(i);
        arma::uvec data_index_iq = find(data_index.tube(i,q) == 1);
        int Iq_j = data_index_iq.size();

        arma::vec y_tilde_iq;

        arma::mat B_iq_raw = B(i).subcube(q,0,1,q, J_max-1,L-1);
        // B_iq_raw now is dim J_max*(L-1)
        arma::mat B_iq = B_iq_raw.rows(data_index_iq);
        arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
        arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
        arma::vec y_vec = arma::vectorise(y.tube(i, q));
        arma::vec y_selected = y_vec.elem(data_index_iq);
        arma::vec alpha_q = alpha.row(q).t();
        arma::vec omega_vec(Iq_j, arma::fill::zeros);
        omega_vec.fill(omega(i, q));

        arma::vec beta_intcp_vec(Iq_j, arma::fill::zeros);
        beta_intcp_vec.fill(beta_kq0(k,q));


        arma::vec beta_kq = beta.tube(k,q);

        if (k == K-1){
          y_tilde_iq = y_selected - Z_iq*alpha_q - omega_vec;
          //if((i> n-21) && q == 1){
          //  Rcpp::Rcout << "i: " << i << "\n";
          //  Rcpp::Rcout << "k: " << k << "\n";
          //  Rcpp::Rcout << "y_selected: " << y_selected << "\n";
          //  Rcpp::Rcout << "Z_iq*alpha_q: " << Z_iq*alpha_q << "\n";
          //  Rcpp::Rcout << "omega_vec: " << omega_vec << "\n";
          //  Rcpp::Rcout << "y_tilde_iq: " << y_tilde_iq << "\n";
          //
          //}
        } else{

          y_tilde_iq = y_selected - Z_iq*alpha_q - B_iq*beta_kq - beta_intcp_vec - omega_vec;
          //if((i> n-21) && q == 1){
          //  Rcpp::Rcout << "i: " << i << "\n";
          //  Rcpp::Rcout << "k: " << k << "\n";
          //  Rcpp::Rcout << "y_tilde_iq: " << y_tilde_iq << "\n";
          //}
        }

        arma::mat B_iqj = B_iq;
        arma::vec B_iqj_wcoeff = trans(B_iqj) * y_tilde_iq;
        arma::vec xi_Kq = xi.tube(K-1,q);

        mu_eta_sum_K = mu_eta_sum_K + dot(B_iqj_wcoeff,xi_Kq);


        arma::vec vi = B_iqj*xi_Kq;
        Sigma_eta_sum_K = Sigma_eta_sum_K + dot(vi,vi);


        //if((i == n-2 || i == n-1) && q == 1){
        //  Rcpp::Rcout << "i: " << i << "\n";
        //  Rcpp::Rcout << "k: " << k << "\n";
        //  Rcpp::Rcout << "y_tilde_iq: " << y_tilde_iq << "\n";
        //}
      }
      //if (q == 1) {
      //  Rcpp::Rcout << "Sigma_eta_sum_K"<<Sigma_eta_sum_K << std::endl;
      //  cout<< "mu_eta_sum_K" << mu_eta_sum_K << endl;
      //}

      //Rcpp::Rcout << "trans_scalar " << trans_scalar << std::endl;

      double sigma_eta = 1/ (Sigma_eta_sum_K/sigma2(q) + 1/(gamma_kq(K-1,q)*nu_kq(K-1,q)));
      //sigma_eta = abs(sigma_eta);

      double mu_eta = sigma_eta*mu_eta_sum_K/sigma2(q);
      //cout << "current k:" << K << "q:" << q + 1 << endl;
      //cout << "mean is:" << mu_eta << endl;
      //cout << "variance is" << sigma_eta << endl;
      //eta_update(K-1,q) = arma::randn( distr_param(mu_eta,sqrt(sigma_eta)) );
      eta_update(K-1,q) = Rcpp::rnorm( 1, mu_eta, sqrt(sigma_eta) )[0];

      //eta_update(K-1,q) = R::rnorm(mu_eta, sqrt(sigma_eta));
    }
  } catch (const std::runtime_error& e) {
    // Catch a std::runtime_error, which is thrown by Armadillo for errors such as dimension mismatches
    Rcpp::Rcout << "Caught a runtime_error exception on eta_first_loops: " << e.what() << std::endl;
  }
  //cout << "finish updating eta in cpp" << endl;
  return eta_update;
  //gc();
}

//' @noRd
//[[Rcpp::export]]
arma::cube update_xi_kq_cpp(const arma::mat& alpha, const arma::cube& beta, const arma::mat& omega,
                            const arma::vec& sigma2, const arma::mat& beta_kq0, const arma::mat& eta,
                            const arma::cube& m,
                            const arma::vec& g_cpp, const arma::cube& data_index,
                            const arma::cube& y, const arma::field<arma::cube>& Z, const arma::field<arma::cube>& B){
  //cout<< "begin update xi in cpp"<< endl;
  int K = beta.n_rows;
  int Q = beta.n_cols;
  int L = B(0).n_slices;
  int I = y.n_rows;
  int J_max = B(0).n_cols;

  arma::cube xi_kq_update(K, Q, L-1, arma::fill::zeros);


  // update xi_kq, k = 1,2,...,K-1
  for (int k = 0; k < K-1; k++){
    for (int q = 0; q < Q; q++){

      arma::vec beta_Kq = beta.tube(K-1, q); // this can be put outside the for loop for i

      // Indexes for subject in k-th class
      arma::uvec index_k = arma::find(g_cpp == k);
      // Initialize sum values
      arma::mat Sigma_xi_sum(L-1, L-1, arma::fill::zeros);
      arma::vec mu_xi_sum(L-1, arma::fill::zeros);
      for (auto i : index_k){
        // subject index for i-th subject at q-th response
        //cout << "i:" << i << "k: " << k <<"q: "<< q <<endl;
        arma::uvec data_index_iq = arma::find(data_index.tube(i,q) == 1);
        int Iq_j = data_index_iq.size();

        arma::mat B_iq_raw = B(i).subcube(q,0,1,q, J_max-1,L-1);
        arma::mat B_iq = B_iq_raw.rows(data_index_iq);

        arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
        arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);

        arma::vec y_vec = arma::vectorise(y.tube(i, q));
        arma::vec y_selected = y_vec.elem(data_index_iq);
        arma::vec alpha_q = alpha.row(q).t();

        arma::vec omega_vec(Iq_j, arma::fill::zeros);
        omega_vec.fill(omega(i, q));

        arma::vec beta_intcp_vec(Iq_j, arma::fill::zeros);
        beta_intcp_vec.fill(beta_kq0(k,q));

        arma::vec y_tilde_iq;
        y_tilde_iq = y_selected - Z_iq * alpha_q -
          B_iq* beta_Kq - beta_intcp_vec -
          omega_vec;

        arma::mat B_iqj = B_iq; //J*(L-1)
        arma::vec B_iqj_wcoeff = trans(B_iqj) * y_tilde_iq;

        arma::mat vi = trans(B_iqj) * B_iqj; //(L-1)*(L-1)
        Sigma_xi_sum = Sigma_xi_sum + vi;

        mu_xi_sum = mu_xi_sum + B_iqj_wcoeff;
      }
      //cout<< "current k:" <<k << "q:" <<q <<endl;
      //cout<< "Sigma_xi_sum" << Sigma_xi_sum << endl;
      //cout<< "eta(k,q)" << eta(k,q) <<endl;

      // posterior mean and covariance matrix
      arma::mat Identity = arma::eye(L-1, L-1);
      arma::mat mat = Sigma_xi_sum*(pow(eta(k,q),2)/sigma2(q)) + Identity;

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
      arma::mat V_xi = arma::inv_sympd(mat);
      arma::vec m_kq = m.tube(k,q);
      arma::vec mu_xi = V_xi * ((eta(k,q)/sigma2(q))*mu_xi_sum + m_kq);
      //cout<< "m_kq" <<m_kq <<endl;

      //Rcpp::Rcout << "this is k: " << k + 1 << " q: " << q << "\n";
      //Rcpp::Rcout << "the mean is: " << mu_xi << "\n";
      //Rcpp::Rcout << "the variance is: " << V_xi << "\n";
      //arma::vec xi_temp_K = arma::mvnrnd(mu_xi, V_xi); // Draw from multivariate normal distribution
      //arma::rowvec xi_temp_K_row = xi_temp_K.t(); // Transpose to a row vector

      xi_kq_update.tube(k,q) = rmvn_rcpp(1, mu_xi, V_xi).row(0); // Store row vector
    }
  }




  try {
    // update xi_Kq
    for (int q = 0; q < Q; q++){
      // Initialize sum values
      arma::mat Sigma_xi_sum_K(L-1, L-1, arma::fill::zeros);
      arma::vec mu_xi_sum_K(L-1, arma::fill::zeros);

      for (int i = 0; i < I; i++){
        int k = g_cpp(i);
        arma::uvec data_index_iq = find(data_index.tube(i,q) == 1);
        int Iq_j = data_index_iq.size();

        arma::vec y_tilde_iq;

        arma::mat B_iq_raw = B(i).subcube(q,0,1,q, J_max-1,L-1);
        // B_iq_raw now is dim J_max*(L-1)
        arma::mat B_iq = B_iq_raw.rows(data_index_iq);
        arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
        arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
        arma::vec y_vec = arma::vectorise(y.tube(i, q));
        arma::vec y_selected = y_vec.elem(data_index_iq);
        arma::vec alpha_q = alpha.row(q).t();
        arma::vec omega_vec(Iq_j, arma::fill::zeros);
        omega_vec.fill(omega(i, q));


        if (k == K-1){
          y_tilde_iq = y_selected - Z_iq*alpha_q - omega_vec;
        } else{
          arma::vec beta_kq = beta.tube(k,q);
          arma::vec beta_intcp_vec(Iq_j, arma::fill::zeros);
          beta_intcp_vec.fill(beta_kq0(k,q));

          y_tilde_iq = y_selected - Z_iq*alpha_q -
            B_iq*beta_kq - beta_intcp_vec - omega_vec;
        }


        arma::mat B_iqj = B_iq;
        arma::vec B_iqj_wcoeff = trans(B_iqj) * y_tilde_iq;

        arma::mat vi = trans(B_iqj) * B_iqj;
        Sigma_xi_sum_K = Sigma_xi_sum_K + vi;

        mu_xi_sum_K = mu_xi_sum_K + B_iqj_wcoeff;
      }
      //cout << "current Q:" << q <<"\n"<<endl;
      //cout << "current Sigma_xi_sum_K: " << Sigma_xi_sum_K <<endl;
      //cout << "current eta(K-1,q):" << eta(K-1,q) << endl;
      // posterior mean and covariance matrix
      arma::mat Identity = arma::eye(L-1, L-1);
      arma::mat V_xi = arma::inv_sympd(Sigma_xi_sum_K*(pow(eta(K-1,q),2)/sigma2(q)) + Identity);
      arma::vec m_Kq = m.tube(K-1,q);
      arma::vec mu_xi = V_xi * ((eta(K-1,q)/sigma2(q))*mu_xi_sum_K + m_Kq);

      //Rcpp::Rcout << "this is k: " << K << " q: " << q+1 << "\n";
      //Rcpp::Rcout << "the mean is: " << mu_xi << "\n";
      //Rcpp::Rcout << "the variance is: " << V_xi << "\n";
      //arma::vec xi_temp_K = arma::mvnrnd(mu_xi, V_xi); // Draw from multivariate normal distribution
      //arma::rowvec xi_temp_K_row = xi_temp_K.t(); // Transpose to a row vector
      xi_kq_update.tube(K-1,q) = rmvn_rcpp(1, mu_xi, V_xi).row(0);

      //xi_kq_update.tube(K-1,q) = mu_xi;

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
  return xi_kq_update;
  //gc();
}
//' @noRd
 // [[Rcpp::export]]
arma::mat update_beta_kq0_cpp(const arma::mat& alpha, const arma::cube& beta, const arma::mat& omega, const arma::vec& sigma2,
                              const arma::mat& gamma_kq0, const arma::mat& nu_kq0,
                              const arma::cube& y, const arma::field<arma::cube>& Z, const arma::field<arma::cube>& B, const arma::vec& g_cpp, const arma::cube& data_index){

  int K = beta.n_rows;
  int Q = beta.n_cols;
  int L = B(0).n_slices;
  //int n = y.n_rows;
  int J_max = B(0).n_cols;
  arma::mat beta_kq0_update(K, Q, arma::fill::zeros);

  for (int k = 0; k < K-1; ++k){
    for (int q = 0; q < Q; ++q){
      double Sigma_beta_kq0_sum = 0;
      double mu_beta_kq0_sum = 0;
      arma::uvec index_k = find(g_cpp == k);
      arma::vec beta_Kq = beta.tube(K-1, q); // this can be put outside the for loop for i
      arma::vec beta_kq = beta.tube(k, q); // this can be put outside the for loop for i

      //Rcpp::Rcout << "k = " << k << ", q = " << q << std::endl;

      for (auto i : index_k){
        arma::uvec data_index_iq = find(data_index.tube(i,q) == 1);
        //Rcpp::Rcout << "data_index_iq: " << data_index_iq << std::endl;
        int Iq_j = data_index_iq.size();
        arma::vec y_tilde_iq;
        arma::mat B_iq_raw = B(i).subcube(q,0,1,q, J_max-1,L-1);
        // B_iq_raw now is dim J_max*(L-1)
        arma::mat B_iq = B_iq_raw.rows(data_index_iq);
        arma::mat Z_iq_raw = Z(i).tube(q,0,q, J_max-1 );
        arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
        arma::vec y_vec = arma::vectorise(y.tube(i, q));
        arma::vec y_selected = y_vec.elem(data_index_iq);
        arma::vec alpha_q = alpha.row(q).t();

        arma::vec omega_vec(Iq_j, arma::fill::zeros);
        omega_vec.fill(omega(i, q));

        y_tilde_iq = y_selected - Z_iq * alpha_q -
          B_iq* beta_Kq - B_iq*beta_kq - omega_vec;

        mu_beta_kq0_sum = mu_beta_kq0_sum + sum(y_tilde_iq);
        //Rcpp::Rcout << "dot result: " << dot(B_iqj_wcoeff, xi_kq) << std::endl;

        Sigma_beta_kq0_sum = Sigma_beta_kq0_sum + Iq_j;
        //Rcpp::Rcout << "Sigma_beta_kq0_sum: " << Sigma_beta_kq0_sum << std::endl;
      }

      double sigma_beta_kq0 = 1/ (Sigma_beta_kq0_sum/sigma2(q) + 1/(gamma_kq0(k,q)*nu_kq0(k,q)));
      //sigma_beta_kq0 = abs(sigma_beta_kq0);
      //Rcpp::Rcout << "sigma_eta " << sigma_eta << std::endl;

      double mu_beta_kq0 = sigma_beta_kq0*mu_beta_kq0_sum/sigma2(q);
      //beta_kq0_update(k,q) = arma::randn( distr_param(mu_beta_kq0,sqrt(sigma_beta_kq0)) );
      beta_kq0_update(k,q) = Rcpp::rnorm( 1, mu_beta_kq0, sqrt(sigma_beta_kq0) )[0];
      //eta_update(k,q) = R::rnorm(mu_eta, sqrt(sigma_eta));
    }
  }
  return beta_kq0_update;
}




// [[Rcpp::export]]
double logpost_omega_cpp(int i, const arma::rowvec& omega_i, const arma::mat& alpha, const arma::cube& beta, const arma::vec& sigma2, const arma::mat& Sigma_omega,
                         const arma::cube& y, const arma::field<arma::cube>& Z, const arma::field<arma::cube>& B, const arma::vec& g,
                         const arma::cube& data_index, int K, int Q, int J_max){
  // calculate the log-posterior for omega_i
  // args: index of i-th subject, i from 0 to I-1 since it call by other Rcpp function
  //    omega_i: omega values for i-subjects, as a rowvector length Q
  double logpost = 0.0;

  // calculate log-likelihood
  int k = g(i) ; // index of class for i-th subject
  //Rcpp::Rcout << "k =  " << k << "\n";
  for (int q = 0; q < Q; ++q){
    //Rcpp::Rcout << "q =  " << q << "\n";
    arma::vec alpha_q = alpha.row(q).t();
    //Rcpp::Rcout << "alpha_q dimensions: " << alpha_q.n_rows << " x " << alpha_q.n_cols << "\n";

    arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1); // index of data for i-th subject at q-th response
    //Rcpp::Rcout << "data_index_iq elements: " << data_index_iq.n_elem << "\n";

    arma::vec y_vec = arma::vectorise(y.tube(i, q));
    //Rcpp::Rcout << "y_vec dimensions: " << y_vec.n_rows << " x " << y_vec.n_cols << "\n";

    arma::vec y_selected = y_vec.elem(data_index_iq);
    //Rcpp::Rcout << "y_selected dimensions: " << y_selected.n_rows << " x " << y_selected.n_cols << "\n";

    arma::mat Z_iq_raw = Z(i).tube(q, 0, q, J_max-1);
    //Rcpp::Rcout << "Z_iq_raw dimensions: " << Z_iq_raw.n_rows << " x " << Z_iq_raw.n_cols << "\n";

    arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
    //Rcpp::Rcout << "Z_iq dimensions: " << Z_iq.n_rows << " x " << Z_iq.n_cols << "\n";

    arma::mat B_iq_raw = B(i).tube(q, 0, q, J_max-1);
    //Rcpp::Rcout << "B_iq_raw dimensions: " << B_iq_raw.n_rows << " x " << B_iq_raw.n_cols << "\n";

    arma::mat B_iq = B_iq_raw.rows(data_index_iq);
    //Rcpp::Rcout << "B_iq dimensions: " << B_iq.n_rows << " x " << B_iq.n_cols << "\n";
    arma::vec y_tilde_iq;

    arma::vec beta_Kq = beta.tube(K-1,q);
    arma::vec beta_kq;
    //Rcpp::Rcout << "beta_Kq dimensions: " << beta_Kq.n_rows << " x " << beta_Kq.n_cols << "\n";

    if(k != K-1){
      beta_kq = beta.tube(k,q);
      //Rcpp::Rcout << "beta_kq dimensions: " << beta_kq.n_rows << " x " << beta_kq.n_cols << "\n";

    }

    if (k == K - 1){
      // if i-th subject in baseline class
      y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq;
    } else {
      // if i-th subject not in baseline class
      y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq - B_iq * beta_kq;
    }
    //Rcpp::Rcout << "y_tilde_iq dimensions: " << y_tilde_iq.n_rows << " x " << y_tilde_iq.n_cols << "\n";
    arma::rowvec y_tilde_iq_rowvec = y_tilde_iq.t();

    int J_iq = data_index_iq.size();
    //Rcpp::Rcout << "J_iq:" << J_iq << std::endl;

    arma::vec mu(J_iq, arma::fill::ones);
    //Rcpp::Rcout << "omega_ilength:" << omega_i.n_cols << std::endl;

    mu *= omega_i(q);
    //Rcpp::Rcout << "mu:" << mu << std::endl;
    arma::rowvec mu_rowvec = mu.t();
    //Rcpp::Rcout << "mu_rowvec dimension: " << mu_rowvec.n_cols << std::endl;


    arma::mat sigma_mat = arma::eye<arma::mat>(J_iq, J_iq) * sigma2(q);
    //Rcpp::Rcout << "sigma_mat dimension: " << sigma_mat.n_rows << " x " << sigma_mat.n_cols << std::endl;

    // print sigma_mat content
    //Rcpp::Rcout << "sigma_mat content:\n" << sigma_mat << std::endl;

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
arma::mat update_omega_cpp(const arma::mat& alpha, const arma::cube& beta, const arma::mat& omega, const arma::vec& sigma2, const arma::mat& Sigma_omega,
                           const arma::cube& y, const arma::field<arma::cube>& Z, const arma::field<arma::cube>& B, const arma::vec& g,
                           const arma::cube& data_index,const double omega_var_update = 0.04){

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
    double ratio = logpost_omega_cpp(i, omega_i_new_row, alpha, beta, sigma2, Sigma_omega, y, Z, B, g, data_index, K, Q, J_max) -
      logpost_omega_cpp(i, omega_i_row, alpha, beta, sigma2, Sigma_omega, y, Z, B, g, data_index, K, Q, J_max);

    // accept or reject
    if (log(R::runif(0, 1)) < ratio){
      omega_update.row(i) = omega_i_new.row(0);
    }
  }

  return(omega_update);
}







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



//' @importFrom Rcpp sourceCpp
//' @noRd


// [[Rcpp::export]]
arma::vec update_sigma2_cpp(arma::mat& alpha, arma::cube& beta, arma::mat& omega,
                            const arma::cube& y, const arma::field<arma::cube>& Z, const arma::field<arma::cube>& B,
                            const arma::cube& data_index, const arma::vec& g, double h_1, double h_2){
  int K = beta.n_rows;
  int Q = beta.n_cols;
  int I = g.n_elem;
  int J_max = Z(0).n_cols;

  arma::vec sigma2_update(Q, arma::fill::zeros);
  // Print beta dimensions
  //Rcpp::Rcout << "beta dimensions: " << beta.n_rows << " x " << beta.n_cols << " x " << beta.n_slices << "\n";

  // Print K
  //Rcpp::Rcout << "K: " << K << "\n";


  for (int q = 0; q < Q; ++q){
    //Rcpp::Rcout << "q : " << q << "\n";

    double y_tilde_sum = 0.0;

    arma::mat data_index_q = data_index.tube(0, q, I - 1, q);
    int n_iq = arma::accu(data_index_q);

    for (int i = 0; i < I; ++i){
      int k = g(i); // index of class for i-th subject
      //Rcpp::Rcout << "i : " << i << "\n";
      //Rcpp::Rcout << "k : " << k << "\n";
      arma::uvec data_index_iq = arma::find(data_index.tube(i, q) == 1); // index of data for i-th subject at q-th response
      //Rcpp::Rcout << "data_index_iq length: " << data_index_iq.n_elem << "\n";

      arma::mat B_iq_raw = B(i).tube(q, 0, q, J_max-1);
      //Rcpp::Rcout << "B_iq_raw dimensions: " << B_iq_raw.n_rows << " x " << B_iq_raw.n_cols << "\n";

      arma::mat B_iq = B_iq_raw.rows(data_index_iq);
      //Rcpp::Rcout << "B_iq dimensions: " << B_iq.n_rows << " x " << B_iq.n_cols << "\n";

      arma::mat Z_iq_raw = Z(i).tube(q, 0, q, J_max-1);
      //Rcpp::Rcout << "Z_iq_raw dimensions: " << Z_iq_raw.n_rows << " x " << Z_iq_raw.n_cols << "\n";

      arma::mat Z_iq = Z_iq_raw.rows(data_index_iq);
      //Rcpp::Rcout << "Z_iq dimensions: " << Z_iq.n_rows << " x " << Z_iq.n_cols << "\n";

      arma::vec y_vec = arma::vectorise(y.tube(i, q));
      //Rcpp::Rcout << "y_vec dimensions: " << y_vec.n_rows << " x " << y_vec.n_cols << "\n";

      arma::vec y_selected = y_vec.elem(data_index_iq);
      //Rcpp::Rcout << "y_selected dimensions: " << y_selected.n_rows << " x " << y_selected.n_cols << "\n";

      arma::vec alpha_q = alpha.row(q).t();
      //Rcpp::Rcout << "alpha_q dimensions: " << alpha_q.n_rows << " x " << alpha_q.n_cols << "\n";

      arma::vec omega_vec(data_index_iq.size(), arma::fill::zeros);
      omega_vec.fill(omega(i, q));
      //Rcpp::Rcout << "omega_vec dimensions: " << omega_vec.n_rows << " x " << omega_vec.n_cols << "\n";

      arma::vec beta_Kq = beta.tube(K-1, q);
      //Rcpp::Rcout << "beta_Kq dimensions: " << beta_Kq.n_rows << " x " << beta_Kq.n_cols << "\n";

      arma::vec y_tilde_iq = y_selected - Z_iq * alpha_q - B_iq * beta_Kq - omega_vec;
      if (k != K) {
        arma::vec beta_kq = beta.tube(k, q);
        //Rcpp::Rcout << "beta_kq dimensions: " << beta_kq.n_rows << " x " << beta_kq.n_cols << "\n";
        y_tilde_iq -= B_iq * beta_kq;
      }

      y_tilde_sum += arma::accu(arma::square(y_tilde_iq));
    }

    double h_1_star = h_1 + n_iq / 2;
    double h_2_star = h_2 + y_tilde_sum / 2;

    sigma2_update(q) = rinvgamma_rcpp(h_1_star, h_2_star);
    // sigma2_update(q) = 1/randg( distr_param(h_1_star, h_2_star) );
  }

  return sigma2_update;
}

























