#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{

  using namespace density;

  Type nll = 0;

  PARAMETER(beta_0);

  DATA_SPARSE_MATRIX(M_obs);

  DATA_SPARSE_MATRIX(Z_spatial);
  DATA_SPARSE_MATRIX(R_spatial);

  PARAMETER_VECTOR(u_spatial_str);
  // PARAMETER(log_prec_spatial);
  
  PARAMETER_VECTOR(u_spatial_iid); 
  PARAMETER(logit_spatial_rho);
  PARAMETER(log_sigma_spatial);

  DATA_SPARSE_MATRIX(Z_interaction2);
  // PARAMETER_ARRAY(eta2);
  // PARAMETER(log_prec_eta2);
  // PARAMETER(lag_logit_eta2_phi_period);
 
  // observations

  DATA_VECTOR(log_offset);
  DATA_VECTOR(births_obs);
  // DATA_VECTOR(pop);

  // DATA_MATRIX(X_tips_dummy);
  // PARAMETER_VECTOR(beta_tips_dummy);
// DATA_SPARSE_MATRIX(Z_tips);
  // DATA_SPARSE_MATRIX(Z_period);

  // DATA_SPARSE_MATRIX(R_tips);

  // DATA_SPARSE_MATRIX(R_period);

// PARAMETER(log_prec_rw_tips);
// PARAMETER(log_prec_rw_period);
// PARAMETER_VECTOR(u_tips);
// PARAMETER_VECTOR(u_period);

 
  nll -= dnorm(beta_0, Type(0), Type(sqrt(1/0.001)), true);

  ///////////////////

  // nll -= dlgamma(log_prec_spatial, Type(1), Type(20000), true);
  // Type prec_spatial = exp(log_prec_spatial);

  // nll -= Type(-0.5) * (u_spatial_str * (R_spatial * u_spatial_str)).sum();
  
  // nll -= dnorm(u_spatial_str.sum(), Type(0), Type(0.01) * u_spatial_str.size(), 1);

  ///////////////////

  // SPATIAL
  // // ICAR
  nll -= Type(-0.5) * (u_spatial_str * (R_spatial * u_spatial_str)).sum();
  nll -= dnorm(u_spatial_str.sum(), Type(0), Type(0.01) * u_spatial_str.size(), 1);
  // // IID
  nll -= dnorm(u_spatial_iid, Type(0), Type(1), true).sum();
  // // Rho
  Type spatial_rho(exp(logit_spatial_rho)/(1+exp(logit_spatial_rho)));
  nll -= log(spatial_rho) +  log(1 - spatial_rho); // Jacobian adjustment for inverse logit'ing the parameter... 
  nll -= dbeta(spatial_rho, Type(0.5), Type(0.5), true);
  
  // // Sigma
  Type sigma_spatial = exp(log_sigma_spatial);
  nll -= dnorm(sigma_spatial, Type(0), Type(2.5), true) + log_sigma_spatial;
  
  vector<Type> spatial = sigma_spatial * (sqrt(1 - spatial_rho) * u_spatial_iid + sqrt(spatial_rho) * u_spatial_str);

  ///////////////////

  DATA_SPARSE_MATRIX(Z_age);
  DATA_SPARSE_MATRIX(R_age);
  PARAMETER(log_prec_rw_age);
  PARAMETER_VECTOR(u_age); 
 
 
  nll -= dlgamma(log_prec_rw_age, Type(1), Type(20000), true);
  Type prec_rw_age = exp(log_prec_rw_age);
  nll += GMRF(R_age)(u_age);
  nll -= dnorm(u_age.sum(), Type(0), Type(0.01) * u_age.size(), true);

  ///////////////////

  // nll -= dnorm(beta_tips_dummy, Type(0), Type(sqrt(1/0.001)), true).sum();
  // nll -= dnorm(beta_urban_dummy, Type(0), Type(1), true).sum();
  

  // nll -= dlgamma(log_prec_rw_tips, Type(1), Type(20000), true);
  // Type prec_rw_tips = exp(log_prec_rw_tips); 

  ///////////////////

  // nll -= dlgamma(log_prec_rw_period, Type(1), Type(20000), true);
  // Type prec_rw_period = exp(log_prec_rw_period);

  // nll -= Type(-0.5) * (u_period * (R_period * u_period)).sum();
  // nll -= dnorm(u_period.sum(), Type(0), Type(0.01) * u_period.size(), true);

  ///////////////////
   // ETA-2 - Space x time interaction

  // nll -= dlgamma(log_prec_eta2, Type(1), Type(20000), true);
  // Type prec_eta2 = exp(log_prec_eta2);

  // nll -= dnorm(lag_logit_eta2_phi_period, Type(0), Type(sqrt(1/0.15)), true);
  // Type eta2_phi_period = 2*exp(lag_logit_eta2_phi_period)/(1+exp(lag_logit_eta2_phi_period))-1;
  
  // nll += SEPARABLE(AR1(Type(eta2_phi_period)), GMRF(R_spatial))(eta2);

  // sum-to-zero on space x time interaction. Ensure each space effects (row) in each year (col) sum to zeo.
  // for (int i = 0; i < eta2.cols(); i++) {
    // nll -= dnorm(eta2.col(i).sum(), Type(0), Type(0.01) * eta2.col(i).size(), true);}

  // vector<Type> eta2_v(eta2);

  vector<Type> log_lambda(
                     beta_0
                     + Z_age * u_age * sqrt(1/prec_rw_age)
                     // + Z_period * u_period * sqrt(1/prec_rw_period)
                     + Z_spatial * spatial                     
                     // + Z_spatial * u_spatial_str * sqrt(1/prec_spatial)
                     // + Z_interaction2 * eta2_v * sqrt(1/prec_eta2)
                     );

  
  vector<Type> mu_obs_pred(M_obs * log_lambda
                          // + Z_tips * u_tips * sqrt(1/prec_rw_tips)  // TIPS RW
                          // + X_tips_dummy * beta_tips_dummy          // TIPS fixed effect
                          + log_offset    
                          );

    
  nll -= dpois(births_obs, exp(mu_obs_pred), true).sum();  

  vector<Type> lambda(exp(log_lambda));
  
  REPORT(lambda);
  // REPORT(log_prec_spatial);
  REPORT(log_sigma_spatial);
  REPORT(logit_spatial_rho);
  // REPORT(log_prec_eta2);
  // REPORT(eta2_phi_period);
  REPORT(log_prec_rw_age);
  // REPORT(log_prec_rw_period);
  // REPORT(log_prec_rw_tips);



  return nll;
  
}
