#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{

  using namespace density;

  Type nll = 0;

  PARAMETER(beta_0);
  // DATA_INTEGER(mics_toggle);
  // DATA_INTEGER(out_toggle);

  DATA_SPARSE_MATRIX(M_obs);



  DATA_SPARSE_MATRIX(Z_interaction3);
  PARAMETER_ARRAY(eta3);
  PARAMETER(log_prec_eta3);
  PARAMETER(lag_logit_eta3_phi_age);

 
  DATA_SPARSE_MATRIX(R_spatial);
  DATA_SPARSE_MATRIX(Z_spatial);

  // observations

  DATA_VECTOR(log_offset);
  DATA_VECTOR(births_obs);
  // DATA_VECTOR(pop);

 
  nll -= dnorm(beta_0, Type(0), Type(sqrt(1/0.001)), true);

  PARAMETER_VECTOR(u_spatial_str);
  PARAMETER_VECTOR(u_spatial_iid); 
  PARAMETER(logit_spatial_rho);
  PARAMETER(log_sigma_spatial);

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

  // ETA-3 - Space x age interaction

  nll -= dlgamma(log_prec_eta3, Type(1), Type(20000), true);
  Type prec_eta3 = exp(log_prec_eta3);

  nll -= dnorm(lag_logit_eta3_phi_age, Type(0), Type(sqrt(1/0.15)), true);
  Type eta3_phi_age = 2*exp(lag_logit_eta3_phi_age)/(1+exp(lag_logit_eta3_phi_age))-1;
  
  nll += SEPARABLE(AR1(Type(eta3_phi_age)), GMRF(R_spatial))(eta3);

  // sum-to-zero on space x time interaction. Ensure each space effects (row) in each year (col) sum to zeo.
  for (int i = 0; i < eta3.cols(); i++) {
    nll -= dnorm(eta3.col(i).sum(), Type(0), Type(0.01) * eta3.col(i).size(), true);}

  vector<Type> eta3_v(eta3);



  vector<Type> log_lambda(
                     // X_mf * beta_mf
                     beta_0
                     // + Z_interaction2 * eta2_v * sqrt(1/prec_eta2)
                     + Z_spatial * spatial
                     + Z_interaction3 * eta3_v * sqrt(1/prec_eta3)
                     );

  
  vector<Type> mu_obs_pred(M_obs * log_lambda
         // TIPS fixed effect
                          // + X_urban_dummy * beta_urban_dummy          // Urban fixed effect
                          + log_offset    
                          );

    
  nll -= dpois(births_obs, exp(mu_obs_pred), true).sum();  

  vector<Type> lambda(exp(log_lambda));
  // vector<Type> births(lambda * pop);


  // if(mics_toggle) {

  //   DATA_SPARSE_MATRIX(M_obs_mics);
    // PARAMETER_VECTOR(beta_tips_dummy_mics);

    // DATA_MATRIX(X_tips_dummy_mics);
    // DATA_SPARSE_MATRIX(Z_tips_mics);

  //   DATA_SPARSE_MATRIX(A_mics);
  //   DATA_VECTOR(log_offset_mics);
  //   DATA_VECTOR(births_obs_mics);

    // nll -= dnorm(beta_tips_dummy_mics, Type(0), Type(1), true).sum();

  //   vector<Type> births_pred_mics(A_mics * births);
  //   vector<Type> pop_mics(A_mics * pop);
  //   vector<Type> lambda_mics(births_pred_mics/pop_mics);

  //   vector<Type> mu_obs_pred_mics(M_obs_mics * log(lambda_mics) +
                                // Z_tips_mics * u_tips * sigma_rw_tips   +     // TIPS RW
                                // X_tips_dummy_mics * beta_tips_dummy_mics +          // TIPS fixed effect
  //                               log_offset_mics

  //               );

  //   nll -= dpois(births_obs_mics, exp(mu_obs_pred_mics), true).sum();  

  // }

  // if(out_toggle) {

  //   DATA_SPARSE_MATRIX(A_out);
  //   // DATA_SPARSE_MATRIX(A_out_restype);

  //   vector<Type> births_out(A_out * births);
  //   vector<Type> population_out(A_out * pop);
  //   vector<Type> lambda_out(births_out / population_out);

  //   // vector<Type> births_out_restype(A_out_restype * births);
  //   // vector<Type> population_out_restype(A_out_restype * pop);
  //   // vector<Type> lambda_out_restype(births_out_restype / population_out_restype);

  //   REPORT(lambda_out);
  //   REPORT(births_out);

  //   // REPORT(lambda_out_restype);
  //   // REPORT(births_out_restype);
  // }


    
  
  REPORT(lambda);
  // REPORT(births);
  
REPORT(log_prec_eta3);
REPORT(logit_spatial_rho);
REPORT(log_sigma_spatial);

  REPORT(eta3_phi_age);



  return nll;
  
}
