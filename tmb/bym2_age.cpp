#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{

  using namespace density;

  Type nll = 0;

  PARAMETER(beta_0);


  DATA_SPARSE_MATRIX(M_obs);

  DATA_SPARSE_MATRIX(Z_spatial);
  DATA_SCALAR(rankdef_R_spatial); 

  DATA_SPARSE_MATRIX(R_spatial);
  DATA_SPARSE_MATRIX(R_spatial_iid);

  PARAMETER_VECTOR(u_spatial_str);
  PARAMETER_VECTOR(u_spatial_iid); 
  PARAMETER(log_sigma_spatial);
  PARAMETER(logit_spatial_rho);

  DATA_SPARSE_MATRIX(Z_interaction3);
  PARAMETER_ARRAY(eta3a);
  PARAMETER_ARRAY(eta3b);
  PARAMETER(log_prec_eta3);
  PARAMETER(lag_logit_eta3_phi_age);

  // observations

  DATA_VECTOR(log_offset);
  DATA_VECTOR(births_obs);
  // DATA_VECTOR(pop);

 
  nll -= dnorm(beta_0, Type(0), Type(sqrt(1/0.001)), true);



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

  /////////////////

  nll -= dlgamma(log_prec_eta3, Type(1), Type(20000), true);
  Type prec_eta3 = exp(log_prec_eta3);

  nll -= dnorm(lag_logit_eta3_phi_age, Type(0), Type(sqrt(1/0.15)), true);
  Type eta3_phi_age = 2*exp(lag_logit_eta3_phi_age)/(1+exp(lag_logit_eta3_phi_age))-1;

  
  nll += SEPARABLE(AR1(Type(eta3_phi_age)), GMRF(R_spatial))(eta3a);

  // Adjust normalising constant for rank deficience of R_spatial. SEPARABLE calculates the
  // normalizing constant assuming full rank precision matrix. Add the component of the
  // constant back.
  
  Type log_det_Qar1_eta3((eta3a.cols() - 1) * log(1 - eta3_phi_age * eta3_phi_age));
  nll -= rankdef_R_spatial * 0.5 * (log_det_Qar1_eta3 - log(2 * PI));

  vector<Type> eta3a_v(eta3a);

  nll += SEPARABLE(AR1(Type(eta3_phi_age)), GMRF(R_spatial_iid))(eta3b);
  vector<Type> eta3b_v(eta3b);

  
  vector<Type> log_lambda(
                     // X_mf * beta_mf
                     beta_0
                     + Z_age * u_age * sqrt(1/prec_rw_age)
                     + Z_spatial * spatial
                     + Z_interaction3 * eta3a_v * sqrt(1/prec_eta3)
                     + Z_interaction3 * eta3b_v * sqrt(1/prec_eta3)
		                 
                     );

  
  vector<Type> mu_obs_pred(M_obs * log_lambda
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


  REPORT(log_sigma_spatial);
  REPORT(logit_spatial_rho);
  REPORT(u_spatial_str);
  REPORT(u_spatial_iid);

  REPORT(log_prec_eta3);
  REPORT(eta3_phi_age);

  REPORT(log_prec_rw_age);



  return nll;
  
}
