#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{

  using namespace density;

  Type nll = 0;

  PARAMETER(beta_0);

  DATA_SPARSE_MATRIX(M_obs);

  DATA_MATRIX(X_tips_dummy);
  PARAMETER_VECTOR(beta_tips_dummy);

  DATA_SPARSE_MATRIX(Z_tips);
  DATA_SPARSE_MATRIX(R_tips);
  PARAMETER(log_prec_rw_tips);
  PARAMETER_VECTOR(u_tips);

  DATA_SPARSE_MATRIX(Z_spatial);
  DATA_SPARSE_MATRIX(R_spatial);
  DATA_SCALAR(rankdef_R_spatial); // rank deficiency of the R_spatial structure matrix
  
  PARAMETER_VECTOR(u_spatial_str);
  
  DATA_SPARSE_MATRIX(Z_interaction1);
  PARAMETER_ARRAY(eta1);
  PARAMETER(log_prec_eta1);
  PARAMETER(lag_logit_eta1_phi_age);
  PARAMETER(lag_logit_eta1_phi_period);

  DATA_SPARSE_MATRIX(Z_interaction2);
  PARAMETER_ARRAY(eta2);
  PARAMETER(log_prec_eta2);
  PARAMETER(lag_logit_eta2_phi_period);

  DATA_SPARSE_MATRIX(Z_interaction3);
  PARAMETER_ARRAY(eta3);
  PARAMETER(log_prec_eta3);
  PARAMETER(lag_logit_eta3_phi_age);
 
  // observations

  DATA_VECTOR(log_offset);
  DATA_VECTOR(births_obs);
  DATA_VECTOR(pop);
  DATA_INTEGER(mics_toggle);
  DATA_INTEGER(out_toggle);



  // DATA_MATRIX(X_urban_dummy);
  // PARAMETER_VECTOR(beta_urban_dummy);

 
  nll -= dnorm(beta_0, Type(0), Type(sqrt(1/0.001)), true);

  ///////////////////

  nll -= dnorm(beta_tips_dummy, Type(0), Type(sqrt(1/0.001)), true).sum();
  
  nll -= dlgamma(log_prec_rw_tips, Type(1), Type(20000), true);
  Type prec_rw_tips = exp(log_prec_rw_tips); 

  nll -= Type(-0.5) * (u_tips * (R_tips * u_tips)).sum();
  nll -= dnorm(u_tips.sum(), Type(0), Type(0.01) * u_tips.size(), true);

  /////////////////

  // nll -= dnorm(beta_urban_dummy, Type(0), Type(sqrt(1/0.001)), true).sum();

  ///////////////////

  PARAMETER_VECTOR(u_spatial_iid); 
  PARAMETER(log_sigma_spatial);
  PARAMETER(logit_spatial_rho);

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

  // PARAMETER(lag_logit_phi_age);
  // nll -= dnorm(lag_logit_phi_age, Type(0), Type(sqrt(1/0.15)), true);
  // Type phi_age = 2*exp(lag_logit_phi_age)/(1+exp(lag_logit_phi_age))-1;
  // nll += AR1(Type(phi_age))(u_age);

  ///////////////////

  DATA_SPARSE_MATRIX(Z_period);
  DATA_SPARSE_MATRIX(R_period);
  PARAMETER(log_prec_rw_period);
  PARAMETER_VECTOR(u_period); 
 
 
  nll -= dlgamma(log_prec_rw_period, Type(1), Type(20000), true);
  Type prec_rw_period = exp(log_prec_rw_period);

  nll -= Type(-0.5) * (u_period * (R_period * u_period)).sum();
  nll -= dnorm(u_period.sum(), Type(0), Type(0.01) * u_period.size(), true);

  // PARAMETER_VECTOR(lag_logit_ar2_phi_period);
  // PARAMETER(lag_logit_phi_period);
  // nll += AR1(Type(phi_period))(u_period);

  // nll -= dnorm(lag_logit_phi_period, Type(0), Type(sqrt(1/0.15)), true);
  // Type phi_period = 2*exp(lag_logit_phi_period)/(1+exp(lag_logit_phi_period))-1;

  // nll -= dnorm(lag_logit_ar2_phi_period, Type(0), Type(sqrt(1/0.15)), true).sum();
  // vector<Type> ar2_phi_period = 2*exp(lag_logit_ar2_phi_period)/(1+exp(lag_logit_ar2_phi_period)) -1;

  // ARk_t<Type> ar2_period(ar2_phi_period);

  // nll += ar2_period(u_period);

  ////////////////////
  // ETA-1 - Age x time interaction

  nll -= dlgamma(log_prec_eta1, Type(1), Type(20000), true);
  Type prec_eta1 = exp(log_prec_eta1);

  nll -= dnorm(lag_logit_eta1_phi_age, Type(0), Type(sqrt(1/0.15)), true);
  Type eta1_phi_age = 2*exp(lag_logit_eta1_phi_age)/(1+exp(lag_logit_eta1_phi_age))-1;

  nll -= dnorm(lag_logit_eta1_phi_period, Type(0), Type(sqrt(1/0.15)), true);
  Type eta1_phi_period = 2*exp(lag_logit_eta1_phi_period)/(1+exp(lag_logit_eta1_phi_period))-1;
  
  nll += SEPARABLE(AR1(Type(eta1_phi_age)), AR1(Type(eta1_phi_period)))(eta1);
  vector<Type> eta1_v(eta1);

  ///////////////////
   // ETA-2 - Space x time interaction

  nll -= dlgamma(log_prec_eta2, Type(1), Type(20000), true);
  Type prec_eta2 = exp(log_prec_eta2);

  nll -= dnorm(lag_logit_eta2_phi_period, Type(0), Type(sqrt(1/0.15)), true);
  Type eta2_phi_period = 2*exp(lag_logit_eta2_phi_period)/(1+exp(lag_logit_eta2_phi_period))-1;
  
  nll += SEPARABLE(AR1(Type(eta2_phi_period)), GMRF(R_spatial))(eta2);

  // Adjust normalising constant for rank deficience of R_spatial. SEPARABLE calculates the
  // normalizing constant assuming full rank precision matrix. Add the component of the
  // constant back.
  
  Type log_det_Qar1_eta2((eta2.cols() - 1) * log(1 - eta2_phi_period * eta2_phi_period));
  nll -= rankdef_R_spatial * 0.5 * (log_det_Qar1_eta2 - log(2 * PI));

  // sum-to-zero on space x time interaction. Ensure each space effects (row) in each year (col) sum to zeo.
  for (int i = 0; i < eta2.cols(); i++) {
    nll -= dnorm(eta2.col(i).sum(), Type(0), Type(0.01) * eta2.col(i).size(), true);}
  
  // array<Type> c; // eta2
  // c = eta2.transpose();

  // for (int i = 0; i < c.cols(); i++) { //eta2
  //   nll -= dnorm(c.col(i).sum(), Type(0), Type(0.01) * c.col(i).size(), true);} // eta2

  vector<Type> eta2_v(eta2);

  ////////////////////

  nll -= dlgamma(log_prec_eta3, Type(1), Type(20000), true);
  Type prec_eta3 = exp(log_prec_eta3);

  nll -= dnorm(lag_logit_eta3_phi_age, Type(0), Type(sqrt(1/0.15)), true);
  Type eta3_phi_age = 2*exp(lag_logit_eta3_phi_age)/(1+exp(lag_logit_eta3_phi_age))-1;

  
  nll += SEPARABLE(AR1(Type(eta3_phi_age)), GMRF(R_spatial))(eta3);

    // Adjust normalising constant for rank deficience of R_spatial. SEPARABLE calculates the
  // normalizing constant assuming full rank precision matrix. Add the component of the
  // constant back.
  
  Type log_det_Qar1_eta3((eta3.cols() - 1) * log(1 - eta3_phi_age * eta3_phi_age));
  nll -= rankdef_R_spatial * 0.5 * (log_det_Qar1_eta3 - log(2 * PI));

  // // sum-to-zero on space x age interaction. Ensure each space effects (row) in each year (col) sum to zeo.
  for (int i = 0; i < eta3.cols(); i++) {
    nll -= dnorm(eta3.col(i).sum(), Type(0), Type(0.01) * eta3.col(i).size(), true);}

  // array<Type> b; //eta3
  // b = eta3.transpose();

  // for (int i = 0; i < b.cols(); i++) { //eta3
  //   nll -= dnorm(b.col(i).sum(), Type(0), Type(0.01) * b.col(i).size(), true);} //eta3

  vector<Type> eta3_v(eta3);

  ///////////////////////

  vector<Type> log_lambda(
                     beta_0
                     + Z_age * u_age * sqrt(1/prec_rw_age)
                     + Z_period * u_period * sqrt(1/prec_rw_period)
                     + Z_spatial * spatial                     
                     // + Z_spatial * u_spatial_str * sqrt(1/prec_spatial)
                     + Z_interaction1 * eta1_v * sqrt(1/prec_eta1)
                     + Z_interaction2 * eta2_v * sqrt(1/prec_eta2)
                     + Z_interaction3 * eta3_v * sqrt(1/prec_eta3)
                     );

  
  vector<Type> mu_obs_pred(M_obs * log_lambda
                          + Z_tips * u_tips * sqrt(1/prec_rw_tips)  // TIPS RW
                          + X_tips_dummy * beta_tips_dummy          // TIPS fixed effect
                          + log_offset    
                          );

    
  nll -= dpois(births_obs, exp(mu_obs_pred), true).sum();  

  vector<Type> lambda(exp(log_lambda));

  vector<Type> births(lambda * pop);


  if(mics_toggle) {

    DATA_SPARSE_MATRIX(M_obs_mics);
    PARAMETER_VECTOR(beta_tips_dummy_mics);

    DATA_MATRIX(X_tips_dummy_mics);
    DATA_SPARSE_MATRIX(Z_tips_mics);

    DATA_SPARSE_MATRIX(A_mics);
    DATA_VECTOR(log_offset_mics);
    DATA_VECTOR(births_obs_mics);

    nll -= dnorm(beta_tips_dummy_mics, Type(0), Type(1), true).sum();

    vector<Type> births_pred_mics(A_mics * births);
    vector<Type> pop_mics(A_mics * pop);
    vector<Type> lambda_mics(births_pred_mics/pop_mics);

    vector<Type> mu_obs_pred_mics(M_obs_mics * log(lambda_mics) +
                                Z_tips_mics * u_tips * sqrt(1/prec_rw_tips)   +     // TIPS RW
                                X_tips_dummy_mics * beta_tips_dummy_mics +          // TIPS fixed effect
                                log_offset_mics

                );

    nll -= dpois(births_obs_mics, exp(mu_obs_pred_mics), true).sum();  

  }

  if(out_toggle) {

    DATA_SPARSE_MATRIX(A_out);
    // DATA_SPARSE_MATRIX(A_out_restype);

    vector<Type> births_out(A_out * births);
    vector<Type> population_out(A_out * pop);
    vector<Type> lambda_out(births_out / population_out);

    // vector<Type> births_out_restype(A_out_restype * births);
    // vector<Type> population_out_restype(A_out_restype * pop);
    // vector<Type> lambda_out_restype(births_out_restype / population_out_restype);

    REPORT(lambda_out);
    REPORT(births_out);

    // REPORT(lambda_out_restype);
    // REPORT(births_out_restype);
  }


  REPORT(log_sigma_spatial);
  REPORT(logit_spatial_rho);

  REPORT(log_prec_eta1);
  REPORT(eta1_phi_age);
  REPORT(eta1_phi_period);

  REPORT(log_prec_eta2);
  REPORT(eta2_phi_period);

  REPORT(log_prec_eta3);
  REPORT(eta3_phi_age);

  REPORT(log_prec_rw_age);
  REPORT(log_prec_rw_period);
  REPORT(log_prec_rw_tips);

  REPORT(beta_tips_dummy);
  // REPORT(beta_urban_dummy);

  REPORT(eta2);
  REPORT(eta3);



  return nll;
  
}
