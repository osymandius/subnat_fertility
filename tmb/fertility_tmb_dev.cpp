#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{

  using namespace density;

  Type nll = 0;

  DATA_MATRIX(X_mf);
  PARAMETER_VECTOR(beta_mf);
  // PARAMETER(beta_0);

  DATA_SPARSE_MATRIX(M_obs);
  // DATA_SPARSE_MATRIX(M_obs_nat);
  
  DATA_MATRIX(X_tips_dummy);
  // DATA_MATRIX(X_tips_dummy_nat);
  PARAMETER_VECTOR(beta_tips_dummy);

  DATA_SPARSE_MATRIX(Z_tips);
  // DATA_SPARSE_MATRIX(Z_tips_nat);

  DATA_SPARSE_MATRIX(Z_age);
  DATA_SPARSE_MATRIX(Z_period);
  DATA_SPARSE_MATRIX(Z_spatial);

  // DATA_SPARSE_MATIX(Z_interaction);
  // PARAMETER_ARRAY(eta);

  // DATA_SPARSE_MATRIX(Z_interaction1);
  // PARAMETER_ARRAY(eta1);

  // DATA_SPARSE_MATRIX(Z_interaction2);
  // PARAMETER_ARRAY(eta2);

  // DATA_SPARSE_MATRIX(Z_interaction3);
  // PARAMETER_ARRAY(eta3);

  DATA_SPARSE_MATRIX(Q_tips);
  DATA_SPARSE_MATRIX(Q_age);
  DATA_SPARSE_MATRIX(Q_period);
  DATA_SPARSE_MATRIX(Q_spatial);

  // PARAMETER(log_sigma_rw_tips);
  // PARAMETER(log_sigma_rw_age);
  // PARAMETER(log_sigma_rw_period);
  // PARAMETER(log_sigma_eta1);

  PARAMETER(log_prec_rw_tips);
  PARAMETER(log_prec_rw_age);
  PARAMETER(log_prec_rw_period);
  // PARAMETER(log_prec_eta1);

  PARAMETER_VECTOR(u_tips);
  PARAMETER_VECTOR(u_age);
  PARAMETER_VECTOR(u_period);

  PARAMETER_VECTOR(u_spatial_str);
  PARAMETER_VECTOR(u_spatial_iid); 
  PARAMETER(logit_spatial_rho);
  PARAMETER(log_sigma_spatial);

  // DATA_SPARSE_MATRIX(A_national);

  DATA_SPARSE_MATRIX(A_out);
  // DATA_SPARSE_MATRIX(A_nat);
  // DATA_SPARSE_MATRIX(A_rate);
  DATA_VECTOR(pop);

  // observations

  DATA_VECTOR(log_offset);
  DATA_VECTOR(births_obs);

  // DATA_VECTOR(log_offset_nat);
  // DATA_VECTOR(births_obs_nat);

  // model
  nll -= dnorm(beta_mf, Type(0), Type(5), true).sum();

  //Fixed effect TIPS dummy
  nll -= dnorm(beta_tips_dummy, Type(0), Type(1), true).sum();

  // RW TIPS

  // Type sigma_rw_tips = exp(log_sigma_rw_tips);
  // Type tau_rw_tips(1/(sigma_rw_tips * sigma_rw_tips));
  // nll -= dnorm(sigma_rw_tips, Type(0), Type(2.5), true) + log_sigma_rw_tips;

  Type prec_rw_tips = exp(log_prec_rw_tips);
  nll -= dgamma(prec_rw_tips, Type(1), Type(50000), true) + log_prec_rw_tips;

  nll -= Type(-0.5) * (u_tips * (Q_tips * u_tips)).sum();
  nll -= dnorm(u_tips.sum(), Type(0), Type(0.01) * u_tips.size(), true);

  //RW AGE
  // Type sigma_rw_age = exp(log_sigma_rw_age);
  // Type tau_rw_age(1/(sigma_rw_age * sigma_rw_age));
  // nll -= dnorm(sigma_rw_age, Type(0), Type(2.5), true) + log_sigma_rw_age;

  Type prec_rw_age = exp(log_prec_rw_age);
  nll -= dgamma(prec_rw_age, Type(1), Type(50000), true) + log_prec_rw_age;

  nll -= Type(-0.5) * (u_age * (Q_age * u_age)).sum();
  nll -= dnorm(u_age.sum(), Type(0), Type(0.01) * u_age.size(), true);

  // RW PERIOD
  // Type sigma_rw_period = exp(log_sigma_rw_period);
  // Type tau_rw_period(1/(sigma_rw_period * sigma_rw_period));
  // nll -= dnorm(sigma_rw_period, Type(0), Type(2.5), true) + log_sigma_rw_period;

  Type prec_rw_period = exp(log_prec_rw_period);
  nll -= dgamma(prec_rw_period, Type(1), Type(50000), true) + log_prec_rw_period;

  nll -= Type(-0.5) * (u_period * (Q_period * u_period)).sum();
  nll -= dnorm(u_period.sum(), Type(0), Type(0.001) * u_period.size(), true);

  //// SPATIAL
  // ICAR
  nll -= Type(-0.5) * (u_spatial_str * (Q_spatial * u_spatial_str)).sum();
  nll -= dnorm(u_spatial_str.sum(), Type(0), Type(0.01) * u_spatial_str.size(), 1);
  
  // IID
  nll -= dnorm(u_spatial_iid, Type(0), Type(1), true).sum();
  
  // Rho
  Type spatial_rho(exp(logit_spatial_rho)/(1+exp(logit_spatial_rho)));
  nll -= log(spatial_rho) +  log(1 - spatial_rho); // Jacobian adjustment for inverse logit'ing the parameter... 
  nll -= dbeta(spatial_rho, Type(0.5), Type(0.5), true);
  
  // Sigma
  Type sigma_spatial = exp(log_sigma_spatial);
  nll -= dnorm(sigma_spatial, Type(0), Type(2.5), true) + log_sigma_spatial;
  
  vector<Type> spatial = sigma_spatial * (sqrt(1 - spatial_rho) * u_spatial_iid + sqrt(spatial_rho) * u_spatial_str);


  // nll += SEPARABLE(GMRF(Q_period), SEPARABLE(GMRF(Q_age), GMRF(Q_spatial)))(eta);
  // vector<Type> eta_v(eta);

  
  // nll += SEPARABLE(GMRF(Q_age), GMRF(Q_period))(eta1); 
  // vector<Type> eta1_v(eta1);
  // Type sigma_eta1 = exp(log_sigma_eta1);
  // nll -= dnorm(sigma_eta1, Type(0), Type(2.5), true) + log_sigma_eta1;
  // Type tau_eta1(1/(sigma_eta1 * sigma_eta1));

  // nll += SEPARABLE(GMRF(Q_period), GMRF(Q_spatial))(eta2);
  // vector<Type> eta2_v(eta2);
  // Type sigma_eta2 = exp(log_sigma_eta2);
  // nll -= dnorm(sigma_eta2, Type(0), Type(2.5), true) + log_sigma_eta2;

  // nll += SEPARABLE(GMRF(Q_age), GMRF(Q_spatial))(eta3);
  // vector<Type> eta3_v(eta3);
  // Type sigma_eta3 = exp(log_sigma_eta3);
  // nll -= dnorm(sigma_eta3, Type(0), Type(2.5), true) + log_sigma_eta3;

  vector<Type> log_lambda(
                     X_mf * beta_mf
                     + Z_age * u_age * prec_rw_age      // Age RW1
                     + Z_period * u_period * prec_rw_period
                     + Z_spatial * spatial
                     // + Z_interaction1 * eta1_v * 1/sigma_eta1
                     // + Z_interaction2 * eta2_v
                     // + Z_interaction3 * eta3_v
                     );

  
  vector<Type> mu_obs_pred(M_obs * log_lambda
                          + Z_tips * u_tips * prec_rw_tips  // TIPS RW
                          + X_tips_dummy * beta_tips_dummy          // TIPS fixed effect
                          + log_offset    
                          );

    
  nll -= dpois(births_obs, exp(mu_obs_pred), true).sum();  

  vector<Type> lambda(exp(log_lambda));
  vector<Type> births(lambda * pop);

  // vector<Type> births_pred_nat(A_nat * births);
  // vector<Type> pop_nat(A_nat * pop);
  // vector<Type> lambda_nat(births_pred_nat/pop_nat);

  // vector<Type> mu_obs_pred_nat(M_obs_nat * log(lambda_nat) +
  //                             Z_tips_nat * u_tips * 1/sigma_rw_tips   +     // TIPS RW
  //                             X_tips_dummy_nat * beta_tips_dummy +          // TIPS fixed effect
  //                             log_offset_nat

  //             );

  // nll -= dpois(births_obs_nat, exp(mu_obs_pred_nat), true).sum();  

  vector<Type> births_out(A_out * births);
  vector<Type> population_out(A_out * pop);
  vector<Type> lambda_out(births_out / population_out);

  REPORT(lambda_out);
  // ADREPORT(prec_eta1);


  return nll;
  
}
