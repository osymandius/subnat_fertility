#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{

  using namespace density;

  Type nll = 0;

  DATA_MATRIX(X_mf);
  PARAMETER_VECTOR(beta_mf);

  DATA_SPARSE_MATRIX(M_all_observations);
  
  DATA_MATRIX(X_tips_dummy);
  PARAMETER_VECTOR(beta_tips_dummy);

  DATA_SPARSE_MATRIX(Z_tips);
  DATA_SPARSE_MATRIX(Z_age);
  DATA_SPARSE_MATRIX(Z_period);
  DATA_SPARSE_MATRIX(Z_spatial);
  
  // DATA_SPARSE_MATRIX(Z_interaction);
  // PARAMETER_ARRAY(eta);
  // DATA_VECTOR(interaction_idx);

  DATA_SPARSE_MATRIX(Q_tips);
  DATA_SPARSE_MATRIX(Q_age);
  DATA_SPARSE_MATRIX(Q_period);
  DATA_SPARSE_MATRIX(Q_spatial);

  PARAMETER(log_sigma_rw_tips);
  PARAMETER(log_sigma_rw_age);
  PARAMETER(log_sigma_rw_period);

  PARAMETER_VECTOR(u_tips);
  PARAMETER_VECTOR(u_age);
  PARAMETER_VECTOR(u_period);

  PARAMETER_VECTOR(u_spatial_str);
  PARAMETER_VECTOR(u_spatial_iid); 
  PARAMETER(logit_spatial_rho);
  PARAMETER(log_sigma_spatial);

  // DATA_SPARSE_MATRIX(A_out);

  // observations

  DATA_VECTOR(log_offset);
  DATA_VECTOR(births_obs);

  // model
  nll -= dnorm(beta_mf, Type(0), Type(5), true).sum();

  // //Fixed effect TIPS dummy
  nll -= dnorm(beta_tips_dummy, Type(0), Type(1), true).sum();

  // RW TIPS
  Type sigma_rw_tips = exp(log_sigma_rw_tips);
  nll -= dnorm(sigma_rw_tips, Type(0), Type(2.5), true) + log_sigma_rw_tips;
  nll -= Type(-0.5) * (u_tips * (Q_tips * u_tips)).sum();
  nll -= dnorm(u_tips.sum(), Type(0), Type(0.001) * u_tips.size(), true);

  // RW AGE
  Type sigma_rw_age = exp(log_sigma_rw_age);
  nll -= dnorm(sigma_rw_age, Type(0), Type(2.5), true) + log_sigma_rw_age;
  nll -= Type(-0.5) * (u_age * (Q_age * u_age)).sum();
  nll -= dnorm(u_age.sum(), Type(0), Type(0.001) * u_age.size(), true);

  // RW PERIOD
  Type sigma_rw_period = exp(log_sigma_rw_period);
  nll -= dnorm(sigma_rw_period, Type(0), Type(2.5), true) + log_sigma_rw_period;
  nll -= Type(-0.5) * (u_period * (Q_period * u_period)).sum();
  nll -= dnorm(u_period.sum(), Type(0), Type(0.001) * u_period.size(), true);

  //// SPATIAL
  // ICAR
  nll -= Type(-0.5) * (u_spatial_str * (Q_spatial * u_spatial_str)).sum();
  nll -= dnorm(u_spatial_str.sum(), Type(0), Type(0.001) * u_spatial_str.size(), 1); // sum to zero constraint
  
  // IID
  nll -=dnorm(u_spatial_iid, Type(0), Type(2.5), true).sum();
  
  // Rho
  Type spatial_rho(exp(logit_spatial_rho)/(1+exp(logit_spatial_rho)));
  nll -= log(spatial_rho) +  log(1 - spatial_rho); // Jacobian adjustment for inverse logit'ing the parameter... 
  nll -= dbeta(spatial_rho, Type(0.5), Type(0.5), true);
  
  // Sigma
  Type sigma_spatial = exp(log_sigma_spatial);
  nll -= dnorm(sigma_spatial, Type(0), Type(2.5), true) + log_sigma_spatial;
  
  vector<Type> spatial = sqrt(1 - spatial_rho) * u_spatial_iid + sqrt(spatial_rho) * u_spatial_str;

  // nll += SEPARABLE(GMRF(Q_period), SEPARABLE(GMRF(Q_age), GMRF(Q_spatial)))(eta);

  // vector<Type> eta_v(eta);
  // // vector<Type> eta_v_clipped(interaction_idx.size());
  
  // nll -= dnorm(eta_v, Type(0), Type(1), true).sum();


  // for (int i = 0; i < interaction_idx.size() -1 ; ++i)
  // {
  //     eta_v_clipped(i) = eta_v(asDouble(interaction_idx(i)));
    
  // }

  vector<Type> mu_mf(X_mf * beta_mf);
  
  vector<Type> mu_obs_pred(M_all_observations * mu_mf +
                          X_tips_dummy * beta_tips_dummy +
                          Z_tips * u_tips * sigma_rw_tips + 
                          Z_age * u_age * sigma_rw_age + 
                          Z_period * u_period * sigma_rw_period + 
                          // Z_interaction * eta_v +
                          Z_spatial * spatial * sigma_spatial +
                          log_offset);
    
  nll -= dpois(births_obs, exp(mu_obs_pred), true).sum();

  vector<Type> omega(exp(mu_mf));

  ADREPORT(omega);

  return nll;
  
}
