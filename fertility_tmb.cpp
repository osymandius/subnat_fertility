#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{
  DATA_MATRIX(dat_m)
  DATA_MATRIX(log_offset)
  
  PARAMETER_VECTOR(alpha_a);
  PARAMETER_VECTOR(beta_a);
  PARAMETER_VECTOR(rw_age);
  PARAMETER_VECTOR(rw_time);
  PARAMETER_MATRIX(rw_interaction);
  PARAMETER(log_prec_epsilon);
  PARAMETER(log_prec_rw_age);
  PARAMETER(log_prec_rw_time);
  PARAMETER(log_prec_rw_interaction);
  
  
  int num_time = dat_m.rows();
  int num_age = dat_m.cols();
  
  Type nll = 0;
  matrix<Type> log_mu(num_time, num_age);
  
  nll -= dnorm(alpha_a, Type(0), Type(1), TRUE).sum();
  nll -= dnorm(beta_a, Type(0), Type(1), TRUE).sum();
  nll -= dnorm(log_prec_epsilon, Type(0), Type(1), TRUE);
  nll -= dnorm(log_prec_rw_age, Type(0), Type(1), TRUE);
  nll -= dnorm(log_prec_rw_time, Type(0), Type(1), TRUE);
  nll -= dnorm(log_prec_rw_interaction, Type(0), Type(1), TRUE);
  
  // nll -= dlgamma(log_prec_epsilon, Type(1), Type(0.00005), TRUE);
  // nll -= dlgamma(log_prec_nu, Type(1), Type(0.00005), 1); 
  
  for (int t = 0; t < num_time; t++) {
    for(int a = 0; a < num_age; a++) {
      log_mu(t, a) = alpha_a(a) + beta_a(a)*t + rw_time(t) + rw_age(a) + rw_interaction(t, a);
    }
  }
  
  for (int a = 0; a < num_age; a++) {
    nll -= dnorm(rw_age(a), a==0 ? 0 : rw_age(a-1), exp(log_prec_rw_age), TRUE);
  }
  
  for (int t = 0; t < num_time; t++) {
    nll -= dnorm(rw_time(t), t==0 ? 0 : rw_time(t-1), exp(log_prec_rw_time), TRUE);
  }
  
  for (int t = 0; t < num_time; t++) {
    for (int a = 0; a < num_age; a++) {
      nll -= dnorm(rw_interaction(t, a), t==0 ? 0 : rw_interaction(t-1, a), exp(log_prec_rw_interaction), TRUE);
    }
  }
  
  
  
  vector<Type> mu_col(num_time);
  
  for (int y = 0; y<num_age; y++) {
    mu_col = log_mu.col(y);
    nll -= dnorm(mu_col, 0, exp(log_prec_epsilon), TRUE).sum();
  }
  
  for (int t = 0; t < num_time; t++) {
    for(int a = 0; a < num_age; a++) {
      nll -= dpois(dat_m(t, a), exp(log_mu(t, a) + log_offset(t, a)), TRUE);
    }
  }
  
  REPORT(log_mu);
  
  return nll;
  
}
