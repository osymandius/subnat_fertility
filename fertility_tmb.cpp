#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{
  
  DATA_ARRAY(dat_m)
  DATA_ARRAY(log_offset)
  DATA_SPARSE_MATRIX(P);
  
  PARAMETER_VECTOR(alpha_a);
  PARAMETER_VECTOR(beta_a);
  PARAMETER_VECTOR(rw_age);
  PARAMETER_VECTOR(rw_time);
  // PARAMETER_MATRIX(rw_interaction);
  PARAMETER(log_prec_epsilon);
  PARAMETER(log_prec_rw_age);
  PARAMETER(log_prec_rw_time);
  // PARAMETER(log_prec_rw_interaction);
  
  PARAMETER_VECTOR(V);
  PARAMETER_VECTOR(W);
  
  PARAMETER(log_sigma2_V);
  PARAMETER(log_sigma2_U);
  
  Type nll = 0;
  
  
  int num_time = 22;
  int num_age = 7;
  int num_dist = 33;
  
  
  Type tau_V = 1 / exp(log_sigma2_V);
  Type tau_U = 1 / exp(log_sigma2_U);
  
  nll -= dgamma(tau_V, Type(0.5), Type(2000), 1);
  nll -= dgamma(tau_U, Type(0.5), Type(2000), 1);
  
  
  // nll -= dnorm(tau_V, Type(0), Type(10), 1);
  // nll -= dnorm(tau_U, Type(0), Type(10), 1);
  
  nll -= dnorm(V, Type(0), exp(0.5 * log_sigma2_V), 1).sum();
  
  
  // This is the pairwise thing I think??
  vector<Type> tmp = P * W;
  nll -= -0.5 * (W * tmp).sum();
  
  
  vector<Type> U = W * exp(0.5 * log_sigma2_U);
  nll -= dnorm(U.sum(), Type(0), Type(0.00001), 1);
  
  
  array<Type> log_mu(num_time, num_age, num_dist);
  
  nll -= dnorm(alpha_a, Type(0), Type(1), TRUE).sum();
  nll -= dnorm(beta_a, Type(0), Type(1), TRUE).sum();
  nll -= dnorm(log_prec_epsilon, Type(0), Type(1), TRUE);
  nll -= dnorm(log_prec_rw_age, Type(0), Type(1), TRUE);
  nll -= dnorm(log_prec_rw_time, Type(0), Type(1), TRUE);
  
  // nll -= dnorm(log_prec_rw_interaction, Type(0), Type(1), TRUE);
  
  // nll -= dlgamma(log_prec_epsilon, Type(1), Type(0.00005), TRUE);
  // nll -= dlgamma(log_prec_nu, Type(1), Type(0.00005), 1); 
  
  for (int t = 0; t < num_time; t++) {
    for(int a = 0; a < num_age; a++) {
      for(int i = 0; i < num_dist; i++) {
        log_mu(t, a, i) = alpha_a(a) + beta_a(a)*t + rw_time(t) + rw_age(a) + V(i) + U(i);
      }
    }
  }
  
  for (int a = 0; a < num_age; a++) {
    nll -= dnorm(rw_age(a), a==0 ? 0 : rw_age(a-1), exp(log_prec_rw_age), TRUE);
  }
  
  
  for (int t = 0; t < num_time; t++) {
    nll -= dnorm(rw_time(t), t==0 ? 0 : rw_time(t-1), exp(log_prec_rw_time), TRUE);
  }
  
  
  // for (int t = 0; t < num_time; t++) {
  //   for (int a = 0; a < num_age; a++) {
  //     nll -= dnorm(rw_interaction(t, a), t==0 ? 0 : rw_interaction(t-1, a), exp(log_prec_rw_interaction), TRUE);
  //   }
  // }
  
  
  
  vector<Type> mu_col(num_time);
  
  for (int y = 0; y<num_age; y++) {
    mu_col = log_mu.col(y);
    nll -= dnorm(mu_col, 0, exp(log_prec_epsilon), TRUE).sum();
  }
  
  
  for (int t = 0; t < num_time; t++) {
    for(int a = 0; a < num_age; a++) {
      for(int i = 0; i < num_dist; i++) {
        std::isnan(log_offset(t,a,i)) ? nll-=0 : nll -= dpois(dat_m(t, a, i), exp(log_mu(t, a, i) + log_offset(t, a, i)), TRUE);
      }
    }
  }
  
  
  REPORT(V);
  REPORT(W);
  REPORT(log_sigma2_V);
  REPORT(log_sigma2_U);
  
  return nll;
  
}
