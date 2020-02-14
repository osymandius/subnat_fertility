#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{
  Type nll = 0;
  
  DATA_ARRAY(dat_m)
  DATA_ARRAY(log_offset)
  DATA_SPARSE_MATRIX(Q);
  DATA_SPARSE_MATRIX(Z_i);
  
  PARAMETER(alpha);
  // PARAMETER_VECTOR(alpha_a);
  // PARAMETER_VECTOR(beta_a);
  PARAMETER_VECTOR(rw_age);
  PARAMETER_VECTOR(rw_time);
  // PARAMETER_MATRIX(rw_interaction);
  PARAMETER(log_prec_epsilon);
  PARAMETER(log_prec_rw_age);
  PARAMETER(log_prec_rw_time);
  // PARAMETER(log_prec_rw_interaction);
  
  // PARAMETER_VECTOR(V);
  // PARAMETER_VECTOR(W);
  
  
  //////// SPATIAL MODEL ////////
  
  PARAMETER_VECTOR(logit_U_str); // ICAR
  vector<Type> U_str = exp(logit_U_str)/(1+exp(logit_U_str));
  
  nll -= Type(-0.5) * (U_str * (Q * U_str)).sum(); // I don't understand this line. Pairwise? But how?
  nll -= dnorm(U_str.sum(), Type(0), Type(0.001) * U_str.size(), 1); // sum to zero constraint
  
  PARAMETER_VECTOR(U_iid); // iid
  nll -=dnorm(U_iid, Type(0), Type(1), true).sum();
  
  PARAMETER(logit_U_rho);
  Type U_rho(exp(logit_U_rho)/(1+exp(logit_U_rho)));
  nll -= log(U_rho) +  log(1 - U_rho); // Jacobian adjustment for inverse logit'ing the parameter... 
  nll -= dbeta(U_rho, Type(0.5), Type(0.5), true);

  
  PARAMETER(log_U_sigma);
  Type U_sigma = exp(log_U_sigma);
  nll -= dnorm(U_sigma, Type(0), Type(2.5), true) + log_U_sigma;
  
  vector<Type> spatial = U_sigma*(sqrt(1 - U_rho) * U_iid + sqrt(U_rho) * U_str);
  vector<Type> spatial_array = Z_i * spatial; 
  
  // PARAMETER(log_sigma2_V);
  // PARAMETER(log_sigma2_U);
  
  int num_time = 22;
  int num_age = 7;
  int num_dist = 33;
  
  
  // Type tau_V = 1 / exp(log_sigma2_V);
  // Type tau_U = 1 / exp(log_sigma2_U);
  // 
  // nll -= dgamma(tau_V, Type(0.5), Type(2000), 1);
  // nll -= dgamma(tau_U, Type(0.5), Type(2000), 1);
  // 
  // 
  // 
  // nll -= dnorm(V, Type(0), exp(0.5 * log_sigma2_V), 1).sum();
  // // This is the pairwise thing I think??
  // vector<Type> tmp = P * W;
  // nll -= -0.5 * (W * tmp).sum();
  // vector<Type> U = W * exp(0.5 * log_sigma2_U);
  // nll -= dnorm(U.sum(), Type(0), Type(0.00001), 1); // sum to zero constraint
  
  
  array<Type> log_mu(num_time, num_age, num_dist);
  
  nll -= dnorm(alpha, Type(0), Type(1), TRUE);
  // nll -= dnorm(alpha_a, Type(0), Type(1), TRUE).sum();
  // nll -= dnorm(beta_a, Type(0), Type(1), TRUE).sum();
  nll -= dnorm(log_prec_epsilon, Type(0), Type(1), TRUE);
  nll -= dnorm(log_prec_rw_age, Type(0), Type(1), TRUE);
  nll -= dnorm(log_prec_rw_time, Type(0), Type(1), TRUE);
  
  
  // alpha_a(a) + beta_a(a)*t + 
  for (int t = 0; t < num_time; t++) {
    for(int a = 0; a < num_age; a++) {
      for(int i = 0; i < num_dist; i++) {
        // log_mu(t, a, i) = alpha + rw_time(t) + rw_age(a) + V(i) + U(i);
        log_mu(t, a, i) = alpha + rw_time(t) + rw_age(a) + U_sigma*(
          sqrt(1 - U_rho) * U_iid(i) + 
          sqrt(U_rho) * U_str(i)
        );
      }
    }
  }
  
  for (int a = 0; a < num_age; a++) {
    nll -= dnorm(rw_age(a), a==0 ? 0 : rw_age(a-1), exp(log_prec_rw_age), TRUE);
  }
  
  
  for (int t = 0; t < num_time; t++) {
    nll -= dnorm(rw_time(t), t==0 ? 0 : rw_time(t-1), exp(log_prec_rw_time), TRUE);
  }

  // Is this necessary?
  vector<Type> mu_col(num_time);
  for (int y = 0; y<num_age; y++) {
    mu_col = log_mu.col(y);
    nll -= dnorm(mu_col, 0, exp(log_prec_epsilon), TRUE).sum();
  }

  
  for (int t = 0; t < num_time; t++) {
    for(int a = 0; a < num_age; a++) {
      for(int i = 0; i < num_dist; i++) {
        if(log_offset(t,a,i) > -100) {
          nll -= dpois(dat_m(t, a, i), exp(log_mu(t, a, i) + log_offset(t, a, i)), TRUE);
        }
      }
    }
  }
  
  ADREPORT(log_mu);
  
  return nll;
  
}
