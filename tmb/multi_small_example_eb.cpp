#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
  
{

  using namespace density;

  Type nll = 0;

  PARAMETER(beta_0);

  DATA_SPARSE_MATRIX(M_obs);

  DATA_SPARSE_MATRIX(R_country);
  DATA_SPARSE_MATRIX(Z_country);
  PARAMETER_VECTOR(u_country);
  PARAMETER(log_prec_country);
  
  // observations

  DATA_VECTOR(log_offset);
  DATA_VECTOR(births_obs);
  DATA_VECTOR(pop);


 
  nll -= dnorm(beta_0, Type(0), Type(sqrt(1/0.001)), true);
  
  nll -= dnorm(log_prec_country, Type(3.91), Type(0.554), true);
  Type prec_country = exp(log_prec_country); 

  nll -= Type(-0.5) * (u_country * (R_country * u_country)).sum();

  

  vector<Type> log_lambda(
                     beta_0 
                     + Z_country * u_country * sqrt(1/prec_country)
                     );

  vector<Type> mu_obs_pred(M_obs * log_lambda
                          + log_offset    
                          );

  nll -= dpois(births_obs, exp(mu_obs_pred), true).sum();  

  vector<Type> lambda(exp(log_lambda));


  vector<Type> births(lambda * pop);


  REPORT(lambda);
  

  REPORT(log_prec_country);
  REPORT(beta_0);


  return nll;
  
}
