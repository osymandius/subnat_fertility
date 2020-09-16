#include <TMB.hpp>	// Links in the TMB libraries

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
	PARAMETER(log_prec_spatial);

	DATA_SPARSE_MATRIX(Z_age);
	DATA_SPARSE_MATRIX(R_age);
	PARAMETER(log_prec_rw_age);
	PARAMETER_VECTOR(u_age);

	DATA_VECTOR(log_offset);
	DATA_VECTOR(births_obs);

	nll -= dnorm(beta_0, Type(0), Type(sqrt(1/0.001)), true);

	nll -= dlgamma(log_prec_spatial, Type(1), Type(20000), true);
	Type prec_spatial = exp(log_prec_spatial);
	nll -= Type(-0.5) * (u_spatial_str * (R_spatial * u_spatial_str)).sum(); 
	nll -= dnorm(u_spatial_str.sum(), Type(0), Type(0.01) * u_spatial_str.size(), 1);

	nll -= dlgamma(log_prec_rw_age, Type(1), Type(20000), true);
	Type prec_rw_age = exp(log_prec_rw_age);
	nll += GMRF(R_age)(u_age);
	nll -= dnorm(u_age.sum(), Type(0), Type(0.01) * u_age.size(), true);

	vector<Type> log_lambda(
				beta_0
				+ Z_age * u_age * sqrt(1/prec_rw_age)
				+ Z_spatial * u_spatial_str * sqrt(1/prec_spatial)
				)

	vector<Type> mu_obs_pred(M_obs * log_lambda
						+ log_offset
						);

	
	nll -= dpois(births_obs, exp(mu_obs_pred), true).sum();	

	REPORT(log_prec_rw_age);
	REPORT(log_prec_spatial);
}