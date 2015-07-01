#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
vec policy_eval_fn(List DP_list, Function u_fn, Function transition_fn, List Q_kernal, List y_kernal, NumericVector I_grid){
	NumericMatrix state = DP_list["state"];
	NumericVector param	= DP_list["param"];
	List policy			= DP_list["policy"];
	IntegerMatrix s_idx = DP_list["state_idx"];
	vec W 				= as<vec>(DP_list["value_fn"]);
	double beta			= DP_list["beta"];
	int	ns				= state.nrow();
	vec W_new			= W;
	
	NumericMatrix Pi1	= transition_fn(state,s_idx,policy,Q_kernal,y_kernal,I_grid);	
	NumericVector u1	= u_fn(state, policy, param, Q_kernal, y_kernal);
	vec u 				= as<vec>(u1);
	mat	Pimat(Pi1.begin(), Pi1.nrow(), Pi1.ncol(), false);	
	mat	B				= eye<mat>(ns,ns) - beta*Pimat;
	if(det(B)>10e-8){
		W_new = inv(B) * u;
	}
	return(W_new);
}

// [[Rcpp::export]]
List policy_iteration_fn(List DP_list, List set_list, Function u_fn, Function transition_fn, Function policy_improve_fn, 
						List Q_kernal, List y_kernal, NumericVector I_grid, NumericVector c_grid, IntegerVector plot_state){
	double tol		= set_list["tol"];
	int max_iter	= set_list["policy_iter_max"], i=0;
	double norm_dif = tol + 1;
	NumericMatrix s_idx	= DP_list["state_idx"];
	List plot_data(plot_state.size());
	
	while(norm_dif>=tol & i<=max_iter){
		List policy			= DP_list["policy"];
		NumericVector policy_k	= policy["k"];
		NumericMatrix policy_c	= policy["c"];
		vec Ew 				= policy_eval_fn(DP_list, u_fn, transition_fn, Q_kernal, y_kernal, I_grid);
		DP_list["value_fn"] = Ew;
		List policy_sol 	= policy_improve_fn(DP_list, Q_kernal, y_kernal, I_grid, c_grid, s_idx, plot_state);
		List policy_new		= policy_sol["policy"];
		NumericVector policy_k_new	= policy_new["k"];
		NumericMatrix policy_c_new	= policy_new["c"];
		norm_dif			= max(abs(policy_c - policy_c_new)) + max(abs(policy_k - policy_k_new));
		DP_list["policy"]	= policy_new;
		plot_data			= policy_sol["plot_data"];
		i++;
	}
	DP_list["Iteration"] = i;
	return Rcpp::List::create(
		Rcpp::Named("DP_list") = DP_list,
		Rcpp::Named("plot_data") = plot_data
	);	
}