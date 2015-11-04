// Model versions;

#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List param_assign(vec param, int nx, int R, int base=0){
	// Assign parameter values; 
	vec beta0_norm = zeros(1), gamma0_norm = ones<vec>(1); 
	vec beta = param.subvec(0, nx-1), beta0_sub = param.subvec(nx, nx+R-2), gamma = param.subvec(nx+R-1, nx+2*R-2);
	vec beta0 = zeros(R);
	if(base == 0){
		beta0 	= join_cols(beta0_norm, beta0_sub);
	}else{
		int j=0;
		for(int i=0; i<R; i++){
			if(i==base - 1){
				continue;
			}else{
				beta0(i) = beta0_sub(j); 
				j++; 
			}
		}
	}
	 
// 	double sigma = exp(param(nx+2*R-1)); 
	double sigma = 1; 
	List out = List::create(_["beta"]=beta, _["beta0"]=beta0, _["gamma"]=gamma, _["sigma"] = sigma) ;
	return(out);
}

// [[Rcpp::export]]
List MDCEV_ll_fnC(vec param, int nx, mat shr, vec y, vec s1_index, mat p, List X_list, int base=0){
	int R = shr.n_cols, N = shr.n_rows;
	List param_out = param_assign(param, nx, R, base);
	vec beta = param_out["beta"], beta0 = param_out["beta0"], gamma = param_out["gamma"];
	double sigma = param_out["sigma"];
	
	// Utility V_tilde and indicator matrix for the original R alternatives;
	mat V_tilde = zeros(N, R), d_tilde = zeros(N, R);
	umat sgn_shr_tilde = shr > 0;
	for(int i=0; i<R; i++){
		mat Xtmp 		= X_list[i];
		d_tilde.col(i)	= shr.col(i) + gamma(i)*p.col(i)/y;
		V_tilde.col(i) 	= Xtmp * beta + beta0(i) - log(d_tilde.col(i)/gamma(i));
	}
	
	// Jacobean matrix determinant; 
	mat d_sgn = d_tilde % sgn_shr_tilde;
	vec J = log(sum(d_sgn, 1)) - sum(log(d_tilde) % sgn_shr_tilde, 1);
	mat sgn_shr = arma::conv_to<arma::mat>::from(sgn_shr_tilde); 
	vec M = sum(sgn_shr, 1);
	
	// Sum up the log likelihood; 
	vec ll = J - (M - ones<vec>(N)) * log(sigma) + sum(V_tilde % sgn_shr_tilde/sigma, 1) 
			- M % log(sum(exp(V_tilde/sigma), 1)); 
	
	List out = List::create(_["ll"]=ll, _["V_tilde"]=V_tilde);
	return(out);
}


// [[Rcpp::export]]
colvec MDCEV_LogLike_fnC(vec param, int nx, vec s1_index, arma::mat shr, vec y, arma::mat p, List X_list, int base=0){
	List ll_out	= MDCEV_ll_fnC(param, nx, shr, y, s1_index, p, X_list, base);
	colvec llvec 	= ll_out["ll"];
	
// 	NumericVector ll = as<NumericVector>(wrap(llvec));
// 	ll.attr("gradient") = grad;
	return(llvec);
}
