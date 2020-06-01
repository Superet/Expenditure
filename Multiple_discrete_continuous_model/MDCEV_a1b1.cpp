// Model versions;
// a: presence of outside good for expenditure budget(1.No vs. 2.Yes); 
// b: quantity constraint (1.No vs. 2.Yes);

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
	 
// 	vec gamma = ones<vec>(R);
	double sigma = exp(param(nx+2*R-1)); 
	List out = List::create(_["beta"]=beta, _["beta0"]=beta0, _["gamma"]=gamma, _["sigma"] = sigma) ;
	return(out);
}

// [[Rcpp::export]]
List MDCEV_ll_fnC(vec param, int nx, mat e, vec e1_index, mat p, List X_list, int base=0){
	int R = e.n_cols, N = e.n_rows;
	List param_out = param_assign(param, nx, R, base);
	vec beta = param_out["beta"], beta0 = param_out["beta0"], gamma = param_out["gamma"];
	double sigma = param_out["sigma"];
	
	// Utility V_tilde and indicator matrix for the original R alternatives;
	mat V_tilde = zeros(N, R), d_tilde = zeros(N, R);
	umat sgn_e_tilde = e > 0;
	for(int i=0; i<R; i++){
		mat Xtmp 		= X_list[i];
		d_tilde.col(i)	= e.col(i) + gamma(i)*p.col(i);
		V_tilde.col(i) 	= Xtmp * beta + beta0(i) - log(d_tilde.col(i)/gamma(i));
	}
	
	// Jacobean matrix determinant; 
	mat d_sgn = d_tilde % sgn_e_tilde;
	vec J = log(sum(d_sgn, 1)) - sum(log(d_tilde) % sgn_e_tilde, 1);
	mat sgn_e = arma::conv_to<arma::mat>::from(sgn_e_tilde); 
	vec M = sum(sgn_e, 1);
	
	// Sum up the log likelihood; 
	vec ll = J - (M - ones<vec>(N)) * log(sigma) + sum(V_tilde % sgn_e_tilde/sigma, 1) 
			- M % log(sum(exp(V_tilde/sigma), 1)); 
	
	List out = List::create(_["ll"]=ll, _["V_tilde"]=V_tilde);
	return(out);
}


// [[Rcpp::export]]
NumericVector MDCEV_LogLike_fnC(vec param, int nx, vec e1_index, arma::mat e, arma::mat p, List X_list, int base=0){
// 	int R = e.n_cols, N = e.n_rows;
	List ll_out	= MDCEV_ll_fnC(param, nx, e, e1_index, p, X_list, base);
	colvec llvec 	= ll_out["ll"];
	
	NumericVector ll = as<NumericVector>(wrap(llvec));
// 	ll.attr("gradient") = grad;
	return(ll);
}
