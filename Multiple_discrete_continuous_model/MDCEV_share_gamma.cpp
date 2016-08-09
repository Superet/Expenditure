// Model versions;

#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List MDCEV_ll_fnC(vec param, int nx, mat shr, vec y, vec s1_index, mat p, List X1, List X2, int base=0){
	int R = shr.n_cols, N = shr.n_rows, nx0 = param.n_elem;
	vec beta = param.subvec(0, nx-1), gamma0 = param.subvec(nx, nx0-2); 
	double sigma = exp(param(nx0-1));
	
	// Utility V_tilde and indicator matrix for the original R alternatives;
	mat V_tilde = zeros(N, R), d_tilde = zeros(N, R), gamma = zeros(N, R);
	umat sgn_shr_tilde = shr > 0;
	for(int i=0; i<R; i++){
		mat Xtmp1 		= X1[i];
		mat Xtmp2		= X2[i]; 
		gamma.col(i)	= exp(Xtmp2 * gamma0); 
		d_tilde.col(i)	= shr.col(i) + gamma.col(i) % p.col(i)/y;
		V_tilde.col(i) 	= Xtmp1 * beta - log(d_tilde.col(i)/gamma.col(i));
	}
	
	// Jacobean matrix determinant; 
	mat d_sgn = d_tilde % sgn_shr_tilde;
	vec J = log(sum(d_sgn, 1)) - sum(log(d_tilde) % sgn_shr_tilde, 1);
	mat sgn_shr = arma::conv_to<arma::mat>::from(sgn_shr_tilde); 
	vec M = sum(sgn_shr, 1);
	
	// Sum up the log likelihood; 
	vec ll = J - (M - ones<vec>(N)) * log(sigma) + sum(V_tilde % sgn_shr_tilde/sigma, 1) 
			- M % log(sum(exp(V_tilde/sigma), 1)) ; 
	
	List out = List::create(_["ll"]=ll, _["V_tilde"]=V_tilde);
	return(out);
}


// [[Rcpp::export]]
colvec MDCEV_LogLike_fnC(vec param, int nx, vec s1_index, arma::mat shr, vec y, arma::mat p, List X1, List X2, int base=0){
	List ll_out	= MDCEV_ll_fnC(param, nx, shr, y, s1_index, p, X1, X2, base);
	colvec llvec 	= ll_out["ll"];
	
// 	NumericVector ll = as<NumericVector>(wrap(llvec));
// 	ll.attr("gradient") = grad;
	return(llvec);
}
