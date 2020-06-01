#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
vec MDCEV_ll_fnC(vec param, mat e, vec y0, List X, vec M, mat d1, mat d2, vec lnJ){
	int R = e.n_cols, N = e.n_rows, nx = param.n_elem;
	
	// Utility V and indicator matrix for the original R alternatives;
	mat V = zeros(N, R), zalpha = zeros(N, R), ez = zeros(N, R);
	vec Vo = -log(y0); 
	umat sgn_e = e > 0;
	
	for(int i=0; i<R; i++){
		mat Z 			= X[i];
		zalpha.col(i) 	= Z * param;
		ez.col(i) 		= exp(zalpha.col(i)); 
		V.col(i)		= zalpha.col(i) + log(d1.col(i));
	}
	mat z_sgn = (zalpha + log(d1) ) % sgn_e;	
	
	// Sum up the log likelihood; 
// 	vec ll = lnJ + sum(z_sgn, 1) - M % log(1/y0 + sum(ez % d1, 1) ) ; 
	vec ll = sum(z_sgn, 1) - M % log(1/y0 + sum(ez % d1, 1) ) ; 
	return(ll);
}
