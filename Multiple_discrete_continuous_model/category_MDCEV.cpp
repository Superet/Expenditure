#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
vec MDCEV_ll_fnC(vec param, mat e, mat e_sgn, vec y0, List X, vec M, mat gamma, vec lnJ){
	int R = e.n_cols, N = e.n_rows, nx = param.n_elem;
	
	// Utility V and indicator matrix for the original R alternatives;
	mat V = zeros(N, R), zalpha = zeros(N, R), eV = zeros(N, R);
	vec Vo = -log(y0 + 1); 
	
	for(int i=0; i<R; i++){
		mat Z 			= X[i];
		zalpha.col(i) 	= Z * param;
// 		V.col(i)		= zalpha.col(i) - log(e.col(i)/gamma.col(i) +1 );
		V.col(i)		= zalpha.col(i) - log(e.col(i) + gamma.col(i) );
		eV.col(i)	 	= exp(V.col(i)); 
	}
	
	// Sum up the log likelihood; 
// 	vec ll = lnJ + sum(z_sgn, 1) - M % log(1/y0 + sum(ez % d1, 1) ) ; 
	vec ll = Vo + sum(V % e_sgn, 1) - M % log(1/(y0+1) + sum(eV, 1) ) ; 
	return(ll);
}
