// Model versions;

#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
vec MDCEV_ll_fnC(vec param, int nx, mat shr, vec y, mat p, List X1, List X2, mat shr_sgn, vec M){
	int R = shr.n_cols, N = shr.n_rows, nx0 = param.n_elem;
	vec beta = param.subvec(0, nx-1), gamma0 = param.subvec(nx, nx0-2); 
	double sigma = exp(param(nx0-1));
	
	// Utility V_tilde and indicator matrix for the original R alternatives;
	mat V_tilde = zeros(N, R), d_tilde = zeros(N, R), gamma = zeros(N, R);
	for(int i=0; i<R; i++){
		mat Xtmp1 		= X1[i];
		mat Xtmp2		= X2[i]; 
		gamma.col(i)	= exp(Xtmp2 * gamma0); 
		d_tilde.col(i)	= shr.col(i) + gamma.col(i) % p.col(i)/y;
		V_tilde.col(i) 	= Xtmp1 * beta - log(d_tilde.col(i)/gamma.col(i));
	}
	
	// Jacobean matrix determinant; 
	mat d_sgn = d_tilde % shr_sgn;
	vec J = log(sum(d_sgn, 1)) - sum(log(d_tilde) % shr_sgn, 1);
	
	// Sum up the log likelihood; 
	vec ll = J - (M - ones<vec>(N)) * log(sigma) + sum(V_tilde % shr_sgn/sigma, 1) 
			- M % log(sum(exp(V_tilde/sigma), 1)) ; 
	
	return(ll);
}
