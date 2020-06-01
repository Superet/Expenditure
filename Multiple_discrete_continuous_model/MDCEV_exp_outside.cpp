// Model versions;

#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List MDCEV_ll_fnC(vec param, int nx, mat e, vec e1, mat p, List X1, List X2, vec M, mat pe_idx, vec d1){
	int R = e.n_cols, N = e.n_rows, nx0 = param.n_elem;
	vec beta = param.subvec(0, nx-1), gamma0 = param.subvec(nx, nx0-2); 
	double sigma = exp(param(nx0-1));
	
	// Utility V_tilde and indicator matrix for the original R alternatives;
	mat V_tilde = zeros(N, R+1), d_tilde = zeros(N, R+1), gamma = zeros(N, R+1);
// 	umat sgn_shr_tilde = shr > 0;
	for(int i=0; i<R; i++){
		mat Xtmp1 		= X1[i];
		mat Xtmp2		= X2[i]; 
		gamma.col(i)	= exp(Xtmp2 * gamma0); 
		d_tilde.col(i)	= e.col(i) + gamma.col(i) % p.col(i);
		V_tilde.col(i) 	= Xtmp1 * beta - log(e.col(i)/(gamma.col(i) % p.col(i)) + ones<vec>(N)) - log(p.col(i));
	}
	// Outside good; 
	gamma.col(R)		= ones<vec>(N); 
	d_tilde.col(R) 		= d1; 	
	V_tilde.col(R)		= -log(e1 + ones<vec>(N)); 
	
	
	// Jacobean matrix determinant; 
	mat d_sgn = d_tilde % pe_idx;
	vec J = log(sum(d_sgn, 1)) - sum(log(d_tilde) % pe_idx, 1);
	
	// Sum up the log likelihood; 
	vec ll = J - (M - ones<vec>(N)) * log(sigma) + sum(V_tilde % pe_idx/sigma, 1) 
			- M % log(sum(exp(V_tilde/sigma), 1)) ; 
	
	List out = List::create(_["ll"]=ll, _["V_tilde"]=V_tilde);
	return(out);
}
