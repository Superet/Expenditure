#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List MDCP_ll_fnC(vec param, int nx, arma::vec y, arma::mat e, arma::mat p, List X_list){
	int R = e.n_cols, N = e.n_rows;
	vec beta = param.subvec(0,nx-1), beta_int = param.subvec(nx, nx+R-1), gamma = param.subvec(nx+R, nx+2*R - 1);
	double psi_R1 = exp(param(nx+2*R)), psi_R2 = exp(param(nx+2*R+1));
	// double sigma = 1.0;
	double sigma = exp(param(nx+2*R+2));

	vec e_R1 = y - sum(e, 1);
	umat sgn_e = e>0;
	mat matone(N, R, fill::ones);
	mat Lambda = zeros(N,R), d = zeros(N,R), V = zeros(N,R), f = zeros(N,R), F = zeros(N,R);

	for(int i=0; i<R; i++){
		mat X = X_list[i];
		Lambda.col(i)	= psi_R1/e_R1 ;
		d.col(i)		= e.col(i) + gamma(i)*p.col(i);
		V.col(i) 		= log(d.col(i)/gamma(i)) - X*beta - beta_int(i) + log(Lambda.col(i) );
		f.col(i)		= exp(-0.5 * square(V.col(i))/pow(sigma,2)) / (sqrt(2*PI)*sigma);
		NumericVector tmp = as<NumericVector>(wrap(vectorise(V.col(i))));
		NumericVector tmp1 = pnorm(tmp, 0.0, sigma);
		F.col(i)		= as<vec>(tmp1);
	}
	mat d_sgn	= d % sgn_e;
	vec J 		= sum(d_sgn, 1) - sum(log(d) % sgn_e, 1); 
	vec ll		= J + sum(log(f) % sgn_e, 1) + sum(log(F) % (matone - sgn_e), 1);
	
	List out = List::create(_["ll"]=ll, _["V"]=V, _["f"]=f, _["F"]=F);
	return(out);
}


// [[Rcpp::export]]
NumericVector MDCP_LogLike_fnC(vec param, int nx, arma::vec y, arma::mat e, arma::mat p, List X_list){
	List ll_out = MDCP_ll_fnC(param, nx, y, e, p, X_list);
	vec llvec = ll_out["ll"];
	NumericVector ll = as<NumericVector>(wrap(llvec));
// 	ll.attr("gradient") = grad;
	return(ll);
}

