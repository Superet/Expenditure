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
List MDCEV_ll_fnC(vec param, int nx,
double c_q, arma::vec y, arma::vec Q, arma::mat e, arma::mat p, List X_list){
	int R = e.n_cols, N = e.n_rows;
	vec beta = param.subvec(0,nx-1), beta_int = param.subvec(nx, nx+R-1), gamma = param.subvec(nx+R, nx+2*R - 1);
	double psi_R1 = exp(param(nx+2*R)), psi_R2 = exp(param(nx+2*R+1));
	// double sigma = 1.0;
	double sigma = exp(param(nx+2*R+2));

	vec e_R1 = y - sum(e, 1);
	vec e_R2 = Q - sum(e/p,1);
	umat sgn_e = e>0;
	mat matone(N, R, fill::ones);
	mat Lambda = zeros(N,R), d = zeros(N,R), V = zeros(N,R);
	vec lambda1 = zeros(N), lambda2 = zeros(N);
	
	lambda1 = psi_R1/square(e_R1);
	lambda2 = psi_R2/square(e_R2 + c_q);
	for(int i=0; i<R; i++){
		mat X = X_list[i];
		Lambda.col(i)	= psi_R1/e_R1 + psi_R2/(p.col(i)%(e_R2+c_q));
		d.col(i)		= e.col(i) + gamma(i)*p.col(i);
		V.col(i) 		= log(d.col(i)/gamma(i)) - X*beta - beta_int(i)+ log(Lambda.col(i) );
	}
	mat d_sgn	= d % sgn_e;
	mat dLambda	= d_sgn / Lambda;
	vec J1 		= (1 + lambda1 % sum(dLambda, 1)) % (1 + lambda2 % sum(dLambda/square(p),1)) - lambda1 % lambda2 % square(sum(dLambda/p,1));
	vec J		= log(J1) - sum(log(d) % sgn_e, 1);
	vec ll		= J - sum(exp(-V/sigma), 1) - sum(sgn_e % V/sigma + sgn_e*log(sigma), 1);
	
	List out = List::create(_["ll"]=ll, _["V"]=V, _["d"] = d,
							_["lambda1"]=lambda1, _["lambda2"]=lambda2, _["Lambda"]=Lambda, _["J1"]=J1);
	return(out);
}

// [[Rcpp::export]]
mat MDCEV_grad_fnC(vec param, int nx,
double c_q, arma::vec y, arma::vec Q, arma::mat e, arma::mat p, List X_list, 
vec ll, mat V, mat d, vec lambda1, vec lambda2, mat Lambda, vec J1){
	int nparam = param.n_elem, N = e.n_rows, R = e.n_cols;
	double psi_R1 = exp(param(nx+R)), psi_R2 = exp(param(nx+R+1));
	double sigma = exp(param(nx+R+2));
	// double sigma = 1.0;
	vec e_R1 = y - sum(e, 1);
	vec e_R2 = Q - sum(e/p,1);
	umat sgn_e = e>0;
	mat matone(N, R, fill::ones);
	mat sgn0 = matone - sgn_e;
	mat out = zeros(N, nparam);
	
	// Gradient w.r.t beta;
	for(int i=0; i<nx; i++){
		for(int r=0; r<R; r++){
			mat X = X_list[r];
			//out.col(i) += dpdV.col(r) % X.col(i);
		}
	}
	
	double eps = 0.001;
	for(int i=nx; i<nparam; i++){
	// for(int i=nx; i<(nparam-2); i++){
		vec param1 = param;
		param1(i) = param1(i) + eps;
		List ll_out1 =  MDCEV_ll_fnC(param1, nx, c_q, y, Q, e, p, X_list);
		vec ll1 = ll_out1["ll"];
		out.col(i) = (ll1 - ll)/eps;
	}	
	return(out);
}

// [[Rcpp::export]]
NumericVector MDCEV_LogLike_fnC(vec param, int nx,
double c_q, arma::vec y, arma::vec Q, arma::mat e, arma::mat p, List X_list){
	List ll_out = MDCEV_ll_fnC(param, nx, c_q, y, Q, e, p, X_list);
	vec llvec = ll_out["ll"], lambda1 = ll_out["lambda1"], lambda2 = ll_out["lambda2"], J1 = ll_out["J1"];
	mat V = ll_out["V"], d = ll_out["d"], Lambda = ll_out["Lambda"];
// 	mat grad = MDCEV_grad_fnC(param, nx, c_q, y, Q, e, p, X_list, llvec, V, d, lambda1, lambda2, Lambda, J1);
	NumericVector ll = as<NumericVector>(wrap(llvec));
// 	ll.attr("gradient") = grad;
	return(ll);
}

