#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List MDCP_ll_fnC(vec param, int nx,
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
	mat Lambda = zeros(N,R), d = zeros(N,R), V = zeros(N,R), f = zeros(N,R), F = zeros(N,R);
	vec lambda1 = zeros(N), lambda2 = zeros(N);
	
	lambda1 = psi_R1/square(e_R1);
	lambda2 = psi_R2/square(e_R2 + c_q);
	for(int i=0; i<R; i++){
		mat X = X_list[i];
		Lambda.col(i)	= psi_R1/e_R1 + psi_R2/(p.col(i)%(e_R2+c_q));
		d.col(i)		= e.col(i) + gamma(i)*p.col(i);
		V.col(i) 		= log(d.col(i)/gamma(i)) - X*beta - beta_int(i) + log(Lambda.col(i) );
		f.col(i)		= exp(-0.5 * square(V.col(i))/pow(sigma,2)) / (sqrt(2*PI)*sigma);
		NumericVector tmp = as<NumericVector>(wrap(vectorise(V.col(i))));
		NumericVector tmp1 = pnorm(tmp, 0.0, sigma);
		F.col(i)		= as<vec>(tmp1);
	}
	mat d_sgn	= d % sgn_e;
	mat dLambda	= d_sgn / Lambda;
	vec J1 		= (1 + lambda1 % sum(dLambda, 1)) % (1 + lambda2 % sum(dLambda/square(p),1)) - lambda1 % lambda2 % square(sum(dLambda/p,1));
	vec J		= log(J1) - sum(log(d) % sgn_e, 1);
	vec ll		= J + sum(log(f) % sgn_e, 1) + sum(log(F) % (matone - sgn_e), 1);
	
	List out = List::create(_["ll"]=ll, _["V"]=V, _["f"]=f, _["F"]=F, _["d"] = d,
							_["lambda1"]=lambda1, _["lambda2"]=lambda2, _["Lambda"]=Lambda, _["J1"]=J1);
	return(out);
}

// [[Rcpp::export]]
mat MDCP_grad_fnC(vec param, int nx,
double c_q, arma::vec y, arma::vec Q, arma::mat e, arma::mat p, List X_list, 
vec ll, mat V, mat f, mat F, mat d, vec lambda1, vec lambda2, mat Lambda, vec J1){
	int nparam = param.n_elem, N = e.n_rows, R = e.n_cols;
	vec beta = param.subvec(0,nx-1), beta_int = param.subvec(nx, nx+R-1), gamma = param.subvec(nx+R, nx+2*R - 1);
	double psi_R1 = exp(param(nx+2*R)), psi_R2 = exp(param(nx+2*R+1));
	double sigma = exp(param(nx+2*R+2));
	
	// double sigma = 1.0;
	vec e_R1 = y - sum(e, 1);
	vec e_R2 = Q - sum(e/p,1);
	umat sgn_e = e>0;
	mat matone(N, R, fill::ones);
	mat sgn0 = matone - sgn_e;
	mat dpdV = V % sgn_e /pow(sigma, 2) - f % sgn0 / F;
	mat out = zeros(N, nparam);
	
	// Gradient w.r.t beta;
	for(int i=0; i<nx; i++){
		for(int r=0; r<R; r++){
			mat X = X_list[r];
			out.col(i) += dpdV.col(r) % X.col(i);
		}
	}
	for(int i=nx; i<nx+R; i++){
		out.col(i) = - sum(dpdV, 1);
	}
	
	// Gradient w.r.t gamma;
	// mat d_sgn	= d % sgn_e;
	// mat dLambda	= d_sgn / Lambda;
	// for(int i=0; i<R; i++){
	// 	vec dgamma_p = (lambda2 % (1 + lambda1 % sum(dLambda, 1))/(Lambda.col(i) % p.col(i)) 
	// 				    + lambda1 % (1 + lambda2 % sum(dLambda/square(p),1)) % p.col(i) / Lambda.col(i)
	// 					- 2*lambda1 % lambda2 % sum(dLambda/p,1)/Lambda.col(i) )/J1 
	// 					- p.col(i) /d.col(i) + V.col(i) % e.col(i) / (gamma(i) * pow(sigma,2) * d.col(i));
	// 	vec dgamma_0 = - f.col(i) % e.col(i) /(gamma(i) * F.col(i) % d.col(i));
	// 	out.col(nx+i) = dgamma_p % sgn_e.col(i) + dgamma_0 % sgn0.col(i);
	// }
	
	double eps = 0.0001;
	for(int i=nx+R; i<nparam; i++){
		vec param1 = param;
		param1(i) = param1(i) + eps;
		List ll_out1 =  MDCP_ll_fnC(param1, nx, c_q, y, Q, e, p, X_list);
		vec ll1 = ll_out1["ll"];
		out.col(i) = (ll1 - ll)/eps;
	}
	
	// Gradient w.r.t beta_q(the utility parameter of quantity constraint);
	// rowvec Rones(R, fill::ones);
	// mat delta2 	= psi_R2/(p % ((e_R2+c_q)*Rones));
	// vec dJ1_daq	= psi_R1/e_R1 % lambda2 % (1 + lambda1 % sum(dLambda, 1)) % sum(dLambda/(Lambda % square(p)), 1) - 
	// 			  lambda1 % (1 + lambda2 % sum(dLambda/square(p),1)) % sum(dLambda % delta2, 1) - 
	// 			  lambda1 % lambda2 % sum(dLambda/p,1) % sum(dLambda % (1-2*delta2)/p, 1);
	// out.col(nparam-1) = dJ1_daq / J1 - sum(dpdV % delta2, 1);
	// out.col(nparam-2) = -out.col(nparam-1) - sum(dpdV, 1);
	
	return(out);
}

// [[Rcpp::export]]
NumericVector MDCP_LogLike_fnC(vec param, int nx,
double c_q, arma::vec y, arma::vec Q, arma::mat e, arma::mat p, List X_list){
	List ll_out = MDCP_ll_fnC(param, nx, c_q, y, Q, e, p, X_list);
	vec llvec = ll_out["ll"], lambda1 = ll_out["lambda1"], lambda2 = ll_out["lambda2"], J1 = ll_out["J1"];
	mat V = ll_out["V"], f = ll_out["f"], F = ll_out["F"], d = ll_out["d"], Lambda = ll_out["Lambda"];
// 	mat grad = MDCP_grad_fnC(param, nx, c_q, y, Q, e, p, X_list, llvec, V, f, F, d, lambda1, lambda2, Lambda, J1);
	NumericVector ll = as<NumericVector>(wrap(llvec));
// 	ll.attr("gradient") = grad;
	return(ll);
}

