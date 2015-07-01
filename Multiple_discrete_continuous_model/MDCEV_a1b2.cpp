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
List param_assign(vec param, int nx, int R){
	// Assign parameter values; 
	vec beta = param.subvec(0,nx-1), beta0 = param.subvec(nx, nx+R-1), gamma = param.subvec(nx+R, nx+2*R - 1);
	double psi_R1 = exp(param(nx+2*R)), psi_R2 = exp(param(nx+2*R+1));
	double sigma = exp(param(nx+2*R+2));

// 	vec beta0_norm = zeros(1); 
// 	vec beta = param.subvec(0, nx-1), beta0_sub = param.subvec(nx, nx+R-2);
// 	vec beta0 = join_cols(beta0_norm, beta0_sub);
// 	vec gamma = ones<vec>(R);
// 	double sigma = 1.0; 
// 	double psi_R2 = 0.01;
	List out = List::create(_["beta"]=beta, _["beta0"]=beta0, _["gamma"]=gamma, _["sigma"] = sigma, _["psi_R2"]=psi_R2) ;
	return(out);
}

// [[Rcpp::export]]
List MDCEV_ll_fnC(vec param, int nx, double c_q, vec Q, mat e, vec e1_index, mat p, List X_list){
	int R = e.n_cols, N = e.n_rows;
	List param_out = param_assign(param, nx, R);
	vec beta = param_out["beta"], beta0 = param_out["beta0"], gamma = param_out["gamma"];
	double psi_R2 = param_out["psi_R2"], sigma = param_out["sigma"];
	
	// Utility V_tilde according to the original order of alternative index; 
	umat sgn_e_tilde = e > 0;
	vec e_R2 = Q - sum(e/p,1);
	mat d_tilde = zeros(N,R), xbeta = zeros(N,R), V_tilde = zeros(N, R);
	d_tilde = e + (ones<vec>(N)*trans(gamma)) % p;
	for(int i=0; i<R; i++){
		mat Xtmp 	= X_list[i];
		xbeta.col(i)= Xtmp * beta + beta0(i);
		V_tilde.col(i) 	= xbeta.col(i) - log(d_tilde.col(i)/gamma(i));
	}
	
	// Adjust for the baseline alternative; 
	mat p_star = zeros(N, R), Lambda = zeros(N, R), V = zeros(N, R);
	vec d1 = zeros(N), gamma1 = zeros(N), xbeta1 = zeros(N), p1 = zeros(N);
	rowvec one_row = ones<rowvec>(R);
	uvec ivec(1);
	IntegerVector seqR_tmp = seq_len(R);
	vec seqR = as<vec>(seqR_tmp);

	for(int i=1; i<=R; i++){
		uvec selr 	= find(e1_index==i);
		uvec selc 	= find(seqR != i);
		if(selr.n_elem > 0){
			ivec.fill(i-1);
			vec oner 	= ones<vec>(selr.n_elem);
			d1(selr)	= d_tilde(selr,ivec);
			gamma1(selr)= gamma(i-1) * oner;
			xbeta1(selr)= xbeta(selr,ivec);
			p1(selr)	= p(selr,ivec);		
		
			p_star.rows(selr)	= 1/p.rows(selr) - 1/(p1(selr) * one_row);
			Lambda.rows(selr)	= gamma(i-1) * (exp(xbeta1(selr))/d1(selr)) * one_row + 
								psi_R2 * p_star.rows(selr)/((e_R2.rows(selr) + c_q)*one_row);
			V(selr, ivec).fill(0); 
			V(selr, selc)		= log(Lambda(selr, selc)) - V_tilde(selr, selc) ;
		}
	}
	
	/// Jacobean matrix determinant; 
	vec lambda2 		= psi_R2/square(e_R2 + c_q);
	mat d_sgn 			= d_tilde % sgn_e_tilde;
	mat d_sgn_Lambda	= d_sgn/Lambda;
	vec J1				= xbeta1 + log(gamma1) - log(d1) + log(lambda2) + 
		log(sum(d_sgn_Lambda,1) % (1/lambda2 + sum(d_sgn_Lambda % square(p_star), 1)) - square(sum(d_sgn_Lambda % p_star, 1)));
	vec J				= J1 - sum(log(d_tilde) % sgn_e_tilde, 1) ; 
	
	// Sum up the log likelihood; 
	vec M = sum(arma::conv_to<arma::mat>::from(sgn_e_tilde), 1);
	vec ll = J - (M - ones<vec>(N))*log(sigma) + sum(-V % sgn_e_tilde/sigma, 1) 
			- M % log(sum(exp(-V/sigma), 1)); 

	List out = List::create(_["ll"]=ll, _["V"]=V, _["p_star"]=p_star, _["Lambda"]=Lambda);
	return(out); 
}

// [[Rcpp::export]]
NumericVector MDCEV_LogLike_fnC(vec param, int nx, double c_q, vec Q, mat e, vec e1_index, mat p, List X_list){
	List ll_out = MDCEV_ll_fnC(param, nx, c_q, Q, e, e1_index, p, X_list);
	vec llvec = ll_out["ll"];
	if(! llvec.is_finite()){
		uvec sel = find_nonfinite(llvec);
		llvec(sel).fill(-100.0);
	}
// 	mat grad = MDCEV_grad_fnC(param, nx, c_q, Q, e, e1_index, p, X_list, llvec);
	NumericVector ll = as<NumericVector>(wrap(llvec));
// 	ll.attr("gradient") = grad;
	return(ll);
}
