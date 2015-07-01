// Model versions;
// a: outside good presence (1.No vs. 2.Yes); 
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
	vec beta0_norm = zeros(1); 
	vec beta = param.subvec(0, nx-1), beta0_sub = param.subvec(nx, nx+R-2);
	vec beta0 = join_cols(beta0_norm, beta0_sub);
// 	vec gamma = param.subvec(nx+R-1, nx+2*R-2);
// 	double sigma  = exp(param(nx+2*R-1));
	vec gamma = ones<vec>(R);
	double sigma = 1.0; 
	List out = List::create(_["beta"]=beta, _["beta0"]=beta0, _["gamma"]=gamma, _["sigma"] = sigma) ;
	return(out);
}

// [[Rcpp::export]]
List eAdjust(vec param, int nx, arma::vec e1_index, arma::mat e, arma::mat p, List X_list){
	int R = e.n_cols, N = e.n_rows;
	List param_out = param_assign(param, nx, R);
	vec beta = param_out["beta"], beta0 = param_out["beta0"], gamma = param_out["gamma"];
	mat e_bar= zeros(N, R-1), gamma_bar = zeros(N, R-1), p_bar = zeros(N, R-1), xbeta = zeros(N, R-1);
	vec e1 = zeros(N), gamma1 = zeros(N), p1 = zeros(N), xbeta1 = zeros(N);
	cube Xvar = zeros(N, nx, R);
	
	// Declare some working variables; 
	IntegerVector seqR_tmp = seq_len(R);
	vec seqR = as<vec>(seqR_tmp);
	mat eyeR = eye(R, R);
	uvec jvec(1);
	rowvec one_nx = ones<rowvec>(nx);
		
	// Loop over each alternatives; 
	for(int i=1; i<=R; i++){
		uvec selr = find(e1_index==i );
		uvec selc = find(seqR != i);	
		if(selr.n_elem > 0){
			vec one_n = zeros(N);
			one_n.elem(selr).ones();
		
			// Matrix associated with the alternatives other than the reference alternative; 
			mat T 				= eyeR.cols(selc);				
			vec onev			= ones<vec>(selr.n_elem);
			e_bar.rows(selr)	= e.rows(selr) * T;
			gamma_bar.rows(selr)= onev * trans(gamma(selc));
			p_bar.rows(selr)	= p.rows(selr) * T;			
			for(int j=0; j<R-1; j++){
				mat Xtmp = X_list[selc(j)];
				jvec.fill(j);
				xbeta(selr, jvec) = Xtmp.rows(selr) * beta + beta0(selc(j));
				Xvar.slice(j+1)	+= Xtmp % (one_n * one_nx);
			}
			
			// Vector for the reference alternative;
			vec e1v				= e.col(i-1);
			vec p1v				= p.col(i-1);
			e1(selr)			= e1v(selr);
			gamma1(selr)		= onev * gamma(i-1);
			p1(selr)			= p1v(selr);
			mat X1				= X_list[i-1];
			xbeta1(selr)		= X1.rows(selr) * beta + beta0(i-1);
			Xvar.slice(0)		+= X1 % (one_n * one_nx);
		}
	}
	cube Xvar1 = reshape(Xvar, N, nx*R, 1);
	vec Xvar_vec = vectorise(Xvar1.slice(0));
	List out = List::create(_["e1"]=e1,_["gamma1"]=gamma1, _["p1"]=p1, _["xbeta1"]=xbeta1,
							_["e_bar"]=e_bar, _["gamma_bar"]=gamma_bar, _["p_bar"]=p_bar, _["xbeta"]=xbeta,
							_["Xvar"] = Xvar_vec);
	return(out);
}

// [[Rcpp::export]]
List MDCEV_ll_fnC(vec param, int nx, mat e, List Adlist){
	int R = e.n_cols, N = e.n_rows;
	List param_out = param_assign(param, nx, R);
	double sigma = param_out["sigma"];
	vec e1 = Adlist["e1"], gamma1 = Adlist["gamma1"], p1 = Adlist["p1"], xbeta1 = Adlist["xbeta1"];
	mat e_bar = Adlist["e_bar"], gamma_bar = Adlist["gamma_bar"], p_bar = Adlist["p_bar"], xbeta = Adlist["xbeta"];
	vec Xvar_vec = Adlist["Xvar"];
	
	rowvec onev	= ones<rowvec>(R-1);
	mat matone(N, R-1, fill::ones);	
	umat sgn_e	= e_bar>0;
	mat sgn_0	= matone - sgn_e;
	
	mat d		= e_bar + gamma_bar % p_bar;
	mat d1		= e1 + gamma1 % p1;
	mat V		= log(d/gamma_bar) - xbeta - (log(d1/gamma1) -xbeta1 ) * onev;
	mat V_sgn 	= V % sgn_e;
	mat d_sgn	= d % sgn_e;
	
	// Compute the Jacobean determinant and log likelihood; 
	vec J		= log(sum(d_sgn,1) + d1) - sum( log(d)%sgn_e, 1)  - log(d1);
	vec ll		= J - sum(exp(-V/sigma) , 1) - sum(V_sgn/sigma + log(sigma)*sgn_e, 1);
	
	List out = List::create(_["ll"]=ll, _["V"]=V, _["d"]=d, _["d1"]=d1, _["sgn_e"]=sgn_e, _["J"]=J);
	return(out);
}

// [[Rcpp::export]]
mat MDCEV_grad_fnC(vec param, int nx, mat e, vec ll, mat V, List Adlist){
	int nparam = param.n_elem, N = e.n_rows, R = e.n_cols;
	List param_out = param_assign(param, nx, R);
	double sigma = param_out["sigma"];
	vec e1 = Adlist["e1"], gamma1 = Adlist["gamma1"], p1 = Adlist["p1"], xbeta1 = Adlist["xbeta1"];
	mat e_bar = Adlist["e_bar"], gamma_bar = Adlist["gamma_bar"], p_bar = Adlist["p_bar"], xbeta = Adlist["xbeta"];
	vec Xvar_vec = Adlist["Xvar"];
	cube Xvar(Xvar_vec.begin(), N, nx, R);
		
	mat matone(N, R-1, fill::ones);	
	umat sgn_e	= e_bar>0;
	mat sgn_0	= matone - sgn_e;
	
	mat out = zeros(N, nparam);
	cube X_dif = zeros(N, R-1, nx);
	for(int i=0; i<R-1; i++){
		X_dif(span(), span(i), span())  = Xvar.slice(0) - Xvar.slice(i+1); 
	}	
	
	for(int i=0; i<nx; i++){
		out.col(i) = (sum(exp(-V/sigma) % X_dif.slice(i), 1) - sum(X_dif.slice(i) % sgn_e, 1) )/sigma;
	}
	double eps = 0.01;
	for(int i=nx; i<nparam; i++){
		vec param1 = param;
		param1(i) = param1(i) + eps;
		List ll_out1 =  MDCEV_ll_fnC(param1, nx, e, Adlist);
		vec ll1 = ll_out1["ll"];
		out.col(i) = (ll1 - ll)/eps;
	}
	return(out);
}

// [[Rcpp::export]]
NumericVector MDCEV_LogLike_fnC(vec param, int nx, vec e1_index, arma::mat e, arma::mat p, List X_list){
	int R = e.n_cols, N = e.n_rows;
	List Adlist	= eAdjust(param, nx, e1_index, e, p, X_list);
	List ll_out	= MDCEV_ll_fnC(param, nx, e, Adlist);
	vec llvec 	= ll_out["ll"];
	mat V		= ll_out["V"];
// 	mat grad 	= MDCEV_grad_fnC( param, nx, e, llvec, V, Adlist);

	NumericVector ll = as<NumericVector>(wrap(llvec));
// 	ll.attr("gradient") = grad;
	return(ll);
}