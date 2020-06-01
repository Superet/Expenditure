// Model versions;

#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double MDCEV_ll_fnC(vec param, int nx, mat shr, vec y, mat p, List X_list, uvec idx, int nh, vec M, int base=0){
	int R = shr.n_cols, N = shr.n_rows;
	vec beta = param.subvec(0, nx-1), gamma = exp(param.subvec(nx+R-1, nx+2*R-2));
	vec beta0 = zeros(R);
	for(int i=0; i<R; i++){
		if(i<base-1){
			beta0(i)	= param(nx+i); 
		}else if(i>base-1){
			beta0(i) 	= param(nx+i-1); 
		}
	}
	double sigma = exp(param(nx+2*R-1)); 
	
	// Organize random effect; 
	mat m1;
	m1.insert_cols(0, param.subvec(nx+2*R, param.n_elem-1));
	m1.reshape(R-1, nh);
	mat m2 = trans(m1);
	m2.insert_cols(base-1, zeros(nh,1)); 
	mat beta0_re	= m2.rows(idx-1);  
		
	// Utility V_tilde and indicator matrix for the original R alternatives;
	mat V_tilde = zeros(N, R), d_tilde = zeros(N, R);
	umat sgn_shr = shr > 0;
	for(int i=0; i<R; i++){
		mat Xtmp 		= X_list[i];
		d_tilde.col(i)	= shr.col(i) + gamma(i)*p.col(i)/y;
		V_tilde.col(i) 	= Xtmp * beta + beta0(i) + beta0_re.col(i) - log(d_tilde.col(i)/gamma(i));
	}
	
	// Jacobean matrix determinant; 
	mat d_sgn = d_tilde % sgn_shr;
	vec J = log(sum(d_sgn, 1)) - sum(log(d_tilde) % sgn_shr, 1);
	
	// Sum up the log likelihood; 
	vec ll = J - (M - ones<vec>(N)) * log(sigma) + sum(V_tilde % sgn_shr/sigma, 1) 
			- M % log(sum(exp(V_tilde/sigma), 1)) ; 
	
	double out = sum(ll); 
	return(-out);
}

// [[Rcpp::export]]
vec MDCEV_grad_fnC(vec param, int nx, mat shr, vec y, mat p, List X_list, uvec idx, int nh, vec M, int base=0){
	int R = shr.n_cols, N = shr.n_rows;
	vec beta = param.subvec(0, nx-1), gamma = exp(param.subvec(nx+R-1, nx+2*R-2));
	vec beta0 = zeros(R); 
	uvec subR(R-1); 
	for(int i=0; i<R; i++){
		if(i<base-1){
			beta0(i)	= param(nx+i); 
			subR(i)		= i; 
		}else if(i>base-1){
			beta0(i) 	= param(nx+i-1); 
			subR(i-1)	= i; 
		}
	}
	double sigma = exp(param(nx+2*R-1)); 
	
	// Organize random effect; 
	mat m1;
	m1.insert_cols(0, param.subvec(nx+2*R, param.n_elem-1));
	m1.reshape(R-1, nh);
	mat m2 = trans(m1);
	m2.insert_cols(base-1, zeros(nh,1)); 
	mat beta0_re	= m2.rows(idx-1);  
		
	// Utility V_tilde and indicator matrix for the original R alternatives;
	mat V_tilde = zeros(N, R), d_tilde = zeros(N, R), eV = zeros(N, R),
		sum_pX = zeros(N, nx), eVx = zeros(N, nx);
	umat sgn_shr = shr > 0;
	mat sgn_shr1 = arma::conv_to<arma::mat>::from(sgn_shr); 
	for(int i=0; i<R; i++){
		mat Xtmp 		= X_list[i];
		d_tilde.col(i)	= shr.col(i) + gamma(i)*p.col(i)/y;
		V_tilde.col(i) 	= Xtmp * beta + beta0(i) + beta0_re.col(i) - log(d_tilde.col(i)/gamma(i));
		eV.col(i) 	 	= exp(V_tilde.col(i)/sigma); 
		for(int j=0; j<nx; j++){
			sum_pX.col(j)	= sum_pX.col(j) + sgn_shr.col(i) % Xtmp.col(j); 
			eVx.col(j) 		= eVx.col(j) + eV.col(i) % Xtmp.col(j); 
		}
	}
	vec sum_eV = sum(eV, 1), sum_pd = sum(sgn_shr % d_tilde, 1); 
	
	vec grad = zeros(param.n_elem); 
	mat dld0 = zeros(N, R-1); 		// Derivative w.r.t beta0; 
	// Gradient w.r.t. beta; 
	for(int i=0; i<nx; i++){
		grad(i) = sum(sum_pX.col(i)/sigma - M/sigma % eVx.col(i) / sum_eV); 
	}
	// Gradient w.r.t. beta0; 
	for(int i=0; i<R; i++){
		if(i<base-1){
			dld0.col(i)	= sgn_shr1.col(i)/sigma  - M % eV.col(i) / sum_eV/sigma; 
			grad(i+nx)	= sum(dld0.col(i)); 
		}else if(i > base-1){
			dld0.col(i-1) = sgn_shr1.col(i)/sigma  - M % eV.col(i) / sum_eV/sigma; 
			grad(i+nx-1)= sum(dld0.col(i-1)); 
		}
	}	
	// Gradient w.r.t. gamma, note that we need to multiply by gamma to compute the gradient w.r.t ln(gamma); 
	for(int i=0; i<R; i++){
		grad(nx+R-1+i)	= sum( gamma(i) * sgn_shr.col(i) % ( p.col(i)/y % (1/sum_pd - 1/d_tilde.col(i)) 
							+ shr.col(i)/(sigma*gamma(i)*d_tilde.col(i)) ) 
							- M/(sigma*gamma(i)) % eV.col(i) % shr.col(i) / (sum_eV % d_tilde.col(i)) ) ; 
	}
	// Gradient w.r.t. sigma, note that we need to multiply by sigma to compute the gradient w.r.t ln(sigma); 
	grad(nx+2*R-1)	= sum((ones<vec>(N) - M) - sum(sgn_shr % V_tilde,1)/sigma 
						+ M/sigma % sum(eV % V_tilde, 1) / sum_eV); 
							
	// Gradient w.r.t. random effect; 
	int	pidx = nx + 2*R; 
	for(int i=0; i<nh; i++){
		uvec idx_ind = (idx == (i+1)); 
		for(int j=0; j<R-1; j++){
// 			grad(pidx) = sum(idx_ind % ( sgn_shr1.col(subR(j))/sigma - M/sigma%eV.col(subR(j))/sum_eV )); 
			grad(pidx) = sum(idx_ind % dld0.col(j)); 
			pidx++; 
		}
	}
	return(-grad);
}

// [[Rcpp::export]]
vec num_grad(vec param, int nx, mat shr, vec y, mat p, List X_list, uvec idx, int nh, vec M, int base=0, double delta = 0.0001){
	vec grad = zeros(param.n_elem); 
	double ll0 = MDCEV_ll_fnC(param, nx, shr, y, p, X_list, idx, nh, M, base); 
	for(int i=0; i<param.n_elem; i++){
		vec param1 = param; 
		param1(i) = param1(i) + delta; 
		double ll1 = MDCEV_ll_fnC(param1, nx, shr, y, p, X_list, idx, nh, M, base); 
		grad(i) = ll1 - ll0;
	}
	return(grad/delta);
}

// [[Rcpp::export]]
double MDCEV_ll_L2(vec param, int nx, mat shr, vec y, mat p, List X_list, uvec idx, int nh, vec M, int base=0, double L2C = 0){
	int R = shr.n_cols; 
	double ll0 = MDCEV_ll_fnC(param, nx, shr, y, p, X_list, idx, nh, M, base); 
	
	// Organize random effect;
	vec re =  param.subvec(nx+2*R, param.n_elem-1); 
	double ll_add = arma::conv_to<double>::from(sum(pow(re, 2)));
	double out = ll0 + L2C*ll_add ; 
	return(out); 
}

// [[Rcpp::export]]
vec MDCEV_grad_L2(vec param, int nx, mat shr, vec y, mat p, List X_list, uvec idx, int nh, vec M, int base=0, double L2C = 0){
	int R = shr.n_cols; 
	vec grad0 = MDCEV_grad_fnC(param, nx, shr, y, p, X_list, idx, nh, M, base); 
	vec re =  param.subvec(nx+2*R, param.n_elem-1); 
	int	pidx = nx + 2*R; 
	for(int i=0; i<re.n_elem; i++){
		grad0(pidx+i) = grad0(pidx+i) + 2*L2C*re(i); 
	}
	return(grad0); 
}
