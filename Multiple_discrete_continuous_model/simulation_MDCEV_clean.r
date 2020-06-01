# MDCEV simulation 

library(reshape2)
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(nloptr)

# Set paraemters 
R		<- 3 									# Number of alternatives
Ra		<- R									# Number of alternatives + number of outside options, here I don't have an outside good
beta0 	<- c(0, -1, -1)							# Intercepts for each alternative
beta	<- c(.5, -.7)							# Coefficients for X
gamma0 	<- c(1, 1, 1)							# Satiation parameters
sigma 	<- 1									# Scale parameters of Extreme type 1

#################
# Simulate data #
#################
set.seed(666666)
nx 		<- length(beta)
N 		<- 500									# Number of observations
gamma	<- rep(1, N) %*% t(gamma0)				# gamma could also have covariates, here I let them to be constant
X_arr 	<- array(rnorm(N*R*nx), c(R, N, nx))	# X dimention: alternative - observation - variable
X_list  <- lapply(1:R, function(i) X_arr[i,,])	
price 	<- matrix(runif(N*R, 2, 4), N, R)
Q		<- runif(N, 1, 20)
y		<- rowSums(price) * Q/R 

eps_draw<- matrix(rgumbel(N*R), N, R)			# Draw epsilon from extreme value distribution
xbeta	<- do.call(cbind, lapply(1:R, function(i) X_list[[i]] %*% beta + beta0[i]))
psi		<- exp(xbeta + eps_draw)

# Use efficient algorithm to solve the allocation
all_vec_fn	<- function(y, psi, gamma, price, R){
	N 	<- length(y)
	bu	<- psi/price
	
	# For each observation(row), sort psi/price
	idx	<- t(apply(bu, 1, order, decreasing = T))		
	sorted.bu	<- t(apply(bu, 1, sort, decreasing = T))
	
	mu	<- rep(0, N)				# Lagrangian multiplier for each observation
	M	<- rep(0, N)				# Number of positive alternatives for each observation
	for(r in 1:R){
		sel		<- mu/y < sapply(1:N, function(i) sorted.bu[i,(M[i]+1)])
		if(sum(sel) > 0){
			M[sel]	<- M[sel] + 1
			pe.ind	<- t(sapply(1:N, function(i) 1*(1:R %in% idx[i,(1:M[i])])))
			mu[sel]	<- rowSums(gamma[sel,]*psi[sel,]*pe.ind[sel,]) / (1 + rowSums(gamma[sel,]*price[sel,]*pe.ind[sel,])/y[sel])
		}
	}
	e	<- gamma*(psi*y/mu - price)*pe.ind
	return(e)
}
e_mat	<- all_vec_fn(y, psi, gamma, price, R)

##############
# Estimation # 
##############
# Cpp source code of likihood function 
cpp.src	<- '
#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List param_assign(vec param, int nx, int R, int base=0){
	// This function assigns parameter values; 
	vec beta0_norm = zeros(1), gamma0_norm = ones<vec>(1); 
	vec beta = param.subvec(0, nx-1), beta0_sub = param.subvec(nx, nx+R-2), gamma = param.subvec(nx+R-1, nx+2*R-2);
	vec beta0 = zeros(R);
	if(base == 0){
		beta0 	= join_cols(beta0_norm, beta0_sub);
	}else{
		int j=0;
		for(int i=0; i<R; i++){
			if(i==base - 1){
				continue;
			}else{
				beta0(i) = beta0_sub(j); 
				j++; 
			}
		}
	}
	 
// 	vec gamma = ones<vec>(R);
	double sigma = exp(param(nx+2*R-1)); 
	List out = List::create(_["beta"]=beta, _["beta0"]=beta0, _["gamma"]=gamma, _["sigma"] = sigma) ;
	return(out);
}

// [[Rcpp::export]]
List MDCEV_ll_fnC(vec param, int nx, mat e, mat p, List X_list, int base=0){
	int R = e.n_cols, N = e.n_rows;
	List param_out = param_assign(param, nx, R, base);
	vec beta = param_out["beta"], beta0 = param_out["beta0"], gamma = param_out["gamma"];
	double sigma = param_out["sigma"];
	
	// Utility V_tilde and indicator matrix for the original R alternatives;
	mat V_tilde = zeros(N, R), d_tilde = zeros(N, R);
	umat sgn_e_tilde = e > 0;
	for(int i=0; i<R; i++){
		mat Xtmp 		= X_list[i];
		d_tilde.col(i)	= e.col(i) + gamma(i)*p.col(i);
		V_tilde.col(i) 	= Xtmp * beta + beta0(i) - log(d_tilde.col(i)/gamma(i));
	}
	
	// Jacobean matrix determinant; 
	mat d_sgn = d_tilde % sgn_e_tilde;
	vec J = log(sum(d_sgn, 1)) - sum(log(d_tilde) % sgn_e_tilde, 1);
	mat sgn_e = arma::conv_to<arma::mat>::from(sgn_e_tilde); 
	vec M = sum(sgn_e, 1);
	
	// Sum up the log likelihood; 
	vec ll = J - (M - ones<vec>(N)) * log(sigma) + sum(V_tilde % sgn_e_tilde/sigma, 1) 
			- M % log(sum(exp(V_tilde/sigma), 1)); 
	
	List out = List::create(_["ll"]=ll, _["V_tilde"]=V_tilde);
	return(out);
}
'
# Compile Cpp code
sourceCpp(code = cpp.src)

# A wrapper function of likelihood 
MDCEV_wrapper <- function(param) {
	MDCEV_ll_fnC(param, nx, e=e_mat, p=price, X_list = X_list)$ll
}

# Estimation with multiple initial values (full set of parameters)
theta_init	<- c(beta, beta0[-1], gamma0, log(sigma))			# NOTE: set the first alternative as base, normalize beta_1 = 0 
names(theta_init)	<- c(paste("beta_",1:nx, sep=""), paste("beta0_",2:R,sep=""), paste("gamma",1:R, sep=""), "ln_sigma")
cat("Initial values are:\n"); print(theta_init); cat("\n")
system.time(tmp <- MDCEV_wrapper(theta_init))

# Estimate all parameters 
pct <- proc.time()
sol	<- maxLik(MDCEV_wrapper, start=theta_init, method="BFGS")
summary(sol)
use.time	<- proc.time() - pct

# Fix scale parameter
myfix 		<- length(theta_init)
sol1	<- maxLik(MDCEV_wrapper, start=theta_init, method="BFGS", fixed = myfix)
summary(sol1)
