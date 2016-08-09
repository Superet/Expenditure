# MDCEV simulation 

library(reshape2)
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(nloptr)

setwd("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model")
model_name <- "MDCEV_share"

# sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Allocation_function.R")

########################
# Estimation functions #
########################
src <- 
'
#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List MDCEV_ll_fnC(vec param, int nx, mat shr, vec y, vec s1_index, mat p, List X_list, int base=0){
	int R = shr.n_cols, N = shr.n_rows;
	vec beta0_norm = zeros(1), gamma0_norm = ones<vec>(1); 
	vec beta = param.subvec(0, nx-1), beta0_sub = param.subvec(nx, nx+R-2), 
		gamma1	= param.subvec(nx+R-1, 2*nx+R-2), gamma0 = param.subvec(2*nx+R-1, 2*(nx+R)-2);
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
	double sigma = exp(param(2*(nx+R)-1));

	// Utility V_tilde and indicator matrix for the original R alternatives;
	mat V_tilde = zeros(N, R), d_tilde = zeros(N, R), gamma = zeros(N, R);
	umat sgn_shr_tilde = shr > 0;
	for(int i=0; i<R; i++){
		mat Xtmp 		= X_list[i];
		gamma.col(i)	= exp(Xtmp * gamma1 + gamma0(i)); 
		d_tilde.col(i)	= shr.col(i) + gamma.col(i)%p.col(i)/y;
		V_tilde.col(i) 	= Xtmp * beta + beta0(i) - log(d_tilde.col(i)/gamma.col(i));
	}
	
	// Jacobean matrix determinant; 
	mat d_sgn = d_tilde % sgn_shr_tilde;
	vec J = log(sum(d_sgn, 1)) - sum(log(d_tilde) % sgn_shr_tilde, 1);
	mat sgn_shr = arma::conv_to<arma::mat>::from(sgn_shr_tilde); 
	vec M = sum(sgn_shr, 1);
	
	// Sum up the log likelihood; 
	vec ll = J - (M - ones<vec>(N)) * log(sigma) + sum(V_tilde % sgn_shr_tilde/sigma, 1) 
			- M % log(sum(exp(V_tilde/sigma), 1)) ; 
	
	List out = List::create(_["ll"]=ll, _["V_tilde"]=V_tilde);
	return(out);
}


// [[Rcpp::export]]
colvec MDCEV_LogLike_fnC(vec param, int nx, vec s1_index, arma::mat shr, vec y, arma::mat p, List X_list, int base=0){
	List ll_out	= MDCEV_ll_fnC(param, nx, shr, y, s1_index, p, X_list, base);
	colvec llvec 	= ll_out["ll"];
	
// 	NumericVector ll = as<NumericVector>(wrap(llvec));
// 	ll.attr("gradient") = grad;
	return(llvec);
}
'
sourceCpp(code = src)

#################
# Simualte data #
#################
# Set paraemters 
R		<- 3 		# Number of alternatives
Ra		<- R		# Number of alternatives + number of outside options
exp_outside <- quant_outside <- FALSE
beta0 	<- c(0, -1, -1)
beta	<- c(.5, -.7)
gamma1	<- c(-1, .6)
gamma0 	<- c(-1, 0, -1)
sigma 	<- 1
qz_cons	<- Inf

# Simulate data 
set.seed(666666)
nx 		<- length(beta)
N 		<- 500		# Number of observations
X_arr 	<- array(rnorm(N*R*nx), c(R, N, nx))
X_list  <- lapply(1:R, function(i) X_arr[i,,])
price 	<- matrix(runif(N*R, 2, 4), N, R)
Q		<- runif(N, 1, 20)
y		<- rowSums(price) * Q/R 

par(mfrow=c(3,1))
hist(y, breaks=100)
hist(Q, breaks=100)
hist(as.vector(price), breaks=100)

eps_draw<- matrix(rgumbel(N*R), N, R)
xbeta	<- do.call(cbind, lapply(1:R, function(i) X_list[[i]] %*% beta + beta0[i]))
psi		<- exp(xbeta + eps_draw)
gamma	<- do.call(cbind, lapply(1:R, function(i) X_list[[i]]%*%gamma1 + gamma0[i]))
gamma	<- exp(gamma)
e_mat <- matrix(NA, N, Ra)
for(i in 1:N){
	tmp	<- Allocation_fn(y = y[i], psi = psi[i,], gamma[i,], Q = Inf, price = price[i,], R, Ra, qz_cons, exp_outside, quant_outside)
	e_mat[i,] <- tmp$e
}
eR 	<- e_mat[,1:R]
mean(eR==0)

sel	<- apply(eR, 1, function(x) !all(x==0))
eR 	<- eR[sel,]
Q	<- Q[sel]
y	<- y[sel]
price <- price[sel,]
X_list <- lapply(X_list, function(x) x[sel,])
e1_index <- apply(eR, 1, function(x) which(x== min(x[x>0])))
shr		<- eR/y

# Esitmation 
MDCEV_wrapper <- function(param){
	MDCEV_LogLike_fnC(param, nx, s1_index = e1_index, shr = shr, y = y, p=price, X_list = X_list)
}

# Estimation with multiple initial values (full set of parameters)
theta_init0	<- c(beta, beta0[-1], gamma1, gamma0, log(sigma))
theta_init2	<- theta_init1 <- theta_init0
theta_init1[1:(nx+R)] <- theta_init0[1:(nx+R)] + rnorm(nx+R)
theta_init2 <- rep(1, length(theta_init2))
theta_init 	<- list(theta_init0, theta_init1, theta_init2)
system.time(tmp <- MDCEV_wrapper(theta_init[[1]]))

tmp_sol <- vector("list", length(theta_init))
pct <- proc.time()
for(i in 1:length(theta_init)){
	names(theta_init[[i]]) <- c(paste("beta_",1:nx, sep=""), paste("beta0_",2:R,sep=""), 
								paste("gamma1_",1:nx, sep=""), paste("gamma0_",1:R, sep=""),
								"ln_sigma")	
	tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS")
}
sel 	<- which.max(sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum)))
sol		<- tmp_sol[[sel]]
lapply(tmp_sol, summary)
cbind(True = theta_init0, sapply(tmp_sol, coef))

