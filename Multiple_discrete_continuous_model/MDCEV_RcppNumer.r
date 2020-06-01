
# MDCEV simulation 

library(reshape2)
library(ggplot2)
library(Rcpp)
library(maxLik)
library(evd)
library(nloptr)
library(RcppEigen)
library(RcppNumerical) 

# Set paraemters 
R		<- 3 									# Number of alternatives
Ra		<- R									# Number of alternatives + number of outside options, here I don't have an outside good
beta	<- c(-1, -1, .5, -.7)							# Coefficients for X
gamma0 	<- c(1, 1, 1)							# Satiation parameters
sigma 	<- 1									# Scale parameters of Extreme type 1

#################
# Simulate data #
#################
set.seed(666666)
nx 		<- length(beta) - R + 1
N 		<- 500									# Number of observations
gamma	<- rep(1, N) %*% t(gamma0)				# gamma could also have covariates, here I let them to be constant
X_arr 	<- array(rnorm(N*R*nx), c(R, N, nx))	# X dimention: alternative - observation - variable
tmp		<- matrix(0, N, R-1)
X_list  <- lapply(1:R, function(i) {tmp1 <- tmp; if(i>1){tmp1[,i-1] <- 1}; cbind(tmp1, X_arr[i,,])} )	
price 	<- matrix(runif(N*R, 2, 4), N, R)
Q		<- runif(N, 1, 20)
y		<- rowSums(price) * Q/R 
Xlong 	<- do.call(rbind, X_list)

eps_draw<- matrix(rgumbel(N*R), N, R)			# Draw epsilon from extreme value distribution
xbeta	<- do.call(cbind, lapply(1:R, function(i) X_list[[i]] %*% beta ))
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

##################
# MLE likelihood # 
##################
src <- '
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>

using namespace Numer;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

class Mdcev: public MFuncGrad
{
private:
    const MapMat Xlong;
    const MapMat e;
    const MapMat e_sgn;
    const MapMat price; 
public:
    Mdcev(const MapMat x_, const MapMat e_, const MapMat esgn_, const MapMat price_) : Xlong(x_), e(e_), e_sgn(esgn_), price(price_){}

    double f_grad(Constvec& beta, Refvec grad)
    {	
    	const int N = e.rows(); 
    	const int nx = Xlong.cols(); 
    	const int R = Xlong.rows()/N; 

    	//const Eigen::VectorXd beta = param; 
    	//beta.setZero(); 
    	//beta.segment(1,nx-1) = param; 
    	const Eigen::VectorXd gammav = Eigen::VectorXd::Ones(R); 
    	//gammav.setOnes(); 
    	const double sigma = 1; 
    	
   		Eigen::VectorXd xbeta = Xlong * beta; 
   		//Eigen::MatrixXd V_tilde(N,R), d_tilde(N,R); 
   		Eigen::ArrayXXd V_tilde(N,R), eV(N, R), d_tilde(N,R), VX(N, nx), eVX(N, nx); 
   		Eigen::ArrayXd M(N); 
   		Eigen::ArrayXXd e_sgn1(N,R); 
   		M = e_sgn.rowwise().sum(); 
   		e_sgn1 = e_sgn; 
   		VX.setZero(); 
   		eVX.setZero(); 
   		
   		for(int i=0; i<R; i++){
   			Eigen::ArrayXd tmp = xbeta.segment(i*N,N);
   			d_tilde.col(i) 	= e.col(i) + price.col(i)*gammav(i); 
   			V_tilde.col(i) 	= tmp - (d_tilde.col(i)/gammav(i)).log();
   			eV.col(i) 		= (V_tilde.col(i)/sigma).exp(); 
   			VX 				= VX + Xlong.array().block(i*N,0,N,nx).colwise() * e_sgn1.col(i);
   			eVX 			= eVX + Xlong.array().block(i*N,0,N,nx).colwise() * eV.col(i); 
   		}
   		Eigen::ArrayXd eVsum = eV.rowwise().sum(); 
        Eigen::ArrayXd J= (d_tilde*e_sgn1).rowwise().sum().log() - 
    						((d_tilde.log()) * e_sgn1).rowwise().sum(); 
    						
    	Eigen::ArrayXd ll = J - (M-1)*log(sigma) + (V_tilde*e_sgn1/sigma).rowwise().sum() - 
    						M * eVsum.log(); 
 
    	const double f = -ll.sum(); 
    	
    	//Gradient;
    	// partial ll/partial beta_k = sum_{i in e_i>0} X_ik/sigma -M/sigma*sum_{i=1}^R e^(V_i/sigma)*X_ik/sum_{i=1}^R e^(V_i/sigma)
    	grad = -VX.colwise().sum() + (eVX.colwise()*(M/eVsum/sigma)).colwise().sum(); 
    	
        return f;
    }
};

// [[Rcpp::export]]
Rcpp::List solvMDCEV(Rcpp::NumericMatrix x, Rcpp::NumericMatrix e, Rcpp::NumericMatrix esgn0, Rcpp::NumericMatrix price, Rcpp::NumericVector inits)
{
    const MapMat xx = Rcpp::as<MapMat>(x);
    const MapMat ee = Rcpp::as<MapMat>(e);
    const MapMat pp = Rcpp::as<MapMat>(price);
    const MapMat esgn = Rcpp::as<MapMat>(esgn0); 
    // Negative log likelihood
    Mdcev nll(xx, ee, esgn, pp);
    // Initial guess
    Eigen::VectorXd beta(xx.cols());
    beta.setZero();
    
    double fopt;
    int status = optim_lbfgs(nll, beta, fopt);
    if(status < 0)
        Rcpp::stop("fail to converge");

  return Rcpp::List::create(
        Rcpp::Named("xopt") = beta,
        Rcpp::Named("fopt") = fopt,
        Rcpp::Named("status") = status
    );
// 	return Rcpp::wrap(beta);
}

// [[Rcpp::export]]
Rcpp::List evalMDCEV(Rcpp::NumericMatrix x, Rcpp::NumericMatrix e, Rcpp::NumericMatrix esgn0, Rcpp::NumericMatrix price, Rcpp::NumericVector inits)
{
    const MapMat xx = Rcpp::as<MapMat>(x);
    const MapMat ee = Rcpp::as<MapMat>(e);
    const MapMat pp = Rcpp::as<MapMat>(price);
    const MapMat esgn = Rcpp::as<MapMat>(esgn0); 
    // Negative log likelihood
    Mdcev nll(xx, ee, esgn, pp);
    // Initial guess
    Eigen::VectorXd beta(xx.cols()), grad(xx.cols() );
    for(int i=0; i<inits.size(); i++){
   	beta(i) = inits(i); 
    }
    Rcpp::Rcout << "beta =" << beta << std::endl;
    
	double res = nll.f_grad(beta, grad); 
	
  return Rcpp::List::create(
        Rcpp::Named("ll") = res,
        Rcpp::Named("grad") =grad
    );
}
'

sourceCpp(code = src)

sol <- solvMDCEV(Xlong, e_mat, esgn0  = 1*(e_mat>0), price = price, inits = beta - .1)
sol

out <- evalMDCEV(Xlong, e_mat, esgn0  = 1*(e_mat>0), price, inits = beta )
out

# MLE in R code 
eval_ll <- function(beta ){
	e_sgn 	<- 1 *(e_mat >0)
	M 		<- rowSums(e_sgn)
	xbeta <- Xlong %*% beta
	d_tilde <- e_mat + price/(rep(1,N) %*% matrix(gamma0, nrow = 1))
	V_tilde <- matrix(xbeta, ncol = R) - log(d_tilde/(rep(1,N) %*% matrix(gamma0, nrow = 1)))
	eV  <- exp(V_tilde/sigma)
	eVsum	<- rowSums(eV)
	tmp	<- lapply(1:R, function(i) e_sgn[,i]*Xlong[seq((i-1)*N+1, i*N),])
	tmp1 <- lapply(1:R, function(i) eV[,i]*Xlong[seq((i-1)*N+1, i*N),])
	VX <- tmp[[1]]; for(i in 2:length(tmp)){VX <- VX + tmp[[i]]}
	eVX <- tmp1[[1]]; for(i in 2:length(tmp)){eVX <- eVX + tmp1[[i]]}
	J 	<- log(rowSums(d_tilde * e_sgn)) - rowSums(log(d_tilde) * e_sgn)
	ll  <- J - (M -1) * log(sigma) + rowSums(V_tilde * e_sgn/sigma) - M *log(eVsum)
	ll 	<- sum(ll)
	grad <- VX/sigma - M*eVX/eVsum
	attr(ll, "gradient") 	<- colSums(grad)
	return(ll)
}

eval_ll(beta)
maxLik(eval_ll, start = rep(0, length(beta)), method = "BFGS")

# Compare the time 
library(microbenchmark)
microbenchmark(solvMDCEV(Xlong, e_mat, esgn0  = 1*(e_mat>0), price = price, inits = beta - .1), 
				maxLik(eval_ll, start = rep(0, length(beta)), method = "BFGS"))
				
				
