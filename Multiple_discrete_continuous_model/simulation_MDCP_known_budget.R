# Simulation of allocation model 

# Model: 
# max_{e_s} \sum_{s=1}^{R} gamma_s*psi_s*log(e_s/p_s + 1) + psi_z*log(e_z) + psi_w*log(e_w)
# s.t. 	sum_{s=1}^{R} e_s + e_z <= y
# 	 	sum_{s=1}^{R} e_s/p_s + e_w <= Q

library(ggplot2)	
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)

setwd("~/Documents/Research/Store switching/Exercise")
sourceCpp("MDCP_function.cpp")

#####################
# Utility functions #
#####################
uP_fn <- function(e, psi, gamma,price, R){
# Negative purchase utility 
# e: R+2 vecotr
# psi, gammais (R+2)-dimension vector 
# price is R-dimension vector
	e0 		<- e[1:R]
	e_R1	<- e[(R+1):length(e)]
	u0		<- c(log(e0/(gamma[1:R]*price) + 1), log(e_R1[1]), log(e_R1[2]+qz_cons))
	out 	<- sum(psi*gamma*u0 ) 
	return(-out)
}

uPGrad_fn <- function(e, psi, gamma,price, R){
	e0 		<- e[1:R]
	e_R1	<- e[(R+1):length(e)]
	u0		<- c(1/(e0+price*gamma[1:R]), 1/e_R1[1], 1/(e_R1[2]+qz_cons) )
	out 	<- psi*gamma*u0 
	return(-out)
	
}

# Allocation function - solve a constrained optimization
Allocation_fn <- function(y, psi, gamma, Q=NULL, price,R,Ra, inits=NULL){
	# We have Ra non-negativity constraints, one budget constraint, and one quantity constraint
	# ui <- rbind(diag(R+2), c(rep(-1, R+1), 0) )
	# inits is a list of initial values
	ui <- rbind(diag(Ra), c(rep(-1, R+1),rep(0, Ra-R-1 )) )
	ci <- c(rep(0,R+1), -qz_cons, -y)
	if(is.null(inits)){
		tmpq  <- min(y/price, Q) * .9
		sel   <- which.min(price)
		tmp	  <- tmpq * min(price)
		inits <- list(c(rep(tmp/(Ra-1), Ra-1), Q - tmpq - 1e-6), 
					  c(rep(1e-8, R), y-(R+2)*1e-8, Q - sum(price)*1e-8 - 1e-6),
					  c(rep(tmp/R, R), y- tmp -1e-6, 1e-8) )
	}
	
	if(is.null(Q)){
		sol <- constrOptim(theta = inits, f=uP_fn, grad=NULL, ui=ui, ci=ci,psi=psi, gamma=gamma, price=price, R=R)
		return(list(e = sol$par, max = -sol$value))
	}else{
		if(Q<=0){
			e <- c(rep(0,R),y)
			if(Ra>R+1){
				e <- c(e, 0)
			}
			return(list(e = e, max = log(y) ) )
		}else{
			ui <- rbind(ui, c(-1/price, 0, -1))
			ci <- c(ci, -Q)
			sol.list <- vector("list", length(inits))
			for(j in 1:length(inits)){
				sol.list[[j]] <- try(constrOptim(theta = inits[[j]], f=uP_fn, grad=uPGrad_fn, ui=ui, ci=ci,psi=psi, 
									gamma=gamma, price=price, R=R) )
			}
			sol.list1 <- sol.list[sapply(sol.list, function(x) !inherits(x, "try-error"))]
			sel <- which.min(sapply(sol.list1, function(x) x$value))
			if(length(sel)==0){
				sol <- list(par = rep(NA, Ra), value=NA, convergence = NA)
			}else{
				sol <- sol.list1[[sel]]
				if(sol$convergence != 0){
					cat("Constrained optimization does not converge at value y=",y,", Q=",Q,"\n")
				}
			}
			return(list(e = sol$par, max = -sol$value, convergence = sol$convergence))
		}
	}
}

############################
# Simulate allocation data # 
############################
# Set parameters 
R 		<- 3
Ra		<- R+2
beta	<- c(1, -1)
beta_int<- c(0,-.5,.5)			# The alternative-specific intercepts
beta_o	<- c(.5, -1)			# The utility parameters of the two outside alternatives
gamma	<- c(1,1,.5)		# Ratiation paramters for only for R alternatives
sigma 	<- 1.5				# Standard deviation of the errors in psi
nx 		<- length(beta)	
qz_cons <- 0

N 		<- 1000
X_arr 	<- array(rnorm(N*R*nx), c(R, N, nx))
X_list  <- lapply(1:R, function(i) X_arr[i,,])
price 	<- matrix(runif(N*R), N, R)
Q		<- rexp(N, .8)
y		<- rowSums(price) * Q/3 + rexp(N)
eps		<- matrix(rnorm(N*R,0, sigma), N, R)
xbeta 	<- apply(X_arr,1,function(xx) xx%*% beta)

e_mat	<- matrix(NA, N, Ra)
for(i in 1:N){
	tmp.psi 	<- c(exp(xbeta[i,] + beta_int + rnorm(R, 0, sigma)), exp(beta_o))
	tmp.gamma	<- c(gamma, 1, 1)
	tmp			<- Allocation_fn(y[i], tmp.psi, tmp.gamma, Q[i], price[i,], R, Ra)
	e_mat[i,]	<- tmp$e
	if(i%%50==0){print(i)}
}
eps <- 1e-4
eR 	<- e_mat[,1:R]
eR[eR<=eps] <- 0
mean(eR==0)

# Esimation 
param <- c(beta, beta_int, gamma, beta_o, sigma)
names(param) <- c(paste("beta",1:nx,sep=""), paste("intercept",1:R,sep=""), paste("gamma",1:R, sep=""), "ln_psi1", "ln_psi2", "sigma")
param_init <- param 
param_init[1] <- param[1] + 1
param_init[2] <- param[2] + 3

system.time(tmp <- MDCP_LogLike_fnC(param_init, nx, qz_cons, y, Q, eR, price, X_list) )
param2 <- param_init
param2[6:7] <- 0

MCDP_wrapper <- function(param){
	return(MDCP_LogLike_fnC(param, nx, qz_cons, y, Q, eR, price, X_list))
}

param2 <- runif(length(param))
summary(sol.free <- maxLik(MCDP_wrapper, start=param2, method="BFGS"))

summary(sol.fixed <- maxLik(MCDP_wrapper, start=param_init, method="BFGS", fixed = c(F,F,F,F,F,T,T,T,F,F,F)))


# Check likihood function along the utility parameters for outside goods. 
my_seq <- seq(-2, 3, by=.2)
tmpname <- c("beta_o1", "beta_o2")
ggtmp <- data.frame(NULL)
for(i in 6:7){
	param1 <- param 
	tmpll <- rep(NA, length(my_seq))
	for(j in 1:length(my_seq)){
		param1[i] <- param[i] + my_seq[j]
		tmpll[j] <- sum(MCDP_wrapper(param1))
	}
	ggtmp <- rbind(ggtmp, data.frame(Variable=tmpname[(i-5)], par = param[i] + my_seq, ll = tmpll))
}

ggplot(ggtmp, aes(par, ll)) + geom_point() + geom_line() + 
		facet_grid(.~ Variable)

