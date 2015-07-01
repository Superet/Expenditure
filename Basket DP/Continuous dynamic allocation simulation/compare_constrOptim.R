# Problem: 
# min f(x) s.t. ui%*%x - ci>=0
# gi(x) = ci - ui %*% x <= 0

# Test
# True answer by constrOptim
# 1. nloptr from nloptr package. 
# 2. Penalty method with hard penalty constant if gi>=0
# 3. Penalty method with quadratic cost function, by updating penalty parameter. 
# 4. Barrier method implemented in Rcpp

library(nloptr)
library(Rcpp)
library(RcppGSL)
setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Continuous dynamic allocation simulation")
sourceCpp("constrOptimC_barrier.cpp")

################################################
# Constrained optimization with penalty method # 
################################################
pen_constrOptim <- function(f, g, theta, ui, ci, method, outer.iter=100, outer.tol=1e-05){
	# NOTE: consistent with constrOptim in R, ui%*% theta - ci >=0 
	# Original: min f(x) s.t. gi(x) <= 0
	# Penalty function: min_x f(x) + mu*sum_i (max(0, gi(x))^2)
	fun <- function(theta, mu){
		gi	<- ci - ui %*% theta						
		p 	<- pmax(0, gi)^2
		f(theta) + mu*sum(p)
	}
	grad <- function(theta, mu){
		# Assume that p(gi+) = 0 if gi<=0 
		gi	<- ci - ui%*%theta
		gi1 <- gi
		gi1[gi<=0] <- 0
		dp 	<- rep(NA, length(gi))
		dp 	<- -2*(t(ui) %*% gi1)
		g(theta) + mu*dp
	}
	mu	<- 5
	iter<- 1
	norm<- outer.tol + 1
	theta_old 	<- theta
	r_old		<- fun(theta_old, mu)
	obj_old		<- f(theta_old)
	while(mu<1e5 & iter <= outer.iter){
		iter		<- iter + 1
		a 			<- optim(theta_old, fun, grad, method=method, mu=mu)
		theta_new 	<- a$par
		r_new		<- a$value
		obj_new		<- f(theta_new)
		norm		<- obj_old - obj_new
		# if(obj_new > obj_old){
		# 	cat("Updated value is greater than the old value");
		# 	break
		# }
		mu			<- mu*5
		theta_old 	<- theta_new
		r_old		<- r_new
		obj_old		<- obj_new
	}
	return(list(par=theta_new, value=obj_new, iter=iter, mu=mu, norm=norm))
}

fQP1 <- function(b){
	cs <- bvec - t(Amat) %*% b
	if(all(cs<=0)){
		fQP(b)
	}else{
		return(1e5)
	}

}

########
# Test #
########
# Explicit constrained optimization 
#  ui %*% theta - ci >= 0
fQP <- function(b) {-sum(c(0,5,0)*b)+0.5*sum(b*b)}
fQP_g <- function(b) { -c(0,5,0) + b }
Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1), 3, 3)
bvec       <- c(-8, 2, 0)
inits 	   <- c(2,-1,-1)
ui 	<- t(Amat)
ci	<- bvec
Ans<- constrOptim(c(2,-1,-1), fQP, NULL, ui = t(Amat), ci = bvec)

# nloptr
eval_g_ineq	<- function(b){
	bvec - t(Amat) %*% b
}
eval_jac_g_ineq <- function(b){
	- Amat
}
nloptr(inits, eval_f = fQP, eval_g_ineq = eval_g_ineq, eval_jac_g_ineq = eval_jac_g_ineq,opts = list("algorithm"="NLOPT_LN_COBYLA"))

# Unconstrained optimization with penalty
system.time(print(optim(inits, fQP1)))
system.time(print(pen_constrOptim(fQP,fQP_g,inits, ui, ci, "BFGS")))

# Translate the constrOptim in R to C++
constrOptimC(inits, ui, ci)
