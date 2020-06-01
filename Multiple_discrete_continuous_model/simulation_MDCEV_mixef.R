# Mixed effect MDCEV model 
# 
# Model: 
# U_ht = sum_r psi_{htr}*gamma_r*log(s_{htr}*y_{hr}/(p_{rt}*gamma_r) + 1)
# where 
# log(psi_{htr}) = beta_0{hr} + Xbeta + epsilon_{htr}
# epsilon_{htr} ~ GEV(0, sigma)


library(evd)
library(nloptr)
library(rstan)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(lbfgs)
library(ggplot2)
library(reshape2)
source('~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R')
sourceCpp(file = "~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/MDCEV_share_ME.cpp")

#################
# Simulate data # 
#################
# Set simulation parameters
set.seed(66)
nh	<- 20
T	<- 30
N	<- nh*T
R	<- 3
nx	<- 2
beta		<- c(.5, -1)
mu			<- c( -.7, .5)
tau			<- c( .2, .2)
sigma		<- .7
gamma		<- rep(1, R)
beta0_base	<- 1			# Set the one beta0 = 0
# mu[beta0_base]	<- 0
# tau[beta0_base]	<- 0

# Simulate data 
idx	<- rep(1:nh, each = T)
X	<- array(runif(R*N*nx, -1, 1), c(R, N, nx))
X_ls<- lapply(1:R, function(i) X[i,,])
p	<- matrix(rexp(N*R, 2), N, R)
y	<- rexp(N, .3)
summary(p)
summary(y)
beta0_re	<- NULL
for(i in 1:(R-1)){ beta0_re	<- cbind(beta0_re, rnorm(nh, mu[i], tau[i])) }
beta0_re.ext	<- kronecker(cbind(0, beta0_re), rep(1, T))
eps		<- matrix(rgev(N*R, scale = sigma), N, R)
psi		<- exp(sapply(1:R, function(i) beta0_re.ext[,i] + X[i,,]%*%beta + eps[,i]))
e  		<- sapply(1:N, function(i) Allocation_fn(y[i], psi[i,], gamma, Q = Inf, price = p[i,], R, Ra =R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE)$e)
e		<- t(e)
shr		<- e/y
shr_sgn	<- 1*(e > 0)
M		<- rowSums(shr_sgn)

# Combine parameters 
param	<- setNames(c(beta, mu, log(gamma), log(sigma), as.vector(t(beta0_re))), 
					c(paste("beta_", 1:nx, sep=''), paste("beta0_", setdiff(1:R, beta0_base), sep=""), 
						paste("lngamma", 1:R, sep=""), "ln_sigma", 
					 	paste("beta_re", rep(1:nh, each = R-1), "_", rep(setdiff(1:R, beta0_base), nh), sep="")) )
					
system.time(print(MDCEV_ll_fnC(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base)))

# Compare numerical gradient and analytical graident # 
g1 	<- MDCEV_grad_fnC(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base)
g2	<- num_grad(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base)
tmp	<- cbind(as.vector(g1), as.vector(g2))
max(abs(tmp[,2] - tmp[,1]))
microbenchmark(MDCEV_grad_fnC(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base), 
		num_grad(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base))

# -------------------------- #
# Estimation with L1 panelty #
ll_L1	<- function(param){
	MDCEV_ll_fnC(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base)
}

grad_L1	<- function(param){
	MDCEV_grad_fnC(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base)
}

est.l1 <- lbfgs::lbfgs(ll_L1, grad_L1, vars = rep(1, length(param)), 
						orthantwise_c=10, orthantwise_start = nx+2*R+1, orthantwise_end = length(param))
# compute se
par.est1	<- setNames(est.l1$par, names(param))
summary(abs(param - par.est1))
est.hes		<- optimHess(par.est1, fn = ll_L1, gr = grad_L1)
par.v1		<- chol2inv(chol(est.hes))
par.se1 	<- sqrt(diag(par.v1))

# -------------------------- #
# Estimation with L2 penalty #
system.time(print(MDCEV_ll_L2(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base, L2C = 10)))
g1b 	<- MDCEV_grad_L2(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base, L2C = 10)
L2C		<- 10
ll_L2	<- function(param, L2C = 0){
	MDCEV_ll_L2(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base, L2C = L2C)
}

grad_L2	<- function(param, L2C = 0){
	MDCEV_grad_L2(param, nx, shr, y, p, X_ls, idx, nh, M, beta0_base, L2C = L2C)
}

est.l2	<- lbfgs::lbfgs(ll_L2, grad_L2, vars = rep(1, length(param))  )
par.est2	<- setNames(est.l2$par, names(param))
summary(abs(param - par.est2))
est.hes2 	<- optimHess(par.est2, fn = ll_L2, gr = grad_L2)
par.se2		<- sqrt(diag(solve(est.hes2)))

# ------------------------------ #
# Compare the estimation results #
ggtmp	<- rbind(data.frame(Truev = param, est = par.est1, se = ifelse(is.null(par.se1), 0, par.se1), penalty = "L1"), 
				data.frame(Truev = param, est = par.est2, se = ifelse(is.null(par.se2), 0, par.se2), penalty = "L2"))
ggtmp$var	<- factor(rep(names(param), 2), levels = names(param))				
sel			<- grep("re", as.character(ggtmp$var))
ggplot(ggtmp[-sel,], aes(x=var, y = Truev)) + 
		geom_point() + 
		geom_pointrange(aes(y = est, ymin = est -1.96*se, ymax = est + 1.96*se, col = penalty), position = position_dodge(width = .5)) + 
		coord_flip()
quartz()
ggplot(ggtmp[sel,], aes(x=var, y = Truev)) + 
		geom_point() + 
		geom_pointrange(aes(y = est, ymin = est -1.96*se, ymax = est + 1.96*se, col = penalty), position = position_dodge(width = .5)) + 
		coord_flip()


# ------------------#
# Compare the speed #
microbenchmark(lbfgs::lbfgs(ll_L1, grad_L1, vars = rep(1, length(param)), invisible = 1,
						orthantwise_c=10, orthantwise_start = nx+2*R+1, orthantwise_end = length(param)), 
				lbfgs::lbfgs(ll_L2, grad_L2, vars = rep(1, length(param)), invisible = 1, L2C = L2C),  
				lbfgs::lbfgs(ll_L2, grad_L2, vars = rep(1, length(param)), invisible = 1)
				) 