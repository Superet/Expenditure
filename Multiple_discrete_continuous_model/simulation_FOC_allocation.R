# MDCEV simulation 

library(reshape2)
library(ggplot2)
library(evd)
library(nloptr)
set.seed(687)

# Set simulation parameters
R 	<- 3
N	<- 200
beta	<- c(-6, -8, -5)

# Omega is a list, each element is the first order derivative of expected utilty w.r.t. y
omega.ls	<- vector("list", length = R)
omega.ls[[1]]	<- function(x, par.ret = FALSE){ 
	par <- c(.1, -2, 10); 
	if(par.ret){
		return(par)
	}else{
		return(par[1]*x^2 + par[2]*x + par[3])
	}
}
omega.ls[[2]]	<- function(x, par.ret = FALSE){ 
	par <- c(.2, -1, 5); 
	if(par.ret){
		return(par)
	}else{
		return(par[1]*x^2 + par[2]*x + par[3])
	}
}
omega.ls[[3]]	<- function(x, par.ret = FALSE){ 
	par <- c(.3, -3, 8); 
	if(par.ret){
		return(par)
	}else{
		return(par[1]*x^2 + par[2]*x + par[3])
	}
}
tmp	<- seq(min(Inc), max(Inc), length = 100)
ggtmp	<- lapply(1:R, function(i) omega.ls[[i]](tmp))
ggtmp	<- data.frame(x = rep(tmp, R), omega = unlist(ggtmp), k = rep(1:R, each = length(tmp)))
ggplot(ggtmp, aes(x, omega, col = factor(k))) + geom_line()

# Simulate income
Inc <- rexp(N, .05)
summary(Inc)
eps_draw<- matrix(rgumbel(N*(R+1)), N, R+1)
psi		<- exp(rep(1, N) %*% t(c(beta, 0)) + eps_draw)

# Optimal allocation funciton 
optim_fn1	<- function(psi, omega.ls, Inc){
	# Reverse omega function: psi_i*omega.ls[i](x) = mu
	par.mat	<- sapply(omega.ls, function(x) x(0, TRUE))
	rev.omega <- vector("list", length = length(omega.ls))
	rev.omega[[1]]	<- function(mu, psi){
		par <- omega.ls[[1]](0, TRUE)
		return(-.5*par[2]/par[1] - .5*sqrt(par[2]^2 - 4*par[1]*(par[3] - mu/psi))/par[1])
	}
	rev.omega[[2]]	<- function(mu, psi){
		par <- omega.ls[[2]](0, TRUE)
		return(-.5*par[2]/par[1] - .5*sqrt(par[2]^2 - 4*par[1]*(par[3] - mu/psi))/par[1])
	}
	rev.omega[[3]]	<- function(mu, psi){
		par <- omega.ls[[2]](0, TRUE)
		return(-.5*par[2]/par[1] - .5*sqrt(par[2]^2 - 4*par[1]*(par[3] - mu/psi))/par[1])
	}
	
	# Function that solves for mu: 
	# sum_c y_c^*(mu) + y0(mu) = Inc
	solve_mu <- function( pind, psi, Inc){
		f <- function(mu){
			ystar	<- sapply(1:R, function(i) rev.omega[[i]](mu, psi[i]))
			out	<- sum(ystar * pind) + 1/mu - Inc
			return(out)
		}
		mu.star	<- uniroot(f, interval = c(min.mu, max.mu))
		return(mu.star)
	}
	
	bu 	<- sapply(1:R, function(i) psi[i]*omega.ls[[i]](0) )
	idx	<- order(bu, decreasing = T)
	sorted.bu	<- bu[idx]
	pind	<- 1 * (bu >= 1/Inc)
	M 		<- sum(pind > 0)
	min.mu	<- max(psi* (par.mat[3,] - .25*par.mat[2,]^2/par.mat[1,]))
	max.mu	<- max(psi* (par.mat[1,]*Inc^2 + par.mat[2,]*Inc^2 + par.mat[3,]))
	max.mu 	<- ifelse(max.mu <= min.mu, min.mu + 5000, max.mu)
	mu		<- solve_mu(pind, psi, Inc)$root
	while( M > 0 & sorted.bu[M] < mu ){
		print(M)
		pind[idx[M]]	<- 0
		M		<- M - 1
		if(M==0){
			mu		<- 1/Inc
		}else{
			mu		<- solve_mu(pind, psi, Inc)$root
		}
		print(mu)
	}
	ystar	<- sapply(1:R, function(i) rev.omega[[i]](mu, psi[i]))
	return(ystar)
}

psi <- c(.5, .6, .1)
