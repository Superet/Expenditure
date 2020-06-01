# MDCEV simulation 

library(reshape2)
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(nloptr)

setwd("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model")
source("0_Allocation_function.R")

# Set paraemters 
R		<- 3 		# Number of alternatives
Ra		<- R		# Number of alternatives + number of outside options
exp_outside <- quant_outside <- FALSE
beta0 	<- c(0, -1, -1)
beta	<- c(.5, -.7)
gamma0 	<- gamma	<- c(1, 1, 1)
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

eps_draw<- matrix(rgumbel(N*R), N, R)
xbeta	<- do.call(cbind, lapply(1:R, function(i) X_list[[i]] %*% beta + beta0[i]))
psi		<- exp(xbeta + eps_draw)

# par(mfrow=c(3,1))
# hist(y, breaks=100)
# hist(Q, breaks=100)
# hist(as.vector(price), breaks=100)

# Use optimization for solve optimal allocation
e_mat <- matrix(NA, N, Ra)
omega	<- rep(NA, N)
pct		<- proc.time()
for(i in 1:N){
	tmp	<- Allocation_fn(y = y[i], psi = psi[i,], gamma, Q = Q[i], price = price[i,], R, Ra, qz_cons, exp_outside, quant_outside)
	e_mat[i,] 	<- tmp$e
	omega[i] 	<- tmp$max
}
use.time1	<- proc.time() - pct

# Use efficient algorithm for solution
allc_fn <- function(y, psi, gamma, price, R){
	bu	<- psi/price
	idx	<- order(bu, decreasing = T)
	sorted.bu	<- bu[idx]
	mu	<- 0
	M 	<- 0
	e	<- rep(0, R)
	mu.track	<- NULL
	while(mu/y < sorted.bu[(M+1)] & M < R){
		M	<- M + 1
		sel	<- idx[1:M]
		mu	<- sum(gamma[sel]*psi[sel]) / (1 + sum(gamma[sel]*price[sel]/y))
		mu.track	<- c(mu.track, mu)
	}
	sel 	<- idx[1:M]
	e[sel] 	<- gamma[sel]*(psi[sel]*y/mu - price[sel])
	return(list (e = e, mu = mu, mu.track =mu.track))
}

all_vec_fn	<- function(y, psi, gamma, price, R){
	N 	<- length(y)
	bu	<- psi/price
	idx	<- t(apply(bu, 1, order, decreasing = T))
	sorted.bu	<- t(apply(bu, 1, sort, decreasing = T))
	mu	<- rep(0, N)
	M	<- rep(0, N)
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

e_mat1 <- matrix(NA, N, R)
omega1	<- rep(NA, N)
mu		<- rep(NA, N)
pct		<- proc.time()
for(i in 1:N){
	tmp			<- allc_fn(y = y[i], psi = psi[i,], gamma, price = price[i,], R)
	e_mat1[i,] 	<- tmp$e
	mu[i]		<- tmp$mu
	omega1[i]	<- uP_fn(e=e_mat1[i,], psi = psi[i,], gamma, price=price[i,], R, Ra, qz_cons, exp_outside, quant_outside)
}
use.time2 <- proc.time() - pct

summary(omega - omega1)

pct	<- proc.time()
e_mat2	<- all_vec_fn(y, psi, rep(1, N) %*% t(gamma), price, R)
use.time3 <- proc.time() - pct
max(abs(e_mat1 - e_mat2))

# Compare timing 
c(use.time1[3], use.time2[3], use.time3[3])

################################
# Add a new choice alternative #
################################
beta0.new	<- beta0[c(1:R,1)]
gamma.new	<- gamma[c(1:R,1)]
X.list.new	<- vector("list", length = R+1)
X.list.new	<- X_list[c(1:R, 1)]
eps_draw_new<- cbind(eps_draw, rgumbel(N))
xbeta.new	<- do.call(cbind, lapply(1:(R+1), function(i) X.list.new[[i]] %*% beta + beta0.new[i]))
psi.new		<- exp(xbeta.new + eps_draw_new)
price.new	<- price[,c(1:R,1)]

e_mat1.new 	<- matrix(NA, N, R+1)
omega1.new	<- rep(NA, N)
mu.new		<- rep(NA, N)
pct		<- proc.time()
for(i in 1:N){
	tmp				<- allc_fn(y = y[i], psi = psi.new[i,], gamma.new, price = price.new[i,], R+1)
	e_mat1.new[i,] 	<- tmp$e
	mu.new[i]		<- tmp$mu
	omega1.new[i]	<- uP_fn(e=e_mat1.new[i,], psi = psi.new[i,], gamma.new, price=price.new[i,], R+1, Ra+1, qz_cons, exp_outside, quant_outside)
}
use.time2 <- proc.time() - pct

summary(mu.new - mu)
summary(e_mat1)
summary(e_mat1.new)

cat("Average share before:\n"); print(colMeans(e_mat1/y))
cat("Average share after:\n"); print(colMeans(e_mat1.new/y))

###########################################################
# Sensitivity analysis of the correlation in random terms # 
###########################################################
xi		<- rgumbel(N)
ggtmp	<- data.frame()
rho.vec	<- seq(-1, 1, .1)

for(i in 1:length(rho.vec)){
	eps_draw.new<- cbind(eps_draw, rho.vec[i]*eps_draw[,1]+sqrt(1-rho.vec[i]^2)*xi)
	psi.new		<- exp(xbeta.new + eps_draw.new)
	tmp			<- all_vec_fn(y, psi.new, rep(1, N) %*% t(gamma.new), price.new, R+1)
	ggtmp		<- rbind(ggtmp, c(rho = rho.vec[i], colMeans(tmp/y)))
}
names(ggtmp)	<- c("rho", "a", "b", "c", "new.a")
tmp				<- setNames(c(rho = NA, colMeans(e_mat1/y), 0), c("rho", "a", "b", "c", "new.a") )
ggtmp			<- rbind(ggtmp, tmp)
ggtmp1			<- melt(ggtmp, id.vars = "rho")

ggplot(ggtmp, aes(rho, a + new.a)) + geom_line()
ggplot(ggtmp1, aes(rho, value)) + geom_line() + 
		geom_hline(data = subset(ggtmp1,is.na(rho)), aes(yintercept = value), col = "red", linetype = 2) + 
		facet_grid(. ~ variable) 

