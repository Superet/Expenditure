# MDCEV simulation 

library(reshape2)
library(ggplot2)
library(maxLik)
library(evd)
library(nloptr)
library(gtools)

# Functions 
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

S		<- 2
theta	<- c(.6, 1)
nest	<- matrix(c(1,1,0,0,0,1), R, 2)
# eps_draw<- cbind(rmvevd(N, dep = theta[1], model = "log", d = 2), rgev(N))
eps_draw<- matrix(rgumbel(N*R), N, R)
eps_draw[,2]	<- (1-theta[1])*eps_draw[,1] + sqrt(1-(1-theta[1])^2)*eps_draw[,2]
xbeta	<- do.call(cbind, lapply(1:R, function(i) X_list[[i]] %*% beta + beta0[i]))
psi		<- exp(xbeta + eps_draw)
e_mat	<- all_vec_fn(y, psi, rep(1, N) %*% t(gamma), price, R)
e_sgn	<- 1*(e_mat > 0)

# Estimation # 
sumX_fn	<- function(q, r, theta){
	if(r == 1){
		return(1)
	}else if( r == 2){
		return((1-theta)/theta*q*(q-1)/2 )
	}else{
		A1	<- (1-theta)/theta * seq(q-1, 1, -1)
		A2	<- seq(r-2, 1, -1)
		idx	<- combinations(q-2, r-2, 1:(length(A1)-1))
		out	<- NULL
		for(i in 1:nrow(idx)){
			A3	<- A1[idx[i,]]
			A	<- A3 + A2
			B	<- A1[(idx[i,(r-2)]+1):length(A1)]
			out	<- c(out, prod(A) * sum(B))
		}
		return(sum(out))
	}
}

llfun	<- function(param, y, e_mat, nest, X_list, price){
	N		<- length(y)
	R		<- ncol(e_mat)
	nx		<- ncol(X_list[[1]])
	beta	<- param[1:nx]
	beta0	<- param[(nx+1):(nx+R)]
	gamma	<- exp(param[(nx+R+1):(nx+R*2)])
	theta	<- param[(nx+R*2+1):length(param)]
	theta	<- c(exp(theta)/(1 + exp(theta)), 1)
	
	# nest_idx<- apply(nest, 1, function(x) which(x==1))
	theta1	<- t(theta) %*% t(nest)
	V		<- do.call(cbind, lapply(1:R, function(i) X_list[[i]] %*% beta + beta0[i])) - 
				log(e_mat/(rep(1,N) %*% t(gamma)) + price)
	V		<- V / (rep(1, N) %*% theta1)
	eV		<- exp(V)
	eVs		<- eV %*% nest
	es		<- eVs^(rep(1, N) %*% t(theta))
	ps		<- es/rowSums(es)	
	q		<- e_sgn %*% nest
	A		<- e_mat + price * (rep(1, N) %*% t(gamma))
		
	# Pre-calculate sumX
	sumX	<- vector("list", S)
	q_max	<- colSums(nest)	
	for(i in 1:S){
		sumX[[i]]	<- matrix(0, q_max[i], q_max[i])
		for(j in 1:q_max[i]){
			for(k in 1:j){
				sumX[[i]][j,k]	<- sumX_fn(q=j, r=k, theta[i])
			}
		}
	}
	
	# Calculate matrix b, a N \times (sum_n q_max) matrix
	# b_ij = prod_{n=1}^N ps^{q_n-r_n+1} * (prod_{n=1}^N sum(X_{nr})) * (sum_{n=1}^N (q_n-r_n+1) - 1)!
	r_mat	<- as.matrix(do.call(expand.grid, lapply(q_max, function(x) 0:x)))[-1,]
	b		<- matrix(0, N, nrow(r_mat))
	for(i in 1:nrow(r_mat)){
		b3	<- q - rep(1,N) %*% t(r_mat[i,]) + 1
		b1	<- ps^b3
		# b2	<- sapply(1:S, function(j) ifelse(q[,j]==0 | q[,j]<r_mat[i,j], NA, sumX[[j]][q[,j],r_mat[i,j]] ))
# 		b2		<- sapply(1:S, function(j) 
# 					{out <- rep(1, N); sel <- q[,j]==0 | q[,j]<r_mat[i,j]; out[!sel] <- sumX[[j]][q[!sel,j],r_mat[i,j]]; out})		
		b2		<- sapply(1:S, function(j) 
					{out <- rep(1, N); sel <- q[,j]==0 | q[,j]<r_mat[i,j]|r_mat[i,j]==0; out[!sel] <- sumX[[j]][q[!sel,j],r_mat[i,j]]; out})		
		sel	<- b3 - 1 < 0 | q == 0					# For the observations that have q_n < r_n or q_n = 0
		sel1	<- rowSums(b3 - 1 < 0) > 0			# For the observations that have q_n < r_n, b = 0
		if(any(r_mat[i,] == 0)){					# In addition, for the observations that q_n > 0 and r_n = 0, then b = 0
			sel2	<- q[,r_mat[i,]==0] != 0
			if(sum(sel2) > 0){
				sel1	<- sel1 | sel2
			}
		}
		b1[sel]	<- 1
		b2[sel]	<- 1
		b3[sel]	<- 0
		b[!sel1,i]	<- apply(b1[!sel1,], 1, prod) * apply(b2[!sel1,], 1, prod, na.rm=T) * factorial(rowSums(b3[!sel1,]) - 1)
	}
	ll	<- log(rowSums(A*e_sgn)) - rowSums(log(A)*e_sgn) + rowSums(V*e_sgn) - rowSums(q*log(eVs)) + log(rowSums(b))
	return(ll)
}

param.init	<- c(beta, beta0, log(gamma), log(theta[1]/(1-theta[1])))
system.time(tmp <- llfun(param.init, y = y, e_mat = e_mat, nest = nest, X_list = X_list, price = price))
myfix		<- c(3, 6)		#length(param.init)
sol			<- maxLik(llfun, start=param.init, method="BFGS", fixed = myfix,
					y = y, e_mat = e_mat, nest = nest, X_list = X_list, price = price)
summary(sol)

